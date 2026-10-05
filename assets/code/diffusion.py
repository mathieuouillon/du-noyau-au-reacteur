"""
diffusion.py -- 2-D multigroup neutron diffusion eigenvalue solver.

Solves the steady-state multigroup diffusion equation

    -div( D_g grad phi_g ) + Sigma_r,g phi_g
        = sum_{g'!=g} Sigma_s,g'->g phi_g'  +  (chi_g / k) sum_g' nu_Sigma_f,g' phi_g'

on a structured (possibly non-uniform) Cartesian mesh, using a mesh-centered
finite-volume discretization and power iteration for the k-eigenvalue.

Conventions
-----------
  * Lengths in cm, cross sections in 1/cm, D in cm.
  * Energy group 0 is the FASTEST group (standard reactor-physics ordering).
  * sigma_s[g_from, g_to]; the diagonal (self-scatter) is ignored, since
    self-scatter cancels out of the removal term.
  * 2-D means "per unit height in z"; leakage in z can be folded in with a
    buckling term via Material.buckling_b2.
"""

from dataclasses import dataclass, field
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# --------------------------------------------------------------------------
# Materials
# --------------------------------------------------------------------------
@dataclass
class Material:
    """Macroscopic multigroup cross-section set for one material."""
    name: str
    D: np.ndarray            # (G,)   diffusion coefficient  [cm]
    sigma_a: np.ndarray      # (G,)   absorption             [1/cm]
    sigma_s: np.ndarray      # (G,G)  scatter g_from -> g_to [1/cm]
    nu_sigma_f: np.ndarray   # (G,)   nu * fission           [1/cm]
    chi: np.ndarray          # (G,)   fission spectrum, sums to 1
    buckling_b2: float = 0.0  # transverse (axial) buckking, [1/cm^2]

    def __post_init__(self):
        self.D = np.asarray(self.D, float)
        self.sigma_a = np.asarray(self.sigma_a, float)
        self.sigma_s = np.asarray(self.sigma_s, float)
        self.nu_sigma_f = np.asarray(self.nu_sigma_f, float)
        self.chi = np.asarray(self.chi, float)
        G = self.D.size
        for arr, shape in ((self.sigma_a, (G,)), (self.nu_sigma_f, (G,)),
                           (self.chi, (G,)), (self.sigma_s, (G, G))):
            if arr.shape != shape:
                raise ValueError(f"{self.name}: expected shape {shape}, got {arr.shape}")
        if self.nu_sigma_f.any() and not np.isclose(self.chi.sum(), 1.0):
            raise ValueError(f"{self.name}: chi must sum to 1")

    @property
    def G(self) -> int:
        return self.D.size

    @property
    def sigma_r(self) -> np.ndarray:
        """Removal cross section: absorption + out-scatter (excl. self-scatter),
        plus the transverse-leakage pseudo-absorption D*B^2."""
        out_scatter = self.sigma_s.sum(axis=1) - np.diag(self.sigma_s)
        return self.sigma_a + out_scatter + self.D * self.buckling_b2


# --------------------------------------------------------------------------
# Geometry / mesh
# --------------------------------------------------------------------------
BC_TYPES = ("vacuum", "reflective", "zero")


@dataclass
class Mesh:
    """Structured Cartesian mesh with a per-cell material assignment.

    hx, hy   : cell widths, shape (nx,) and (ny,)
    mat_map  : (ny, nx) integer array indexing into `materials`
    bc       : dict with keys 'xmin','xmax','ymin','ymax'
    """
    hx: np.ndarray
    hy: np.ndarray
    mat_map: np.ndarray
    materials: list
    bc: dict = field(default_factory=lambda: dict(xmin="reflective", xmax="vacuum",
                                                  ymin="reflective", ymax="vacuum"))

    def __post_init__(self):
        self.hx = np.asarray(self.hx, float)
        self.hy = np.asarray(self.hy, float)
        self.mat_map = np.asarray(self.mat_map, int)
        if self.mat_map.shape != (self.ny, self.nx):
            raise ValueError(f"mat_map shape {self.mat_map.shape} != {(self.ny, self.nx)}")
        for side, kind in self.bc.items():
            if kind not in BC_TYPES:
                raise ValueError(f"bad BC '{kind}' on {side}; use one of {BC_TYPES}")
        Gs = {m.G for m in self.materials}
        if len(Gs) != 1:
            raise ValueError("all materials must have the same number of groups")

    @property
    def nx(self): return self.hx.size

    @property
    def ny(self): return self.hy.size

    @property
    def ncells(self): return self.nx * self.ny

    @property
    def G(self): return self.materials[0].G

    @property
    def volume(self):
        """(ny, nx) cell volumes, per unit height."""
        return np.outer(self.hy, self.hx)

    def cell_centers(self):
        xc = np.cumsum(self.hx) - self.hx / 2.0
        yc = np.cumsum(self.hy) - self.hy / 2.0
        return xc, yc


# --------------------------------------------------------------------------
# Discretization
# --------------------------------------------------------------------------
def _interior_coupling(D1, h1, D2, h2):
    """Finite-volume coupling coefficient [1/cm] between two adjacent cells.

    Derived by requiring continuity of flux and current at the interface:
    the net current is  J = coupling * (phi_1 - phi_2) * area.
    """
    return 2.0 * D1 * D2 / (D1 * h2 + D2 * h1)


def _boundary_coupling(D, h, kind):
    """Coupling to the outside world at a domain boundary [1/cm].

    vacuum     : Marshak / linear-extrapolation BC, flux vanishes 2D beyond
                 the surface  ->  J = 2 D phi / (h + 4 D)
    zero       : flux vanishes exactly at the surface -> J = 2 D phi / h
    reflective : no net current
    """
    if kind == "reflective":
        return 0.0
    if kind == "vacuum":
        return 2.0 * D / (h + 4.0 * D)
    if kind == "zero":
        return 2.0 * D / h
    raise ValueError(kind)


def build_operators(mesh):
    """Assemble the loss matrix A and fission-production matrix F.

    Unknowns are ordered group-major: index = g * ncells + (j * nx + i).
    Returns (A, F) as CSR sparse matrices of size (G*ncells, G*ncells).
    """
    nx, ny, G, N = mesh.nx, mesh.ny, mesh.G, mesh.ncells
    hx, hy = mesh.hx, mesh.hy
    V = mesh.volume
    mats = mesh.materials

    # Per-cell, per-group property lookups
    Dc = np.empty((G, ny, nx))
    SRc = np.empty((G, ny, nx))
    for j in range(ny):
        for i in range(nx):
            m = mats[mesh.mat_map[j, i]]
            Dc[:, j, i] = m.D
            SRc[:, j, i] = m.sigma_r

    rows, cols, vals = [], [], []
    frows, fcols, fvals = [], [], []

    def add(r, c, v):
        rows.append(r); cols.append(c); vals.append(v)

    for g in range(G):
        off = g * N
        for j in range(ny):
            for i in range(nx):
                c = j * nx + i
                r = off + c
                diag = SRc[g, j, i] * V[j, i]

                # ---- x-direction faces (area = hy[j]) ----
                if i > 0:
                    cc = _interior_coupling(Dc[g, j, i], hx[i], Dc[g, j, i - 1], hx[i - 1]) * hy[j]
                    diag += cc
                    add(r, off + c - 1, -cc)
                else:
                    diag += _boundary_coupling(Dc[g, j, i], hx[i], mesh.bc["xmin"]) * hy[j]

                if i < nx - 1:
                    cc = _interior_coupling(Dc[g, j, i], hx[i], Dc[g, j, i + 1], hx[i + 1]) * hy[j]
                    diag += cc
                    add(r, off + c + 1, -cc)
                else:
                    diag += _boundary_coupling(Dc[g, j, i], hx[i], mesh.bc["xmax"]) * hy[j]

                # ---- y-direction faces (area = hx[i]) ----
                if j > 0:
                    cc = _interior_coupling(Dc[g, j, i], hy[j], Dc[g, j - 1, i], hy[j - 1]) * hx[i]
                    diag += cc
                    add(r, off + c - nx, -cc)
                else:
                    diag += _boundary_coupling(Dc[g, j, i], hy[j], mesh.bc["ymin"]) * hx[i]

                if j < ny - 1:
                    cc = _interior_coupling(Dc[g, j, i], hy[j], Dc[g, j + 1, i], hy[j + 1]) * hx[i]
                    diag += cc
                    add(r, off + c + nx, -cc)
                else:
                    diag += _boundary_coupling(Dc[g, j, i], hy[j], mesh.bc["ymax"]) * hx[i]

                add(r, r, diag)

                # ---- in-scattering (moved to the LHS, so it is negative) ----
                m = mats[mesh.mat_map[j, i]]
                for gp in range(G):
                    if gp == g:
                        continue
                    s = m.sigma_s[gp, g]
                    if s != 0.0:
                        add(r, gp * N + c, -s * V[j, i])

                # ---- fission production ----
                if m.chi[g] != 0.0:
                    for gp in range(G):
                        if m.nu_sigma_f[gp] != 0.0:
                            frows.append(r)
                            fcols.append(gp * N + c)
                            fvals.append(m.chi[g] * m.nu_sigma_f[gp] * V[j, i])

    n = G * N
    A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    F = sp.csr_matrix((fvals, (frows, fcols)), shape=(n, n))
    return A, F


# --------------------------------------------------------------------------
# Eigenvalue solve
# --------------------------------------------------------------------------
def _source_error(s_new, s_old):
    """Max relative change of the normalized fission source, over cells that
    actually carry a source (reflector cells have none and must be skipped)."""
    active = s_old > 1e-12 * s_old.max()
    if not active.any():
        return 0.0
    return float(np.max(np.abs(s_new[active] - s_old[active]) / s_old[active]))


@dataclass
class Solution:
    k_eff: float
    flux: np.ndarray          # (G, ny, nx)
    fission_source: np.ndarray  # (ny, nx)
    iterations: int
    converged: bool
    residual: float
    history: list


def solve(mesh, tol_k=1e-10, tol_s=1e-8, max_iter=2000,
          shift=0.05, n_warmup=15, verbose=False):
    """Power iteration for the fundamental-mode k-eigenvalue of

        A phi = (1/k) F phi

    Scattering lives inside A rather than in a group-sweep loop, so
    up-scattering is handled with no extra machinery.

    Plain power iteration converges at the dominance ratio k1/k0, which for a
    large loosely-coupled core is often 0.95-0.99 -- hundreds of iterations.
    The Wielandt shift fixes this: iterate instead on

        (A - F/ke) phi = lambda F phi,     lambda = 1/k - 1/ke

    with the shift ke placed just above k. That pushes the shifted dominance
    ratio far below the unshifted one, at the cost of refactorizing whenever
    the shift is retargeted. Set shift=None for textbook power iteration.
    """
    nx, ny, G, N = mesh.nx, mesh.ny, mesh.G, mesh.ncells
    A, F = build_operators(mesh)
    Acsc = A.tocsc()
    Fcsc = F.tocsc()

    phi = np.ones(G * N)
    src = F @ phi
    if src.sum() <= 0:
        raise ValueError("no fissile material in the problem")
    src = src / src.sum()
    k = 1.0
    history = []
    converged = False
    n_factor = 0

    if shift is None:
        # ---------------- textbook power iteration ----------------
        # A^-1 F has eigenvalue k, so with a unit-sum source the growth
        # factor of the source IS k -- no multiplicative bookkeeping needed.
        lu = spla.splu(Acsc)
        n_factor = 1
        for it in range(1, max_iter + 1):
            phi = lu.solve(src)
            src_raw = F @ phi
            k_new = src_raw.sum()
            src_new = src_raw / k_new

            k_err = abs(k_new - k) / k_new
            s_err = _source_error(src_new, src)
            history.append((it, k_new, k_err, s_err))
            if verbose:
                print(f"  it {it:4d}  k = {k_new:.8f}  dk = {k_err:.2e}  dS = {s_err:.2e}")

            src, k = src_new, k_new
            if k_err < tol_k and s_err < tol_s:
                converged = True
                break
    else:
        # ---------------- Wielandt-shifted iteration ----------------
        # (A - F/ke)^-1 F has eigenvalue 1/lambda with lambda = 1/k - 1/ke.
        # The fundamental mode has the k nearest below ke, hence the smallest
        # lambda, hence the largest 1/lambda -- still the dominant mode, so
        # ordinary power iteration on the shifted operator finds it.
        #
        # DANGER: this only holds if ke sits ABOVE the true k. If the shift is
        # applied while the estimate is still far too low, modes above ke get
        # negative lambda, the iteration can lock onto the wrong mode, and it
        # will then "converge" quietly to a wrong answer. So: take a few plain
        # power iterations first to get a usable estimate, and verify the
        # eigenpair against the unshifted residual at the end regardless.
        lu0 = spla.splu(Acsc)
        n_factor = 1
        for it in range(1, n_warmup + 1):
            phi = lu0.solve(src)
            src_raw = F @ phi
            k = src_raw.sum()
            src = src_raw / k
            history.append((it, k, np.nan, np.nan))
        del lu0

        ke = k + shift
        lu = spla.splu((Acsc - Fcsc / ke).tocsc())
        n_factor += 1

        for it in range(n_warmup + 1, max_iter + 1):
            phi = lu.solve(src)
            src_raw = F @ phi
            rho = src_raw.sum()              # -> 1 / lambda
            k_new = 1.0 / (1.0 / rho + 1.0 / ke)
            src_new = src_raw / rho

            k_err = abs(k_new - k) / k_new
            s_err = _source_error(src_new, src)
            history.append((it, k_new, k_err, s_err))
            if verbose:
                print(f"  it {it:4d}  k = {k_new:.8f}  dk = {k_err:.2e}  "
                      f"dS = {s_err:.2e}  ke = {ke:.4f}")

            src, k = src_new, k_new
            if k_err < tol_k and s_err < tol_s:
                converged = True
                break

            # retarget the shift once it has drifted from the eigenvalue
            if abs((k + shift) - ke) > 0.25 * shift:
                ke = k + shift
                lu = spla.splu((Acsc - Fcsc / ke).tocsc())
                n_factor += 1

    # ---- verify the eigenpair against the ORIGINAL, unshifted problem ----
    # This is independent of whatever iteration scheme produced it, so it
    # catches wrong-mode convergence, a bad shift, or an assembly bug.
    rhs = (F @ phi) / k
    residual = np.linalg.norm(A @ phi - rhs) / np.linalg.norm(rhs)
    positive = bool(np.all(phi > 0))   # the fundamental mode has no sign change

    if verbose:
        print(f"  {n_factor} LU factorization(s), {len(history)} iterations")
        print(f"  residual = {residual:.3e}   strictly positive flux: {positive}")

    if residual > 1e-6 or not positive:
        raise RuntimeError(
            f"solve() did not find the fundamental mode "
            f"(k={k:.6f}, residual={residual:.2e}, positive={positive}). "
            f"Try shift=None or a larger shift.")

    # normalize so the mean flux is 1
    phi = phi / phi.mean()
    flux = phi.reshape(G, ny, nx)
    fsrc = (F @ phi)[:N].reshape(ny, nx)  # group-0 row carries chi_0 * total production
    return Solution(k, flux, fsrc, len(history), converged, residual, history)
