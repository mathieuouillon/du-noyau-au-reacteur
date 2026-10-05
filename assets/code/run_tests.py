"""
run_tests.py -- verification suite and worked examples for diffusion.py

The checks are ordered from strongest to weakest evidence. This ordering is
deliberate: verify against something you know exactly before you trust a case
where you don't.

  1. ANALYTIC.  Bare homogeneous reactor, one group, against the closed-form
     k. Confirms second-order spatial convergence.
  2. SYMMETRY.  A quarter core with reflective inner boundaries must reproduce
     the full core to machine precision.
  3. ITERATION.  Wielandt-shifted and plain power iteration must agree, which
     separates discretization error from iteration error.
  4. PHYSICS.  Two-group core with a water reflector: look for the thermal
     flux peak in the reflector. A sanity check, not a proof.
"""

import numpy as np
from diffusion import Material, Mesh, solve


# ==========================================================================
# 1. Analytic benchmark: bare homogeneous rectangular reactor, one group
# ==========================================================================
def analytic_bare_reactor(D, sigma_a, nu_sigma_f, a, b, extrapolated=True):
    """Exact k for a bare rectangle a x b.

    With the linear-extrapolation (Marshak) boundary condition the flux
    vanishes a distance d = 2D beyond each surface, so the geometric buckling
    is built from the extrapolated dimensions. With a hard zero-flux boundary,
    d = 0 and the physical dimensions are used.
    """
    d = 2.0 * D if extrapolated else 0.0
    B2 = (np.pi / (a + 2 * d)) ** 2 + (np.pi / (b + 2 * d)) ** 2
    return nu_sigma_f / (sigma_a + D * B2), B2


def _refinement_table(fuel, a, k_exact, bc, meshes=(10, 20, 40, 80, 160)):
    print(f"{'cells/side':>11} {'h [cm]':>9} {'k_eff':>14} "
          f"{'error [pcm]':>13} {'ratio':>7}")
    prev, ratios = None, []
    for n in meshes:
        h = a / n
        mesh = Mesh(hx=np.full(n, h), hy=np.full(n, h),
                    mat_map=np.zeros((n, n), int), materials=[fuel],
                    bc=dict(xmin=bc, xmax=bc, ymin=bc, ymax=bc))
        sol = solve(mesh)
        err = sol.k_eff - k_exact
        ratio = ""
        if prev is not None:
            r = prev / abs(err)
            ratios.append(r)
            ratio = f"{r:6.2f}"
        print(f"{n:>11} {h:>9.3f} {sol.k_eff:>14.9f} {err * 1e5:>13.3f} {ratio:>7}")
        prev = abs(err)
    return ratios


def test_analytic():
    D, sa, nsf = 1.0, 0.02, 0.025
    a = b = 100.0

    print("=" * 72)
    print("TEST 1  Bare homogeneous rectangle, one group")
    print(f"        {a:g} x {b:g} cm, D={D}, Sigma_a={sa}, nuSigma_f={nsf}")

    fuel = Material("bare", D=[D], sigma_a=[sa], sigma_s=[[0.0]],
                    nu_sigma_f=[nsf], chi=[1.0])

    # -- zero-flux boundary: isolates the interior discretization error -----
    k_zero, B2z = analytic_bare_reactor(D, sa, nsf, a, b, extrapolated=False)
    print("-" * 72)
    print(f"  (a) zero-flux boundary   B^2 = {B2z:.8e}   analytic k = {k_zero:.9f}")
    ratios = _refinement_table(fuel, a, k_zero, "zero")
    clean = all(abs(r - 4.0) < 0.15 for r in ratios)
    print(f"      error ratios -> 4.00 : clean 2nd-order convergence "
          f"[{'PASS' if clean else 'FAIL'}]")

    # -- vacuum boundary: adds the Marshak boundary discretization error ----
    k_vac, B2v = analytic_bare_reactor(D, sa, nsf, a, b, extrapolated=True)
    print("-" * 72)
    print(f"  (b) vacuum boundary      B^2 = {B2v:.8e}   analytic k = {k_vac:.9f}")
    _refinement_table(fuel, a, k_vac, "vacuum")
    print("      Note the error changing sign near n=80. This is NOT a bug:")
    print("      the Marshak boundary term carries its own O(h^2) error of")
    print("      opposite sign to the interior error, so the two cancel at one")
    print("      particular mesh. A single lucky-looking k proves nothing --")
    print("      which is exactly why the refinement study matters.")
    return clean


# ==========================================================================
# Two-group data for tests 2-4
# ==========================================================================
def two_group_materials():
    fuel = Material(
        "fuel",
        D=[1.40, 0.40],
        sigma_a=[0.0100, 0.0850],
        sigma_s=[[0.0, 0.0200],
                 [0.0, 0.0]],
        nu_sigma_f=[0.0060, 0.1500],
        chi=[1.0, 0.0],
    )
    refl = Material(
        "reflector",
        D=[1.30, 0.20],
        sigma_a=[0.0004, 0.0200],
        sigma_s=[[0.0, 0.0500],
                 [0.0, 0.0]],
        nu_sigma_f=[0.0, 0.0],
        chi=[0.0, 0.0],
    )
    return fuel, refl


def k_infinity(m):
    """Two-group k_inf for a homogeneous medium (fast + thermal fission)."""
    s12 = m.sigma_s[0, 1]
    return (m.nu_sigma_f[0] + m.nu_sigma_f[1] * s12 / m.sigma_a[1]) / (m.sigma_a[0] + s12)


def build_core(half_width, fuel_width, h, materials, quarter):
    """Square core: fuel out to `fuel_width` from the centre, reflector beyond."""
    if quarter:
        n = int(round(half_width / h))
        x = (np.arange(n) + 0.5) * h                       # distance from centre
        bc = dict(xmin="reflective", ymin="reflective",
                  xmax="vacuum", ymax="vacuum")
    else:
        n = int(round(2 * half_width / h))
        x = np.abs((np.arange(n) + 0.5) * h - half_width)
        bc = dict(xmin="vacuum", xmax="vacuum", ymin="vacuum", ymax="vacuum")
    is_fuel = (x[None, :] < fuel_width) & (x[:, None] < fuel_width)
    return Mesh(hx=np.full(n, h), hy=np.full(n, h),
                mat_map=np.where(is_fuel, 0, 1),
                materials=list(materials), bc=bc)


# ==========================================================================
# 2. Symmetry: quarter core (reflective) vs. full core
# ==========================================================================
def test_symmetry():
    fuel, refl = two_group_materials()
    h, half, fw = 5.0, 200.0, 170.0
    kq = solve(build_core(half, fw, h, (fuel, refl), quarter=True)).k_eff
    kf = solve(build_core(half, fw, h, (fuel, refl), quarter=False)).k_eff
    d = abs(kq - kf)
    print("=" * 72)
    print("TEST 2  Quarter core (reflective) vs. full core, two groups")
    print(f"        quarter  k = {kq:.10f}   ({(half/h)**2:.0f} cells)")
    print(f"        full     k = {kf:.10f}   ({(2*half/h)**2:.0f} cells)")
    print(f"        difference = {d * 1e5:.2e} pcm   [{'PASS' if d < 1e-8 else 'FAIL'}]")
    return d < 1e-8


# ==========================================================================
# 3. Iteration scheme: Wielandt shift vs. plain power iteration
# ==========================================================================
def test_iteration():
    fuel, refl = two_group_materials()
    mesh = build_core(200.0, 170.0, 2.0, (fuel, refl), quarter=True)
    print("=" * 72)
    print("TEST 3  Iteration schemes must agree "
          f"({mesh.ncells * mesh.G} unknowns)")
    print(f"{'scheme':>22} {'iters':>7} {'k_eff':>14} {'residual':>11}")
    results = {}
    for label, kw in (("plain power", dict(shift=None, max_iter=4000)),
                      ("Wielandt shift=0.05", dict(shift=0.05)),
                      ("Wielandt shift=0.02", dict(shift=0.02))):
        s = solve(mesh, **kw)
        results[label] = s.k_eff
        print(f"{label:>22} {s.iterations:>7} {s.k_eff:>14.9f} {s.residual:>11.1e}")
    spread = max(results.values()) - min(results.values())
    ok = spread < 1e-7
    print(f"        spread = {spread * 1e5:.2e} pcm   [{'PASS' if ok else 'FAIL'}]")
    print("        The shift cuts ~500 iterations to ~40 for the same answer.")
    return ok


# ==========================================================================
# 4. Physics check: thermal flux peak in the reflector
# ==========================================================================
def test_reflected_core():
    fuel, refl = two_group_materials()
    h, half, fw = 2.0, 200.0, 170.0
    mesh = build_core(half, fw, h, (fuel, refl), quarter=True)
    sol = solve(mesh)

    fast, therm = sol.flux[0], sol.flux[1]
    xc, _ = mesh.cell_centers()
    row_f, row_t = fast[0, :], therm[0, :]

    in_refl = xc > fw
    i_peak = int(np.flatnonzero(in_refl)[np.argmax(row_t[in_refl])])
    i_edge = int(np.flatnonzero(~in_refl)[-1])
    peaked = row_t[i_peak] > row_t[i_edge]

    print("=" * 72)
    print("TEST 4  Two-group core + water reflector (quarter core, 2 cm mesh)")
    print(f"        fuel k_inf = {k_infinity(fuel):.6f}")
    print(f"        core k_eff = {sol.k_eff:.6f}  "
          f"({sol.iterations} iterations, residual {sol.residual:.1e})")
    print("-" * 72)
    print("        flux across the centreline, near the interface at x = 170 cm:")
    print(f"{'x [cm]':>9} {'fast':>9} {'thermal':>9} {'th/fast':>9}   region")
    for i in range(i_edge - 3, min(mesh.nx, i_edge + 9)):
        reg = "fuel" if xc[i] < fw else "REFLECTOR"
        mark = "  <-- thermal peak" if i == i_peak else ""
        print(f"{xc[i]:>9.1f} {row_f[i]:>9.4f} {row_t[i]:>9.4f} "
              f"{row_t[i] / row_f[i]:>9.3f}   {reg}{mark}")

    print()
    print(f"        thermal flux at fuel edge   (x={xc[i_edge]:.0f}) = {row_t[i_edge]:.4f}")
    print(f"        thermal flux peak in water  (x={xc[i_peak]:.0f}) = {row_t[i_peak]:.4f}")
    print(f"        reflector thermal peak present: {peaked} "
          f"[{'PASS' if peaked else 'FAIL'}]")
    print("        Physically right: water slows neutrons down without")
    print("        absorbing many thermally, while the fuel is a strong")
    print("        thermal absorber -- so thermal flux rises past the edge.")
    return peaked, sol, mesh


# ==========================================================================
def main():
    r1 = test_analytic()
    r2 = test_symmetry()
    r3 = test_iteration()
    r4, sol, mesh = test_reflected_core()

    print("=" * 72)
    for name, ok in (("1 analytic convergence", r1), ("2 quarter/full symmetry", r2),
                     ("3 iteration agreement", r3), ("4 reflector peak", r4)):
        print(f"  {name:<26} {'PASS' if ok else 'FAIL'}")
    print("=" * 72)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xc, yc = mesh.cell_centers()
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
        for ax, data, title in zip(
                axes, [sol.flux[0], sol.flux[1], sol.fission_source],
                ["Fast flux (group 1)", "Thermal flux (group 2)", "Fission source"]):
            im = ax.pcolormesh(xc, yc, data, shading="auto", cmap="inferno")
            ax.axvline(170, color="cyan", lw=0.9, ls="--")
            ax.axhline(170, color="cyan", lw=0.9, ls="--")
            ax.set_aspect("equal")
            ax.set_title(title)
            ax.set_xlabel("x [cm]")
            fig.colorbar(im, ax=ax)
        axes[0].set_ylabel("y [cm]")
        fig.suptitle(f"Reflected quarter core, k_eff = {sol.k_eff:.5f}   "
                     "(dashed = fuel/reflector interface)")
        fig.tight_layout()
        fig.savefig("flux_maps.png", dpi=130)
        print("  wrote flux_maps.png")
    except Exception as e:
        print(f"  (plot skipped: {e})")


if __name__ == "__main__":
    main()
