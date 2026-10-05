---
title: "Le solveur de diffusion multigroupe (en anglais)"
date: 2026-10-02 09:00:00 +0200
categories: [Neutronique, Méthodes numériques]
tags: [diffusion, valeur propre, Wielandt, vérification]
description: "Physique, discrétisation et itération d'un solveur de diffusion 2D vérifié contre une solution analytique."
image:
  path: /assets/img/nucleaire/flux_maps.png
  alt: "Le solveur de diffusion multigroupe (en anglais)"
---

> Cet article a été écrit en anglais, au tout début du projet. Il détaille le solveur sur lequel repose [Neutronique du cœur]({{ '/posts/neutronique-du-coeur/' | relative_url }}).
{: .prompt-info }

A walkthrough of [`diffusion.py`]({{ '/assets/code/diffusion.py' | relative_url }}), from the neutron balance to the last line of
`solve()`. Read it alongside the code.


![Flux rapide, flux thermique et source de fission du cœur réfléchi]({{ '/assets/img/nucleaire/flux_maps.png' | relative_url }})
_Flux rapide, flux thermique et source de fission du cœur réfléchi._

---

## Part 1 — The physics

### 1.1 What the flux actually is

Everything is built on the **scalar flux** `phi(r, E)`. The name is misleading:
it is not a flow of anything, and it has no direction. It is

```
phi = n * v          [neutrons/cm^2/s]
```

the neutron density times their speed. The reason it is the useful variable is
that reaction rates are linear in it:

```
reaction rate density = Sigma_x * phi     [reactions/cm^3/s]
```

where `Sigma_x` is a **macroscopic cross section** in 1/cm — the probability of
reaction type `x` per unit path length. A better mental model for `phi` is
*total neutron track length per unit volume per unit time*. Reaction rate =
track length × probability per length. That framing makes the units obvious and
makes Monte Carlo track-length estimators intuitive later.

Consequence worth internalizing: **flux is not power**. Power follows the
fission rate `Sigma_f * phi`. In a core with a water reflector the thermal flux
peaks *outside* the fuel, where the power is exactly zero. You will see this in
Test 4.

### 1.2 Where diffusion comes from

The honest equation is the **neutron transport equation**, which tracks
`psi(r, E, Omega)` — flux resolved by direction `Omega`, a 6-D problem. Nobody
wants to solve that if they can avoid it.

Diffusion theory is what you get by expanding the angular dependence in
spherical harmonics and truncating after the linear term (the "P1"
approximation). The result is **Fick's law**:

```
J = -D grad(phi)
```

Net current flows down the flux gradient. `D` is the diffusion coefficient in
cm, roughly `1/(3 Sigma_tr)`.

This is an *approximation*, and knowing when it breaks is most of the skill:

| Diffusion is fine | Diffusion fails |
|---|---|
| Deep inside a large homogeneous region | Within ~2 mean free paths of a boundary |
| Weak absorption (scattering dominates) | Strong absorbers: control rods, burnable poison |
| Smooth, near-isotropic angular flux | Voids and streaming paths (D → infinity is nonsense) |
| Away from sharp material interfaces | Strongly anisotropic flux near interfaces |

This is why real core simulators do not use raw finite-difference diffusion —
they use nodal diffusion with *discontinuity factors* calibrated by transport
calculations, so the diffusion solution reproduces transport answers at the
assembly scale. Diffusion is a cheap engine that has been taught to mimic a
correct one.

### 1.3 The balance equation

Write the neutron balance for one energy group `g`, in steady state. Every term
is neutrons/cm³/s:

```
-div(D_g grad phi_g)  +  Sigma_r,g phi_g  =  sum_{g'!=g} Sigma_s,g'->g phi_g'  +  chi_g/k * sum_g' nu_Sigma_f,g' phi_g'
\__________________/     \______________/     \_________________________/     \___________________________________/
     leakage                 removal                  in-scattering                      fission
```

- **Leakage** — net flow out of a volume, from Fick's law.
- **Removal** — everything that takes a neutron out of group `g` *here*:
  absorption, plus scattering to any other group.
- **In-scattering** — neutrons arriving from other groups.
- **Fission** — neutrons born from fissions in all groups, distributed over
  energy by the fission spectrum `chi_g`.

#### The self-scatter subtlety

Why is the removal term `Sigma_a + Sigma_s,out` and not the total cross section?
Start from the honest form, with total cross section on the left and *all*
in-scattering including `g -> g` on the right:

```
... + Sigma_t,g phi_g = ... + Sigma_s,g->g phi_g + sum_{g'!=g} Sigma_s,g'->g phi_g' + fission
```

The `Sigma_s,g->g phi_g` term appears on both sides — a neutron that scatters
within its own group has not changed group, so it is not a loss. Move it left:

```
(Sigma_t,g - Sigma_s,g->g) phi_g = Sigma_a,g + sum_{g'!=g} Sigma_s,g->g'
```

which is exactly `Material.sigma_r`. This is why the code ignores the diagonal
of `sigma_s`. It genuinely does not matter what you put there.

```python
@property
def sigma_r(self):
    out_scatter = self.sigma_s.sum(axis=1) - np.diag(self.sigma_s)
    return self.sigma_a + out_scatter + self.D * self.buckling_b2
```

The third term is a trick: `D*B^2` is a *pseudo-absorption* representing leakage
in the un-modeled z direction. If the axial flux is a cosine of buckling `B^2`,
the axial leakage is exactly `D B^2 phi`, so a 2-D solve with this term
approximates a 3-D core. Cheap and surprisingly effective.

### 1.4 Multigroup: discretizing energy

Neutrons are born fast (~2 MeV) and, in a thermal reactor, must slow to ~0.025 eV
to fission efficiently. Cross sections vary over that range by orders of
magnitude, with resonances.

The multigroup method chops energy into bins and uses one flux-weighted average
cross section per bin. Two groups suffice for LWR demonstrations:

- **Group 0 — fast.** Where neutrons are born. `chi = [1, 0]` says essentially
  all fission neutrons appear here.
- **Group 1 — thermal.** Where most fission happens: `nu_Sigma_f = [0.006, 0.150]`,
  so the thermal group is 25× more productive.

The only coupling in the example data is `sigma_s[0,1] = 0.020` — slowing down,
fast to thermal. Nothing scatters back up (`sigma_s[1,0] = 0`), which is the
usual approximation above thermal energies.

The catch, and it is a big one: **flux-weighted averages require the flux you
haven't computed yet.** Generating good group constants is a whole discipline
(lattice physics: resonance self-shielding, spatial homogenization). It is why
the README warns that the data pipeline usually outweighs the solver.

### 1.5 Why it is an eigenvalue problem

Look at the balance equation again. There is no external source — fission is
proportional to the flux itself. So the equation is homogeneous:

```
A phi = (1/k) F phi
```

A homogeneous system has only the trivial solution `phi = 0` *unless* the
operator is singular in just the right way. Generically no `phi` balances
exactly: a reactor is either subcritical (dying) or supercritical (growing).

`k` is the fudge factor that forces balance. Divide the fission source by `k`
and there is exactly one value that admits a steady, everywhere-positive
solution. Physically:

```
k = neutron production rate / neutron loss rate
```

- `k = 1` — critical, self-sustaining.
- `k > 1` — supercritical.
- `k < 1` — subcritical.

`k` is reported in **pcm** (per cent mille, 1e-5) because reactor design
cares about differences of a few hundred pcm.

So the code is not solving a linear system. It is finding the **fundamental
eigenpair** — the largest `k` and its associated flux shape. Everything in
`solve()` follows from that.

---

## Part 2 — Discretization

### 2.1 Finite volume: integrate, don't differentiate

Integrate the balance equation over cell `(i,j)` and apply the divergence
theorem to the leakage term:

```
integral of -div(D grad phi) dV  =  -integral of D grad(phi) . n dS  =  sum over 4 faces of J_face * A_face
```

The volume integral of a divergence became a **sum of surface currents**. This
is the whole idea, and it is why finite volume is the right choice here:

**It is conservative by construction.** Whatever leaves cell `i` through the
right face enters cell `i+1` through its left face — the same number, with
opposite sign, because it is literally the same matrix entry. Neutrons cannot
be created or destroyed by discretization error. With a finite-difference
scheme written directly on the derivatives you have to work for that.

Every remaining term is assumed constant over the cell, giving `Sigma * phi * V`.

### 2.2 The interface coupling coefficient

This is the most important derivation in the code. Two adjacent cells with
different materials — how much current crosses the interface?

Let `phi_s` be the (unknown) flux at the interface. Assume linear flux in each
half-cell. Current from the left, using Fick's law over distance `h_1/2`:

```
J = 2 D_1 (phi_1 - phi_s) / h_1        =  a (phi_1 - phi_s),    a = 2 D_1 / h_1
```

And from the right:

```
J = 2 D_2 (phi_s - phi_2) / h_2        =  b (phi_s - phi_2),    b = 2 D_2 / h_2
```

**Current must be continuous** across the interface — neutrons don't pile up on
a surface. Setting these equal:

```
a phi_1 - a phi_s = b phi_s - b phi_2
phi_s = (a phi_1 + b phi_2) / (a + b)
```

Substitute back:

```
J = a (phi_1 - phi_s) = (a b / (a + b)) (phi_1 - phi_2)
```

So the coupling is the **harmonic mean** `ab/(a+b)`. Expanding `a` and `b`:

```python
def _interior_coupling(D1, h1, D2, h2):
    return 2.0 * D1 * D2 / (D1 * h2 + D2 * h1)
```

Two things to notice:

1. **Harmonic, not arithmetic.** Averaging `D` arithmetically across a
   fuel/reflector interface would be flatly wrong. The harmonic mean is
   dominated by the *smaller* `D` — the more diffusion-resistant material
   controls the flow, exactly like resistors in series. This is why the code
   handles material discontinuities correctly with no special-case logic.
2. **The interface flux `phi_s` was eliminated.** It never enters the matrix.
   It was a scaffold used to enforce continuity and then discarded.

### 2.3 Boundary conditions

At a domain boundary there is no neighbor, so `phi_s` is set by physics instead.

**Vacuum (Marshak).** No neutrons come back in from outside. In diffusion theory
the partial inward current is `J- = phi/4 + (D/2) dphi/dx`. Setting it to zero:

```
dphi/dx = -phi / (2D)
```

A flux with that slope hits zero a distance `d = 2D` beyond the surface. So
extrapolate linearly to zero at `h/2 + 2D` from the cell center:

```
J = D phi / (h/2 + 2D) = 2 D phi / (h + 4 D)
```

(Exact transport theory gives `d = 0.7104 * lambda_tr ≈ 2.13 D`, so even this
"exact" BC is a diffusion-theory approximation. It is a boundary, and Part 1.2
warned that diffusion is worst at boundaries.)

**Zero flux:** the same with `d = 0`, giving `2D/h`. Not physical for a real
surface, but it matches the textbook analytic solution exactly, which is why
Test 1 uses it — see Part 4.

**Reflective:** zero net current, so the coupling is `0`. The cell simply has no
term for that face. This is how the quarter-core symmetry in Test 2 works: a
reflective plane is indistinguishable from an identical core mirrored across it.

### 2.4 Assembling the matrix

`build_operators` walks every (cell, group) pair and writes one matrix row.
Unknowns are ordered **group-major**: `index = g * ncells + (j * nx + i)`.

For row `(c, g)` the diagonal accumulates:

```python
diag = SRc[g, j, i] * V[j, i]        # removal
# ... then for each of 4 faces:
if i > 0:
    cc = _interior_coupling(...) * hy[j]   # coupling * face area
    diag += cc                              # leaving this cell
    add(r, off + c - 1, -cc)                # arriving from the neighbor
else:
    diag += _boundary_coupling(...) * hy[j] # leaks out, never comes back
```

The `+cc` on the diagonal and `-cc` off-diagonal are the same number. That
pairing *is* neutron conservation, expressed in linear algebra.

In-scattering goes on the **left** with a minus sign, because everything except
fission belongs in `A`:

```python
add(r, gp * N + c, -s * V[j, i])
```

Note `gp * N + c` — same cell `c`, different group. Scattering couples groups at
a point; leakage couples cells within a group. **`A` is not block diagonal**,
and that is deliberate. Many textbook codes solve group-by-group with an outer
"scattering iteration" loop; by putting scattering inside `A` and factorizing
the whole thing, up-scattering (thermal neutrons gaining energy, essential
below ~1 eV) needs no extra machinery at all. Put entries in the lower triangle
of `sigma_s` and it just works.

Fission goes into a separate matrix `F`:

```python
fvals.append(m.chi[g] * m.nu_sigma_f[gp] * V[j, i])
```

Structurally `F` is, per cell, the outer product `chi ⊗ nu_Sigma_f` — a **rank-1
matrix**. Every fission neutron forgets the energy of the neutron that caused
it; only the total production rate survives, redistributed by `chi`. That rank-1
structure is why a single scalar `k` can characterize the whole reactor.

---

## Part 3 — Solving the eigenproblem

### 3.1 Power iteration

We want the largest `k` in `A phi = (1/k) F phi`. Rearrange:

```
A^-1 F phi = k phi
```

So `k` is an ordinary eigenvalue of `M = A^-1 F`, and the classic algorithm
applies: repeatedly apply `M` and renormalize. Any starting vector is a mix of
eigenvectors; each application multiplies component `j` by `k_j`, so the largest
one takes over.

The code never forms `M`. It factorizes `A` once with a sparse LU and reuses it:

```python
lu = spla.splu(Acsc)
for it in range(1, max_iter + 1):
    phi = lu.solve(src)      # apply A^-1
    src_raw = F @ phi        # apply F
    k_new = src_raw.sum()    # the growth factor IS k
    src = src_raw / k_new
```

Why is `k` just `src_raw.sum()`? Because `src` is normalized to unit sum each
pass, so the growth of the fission source over one application of `M` is the
eigenvalue directly. No multiplicative bookkeeping — and getting that wrong was
one of the two real bugs in this build.

Physically, this is not a numerical trick. It is a **simulated neutron
generation**: start with a guessed fission source, transport those neutrons to
where they get absorbed (`A^-1`), see how many new fissions result (`F`). The
ratio is neutrons this generation over neutrons last generation, which is the
definition of `k`. Power iteration is watching a reactor settle into its natural
flux shape.

#### The dominance ratio

Convergence is geometric at the **dominance ratio** `k_1/k_0`. For the example
core, computed directly:

```
leading k eigenvalues: [1.365785, 1.324577, 1.324577, 1.285605]
dominance ratio = 0.9698
```

That `k_1` appearing **twice** is not a numerical accident. `k_0` is the smooth
fundamental mode peaked at the core center; the first harmonics are the modes
with one sign change, tilted along x or along y. The square core cannot tell
those two apart, so they are exactly degenerate. Geometric symmetry becomes
eigenvalue degeneracy.

A ratio of 0.97 means each iteration only kills 3% of the error — hundreds of
iterations. Large power reactors are worse (loosely coupled: the two ends barely
talk to each other), and that is precisely when you need the answer most.

### 3.2 The Wielandt shift

Standard fix. Move part of the fission term to the left:

```
(A - F/ke) phi = lambda F phi,      lambda = 1/k - 1/ke
```

Same eigenvectors, same `k` values — only the *spacing* changes. The eigenvalues
of the shifted operator are

```
1/lambda_j = k_j * ke / (ke - k_j)
```

As `ke` approaches `k_0` from above, the denominator for `j=0` goes to zero and
that mode blows up relative to everything else. The new dominance ratio is

```
[k_1 (ke - k_0)] / [k_0 (ke - k_1)]
```

Plugging in the real numbers:

| shift | shifted DR | predicted iterations |
|---|---|---|
| none | 0.9698 | ~750 |
| 0.05 | 0.5317 | ~36 |
| 0.02 | 0.3169 | ~20 |

Observed: 507, 41, 31. The right order, and the mechanism is confirmed. (Plain
iteration beats its estimate because the stopping test measures the *change* in
`k` between iterations, not the true error.)

The cost: `A - F/ke` must be refactorized whenever `ke` moves. That is why the
code only retargets after the shift has drifted:

```python
if abs((k + shift) - ke) > 0.25 * shift:
    ke = k + shift
    lu = spla.splu((Acsc - Fcsc / ke).tocsc())
```

Trading a handful of factorizations for hundreds of avoided iterations.

### 3.3 The trap that cost me an hour

Everything above assumes **`ke` sits above the true `k_0`**. If it doesn't,
modes with `k_j > ke` get *negative* `lambda`, the ordering of `|1/lambda|`
scrambles, and power iteration happily converges to the wrong mode — then
reports success, because the iterates genuinely stop changing.

That is exactly what happened. With `shift=0.02` starting from `k=1.0`, the
initial `ke = 1.02` sat far below the true 1.366, and the solver returned
**1.3248 instead of 1.3658** with every convergence flag green. Note that 1.3248
is suspiciously close to `k_1 = 1.324577` — it had locked onto a first harmonic.

Two defenses, both in the code:

**1. Warm up unshifted.** Take 15 plain power iterations first, so `k` is a
decent estimate before any shift is applied.

**2. Verify the eigenpair independently.** This is the one that matters:

```python
rhs = (F @ phi) / k
residual = np.linalg.norm(A @ phi - rhs) / np.linalg.norm(rhs)
positive = bool(np.all(phi > 0))
if residual > 1e-6 or not positive:
    raise RuntimeError(...)
```

The residual is checked against the **original, unshifted** problem. It does not
care how the answer was produced, so it catches wrong-mode convergence, a bad
shift, a broken iteration update, or an assembly bug equally well.

The positivity check is elegant physics-as-diagnostic. By Perron–Frobenius, the
fundamental mode is the **only** eigenvector that is positive everywhere; every
harmonic has a sign change. A negative flux is physically meaningless, so
`np.all(phi > 0)` is a rigorous test for "did I land on the right mode" that
costs one pass over the array.

**General principle:** an iterative solver reporting convergence has told you
its iterates stopped moving. It has not told you they stopped at the right
place. Always verify against the problem you actually meant to solve.

---

## Part 4 — Why the tests are built the way they are

### The analytic benchmark

For a bare homogeneous rectangle, separation of variables gives an exact answer:

```
k = nu_Sigma_f / (Sigma_a + D B^2),    B^2 = (pi/a)^2 + (pi/b)^2
```

`B^2` is the **geometric buckling** — pure geometry, no materials. It measures
how sharply the flux must curve to vanish at the boundary. Small reactor →
large buckling → more leakage → smaller `k`. The equation says criticality is a
race between production and (absorption + leakage), which is the whole of
reactor design in one line.

Testing against this checks the *entire* chain — assembly, boundary conditions,
eigensolver — against a number derived independently of the code.

### Why the convergence study, not a single run

Look at Test 1(b), the vacuum-BC table:

```
   n=  40    error =  4.007 pcm    ratio  4.67
   n=  80    error =  0.331 pcm    ratio 12.11
   n= 160    error = -0.588 pcm    ratio  0.56
```

At `n=80` the error is 0.3 pcm and it looks like the code is superb. It is not.
The interior discretization and the Marshak boundary term each carry `O(h^2)`
error with **opposite signs**; near `n=80` they cancel. Refine further and the
error changes sign and grows again.

Had I validated at a single mesh I might have picked `n=80`, declared 0.3 pcm
accuracy, and shipped a code whose real error is several pcm. This is the most
transferable lesson here: **a single accurate result proves nothing.** What you
verify is the *rate* — the error must fall by 4× per halving for a second-order
scheme. Test 1(a) uses the zero-flux BC precisely to isolate one error source
and get clean 4.00 ratios.

### The symmetry test

Quarter core with reflective planes vs. full core: same physics, completely
different mesh, different boundary code paths, 4× the unknowns. They agree to
2e-11 pcm. This is a **consistency test** — it needs no reference answer, which
makes it the kind you can always afford to write. It would catch an indexing
error in the reflective BC, an asymmetry in assembly, or a mesh-mapping bug.

### Test ordering

The four tests run strongest-evidence-first, deliberately:

1. **Analytic** — compared against exact truth.
2. **Symmetry** — internal consistency, no reference needed.
3. **Iteration** — separates *discretization* error from *iteration* error, so
   that when a number is off you know which half to debug.
4. **Physics** — the reflector thermal peak. The weakest evidence: it confirms
   the answer is not absurd, not that it is right.

Test 4 is the one to be most suspicious of. "The plot looks physical" is how
wrong codes survive for years. It is a smoke test, and it is listed last for
that reason.

---

## Part 5 — What to read next

- **Duderstadt & Hamilton, *Nuclear Reactor Analysis*** — the standard text;
  Chapters 4–5 cover everything in Part 1 far more carefully.
- **Lewis & Miller, *Computational Methods of Neutron Transport*** — for when
  diffusion stops being good enough and you want SN or MOC.
- **Stacey, *Nuclear Reactor Physics*** — stronger on nodal methods and the
  discontinuity factors mentioned in Part 1.2.
- **The OpenMC source** — a well-engineered, readable modern Monte Carlo code;
  worth reading even if you stay deterministic.

---

## Appendix — the code README

A compact, verified k-eigenvalue solver for reactor neutronics. Requires only
`numpy` and `scipy` (`matplotlib` optional, for the flux plots).

```bash
python run_tests.py
```

### Files

| file | contents |
|---|---|
| [`diffusion.py`]({{ '/assets/code/diffusion.py' | relative_url }}) | `Material`, `Mesh`, `build_operators`, `solve` — the whole solver |
| [`run_tests.py`]({{ '/assets/code/run_tests.py' | relative_url }}) | four verification cases + a worked reflected-core example |
| `flux_maps.png` | fast flux, thermal flux, and fission source for the example |

### Physics and discretization

Solves the multigroup diffusion equation

```
-div(D_g grad phi_g) + Sigma_r,g phi_g
    = sum_{g'!=g} Sigma_s,g'->g phi_g'  +  (chi_g / k) sum_g' nu_Sigma_f,g' phi_g'
```

- **Mesh-centered finite volume** on a structured, possibly non-uniform,
  Cartesian grid. Interface coupling `2 D1 D2 / (D1 h2 + D2 h1)` comes from
  requiring continuity of flux and current, so material discontinuities are
  handled correctly without extra treatment.
- **Boundary conditions:** `reflective`, `vacuum` (Marshak / linear
  extrapolation), `zero`.
- **Scattering lives inside the loss matrix**, not in a group-sweep loop, so
  up-scattering needs no extra machinery — put entries in the lower triangle
  of `sigma_s` and it works.
- **Transverse leakage** in z via `Material.buckling_b2`, which turns the 2-D
  solve into an approximate 3-D one.
- **Power iteration** with an optional **Wielandt shift**. The shift matters:
  on the example core it cuts 507 iterations to 41 for an identical answer.

### Verification results

| test | result |
|---|---|
| Analytic bare reactor, zero-flux BC | error ratios 3.99 / 4.00 / 4.00 / 4.00 — clean 2nd order |
| Quarter core vs. full core | agree to 2e-11 pcm |
| Wielandt vs. plain power iteration | agree to 2e-4 pcm |
| Reflector thermal peak | present, 0.176 → 0.328 across the interface |

Two things in here are worth internalizing, because they cost me time:

1. **The vacuum-BC error changes sign near n=80.** The interior and boundary
   discretizations each carry O(h²) error of opposite sign, so they cancel at
   one particular mesh and a single run there looks suspiciously perfect. Only
   a refinement study reveals this. Never validate on one mesh.
2. **A Wielandt shift applied too early converges silently to the wrong
   eigenvalue.** The shift must sit above the true k; if it doesn't, modes
   above it get negative λ and the iteration can lock onto the wrong mode and
   then report convergence. Hence the plain-power warmup, and the unshifted
   residual check at the end that verifies the eigenpair independently of
   however it was produced.

### Usage

```python
import numpy as np
from diffusion import Material, Mesh, solve

fuel = Material("fuel",
    D=[1.40, 0.40], sigma_a=[0.010, 0.085],
    sigma_s=[[0.0, 0.020], [0.0, 0.0]],     # [g_from, g_to]
    nu_sigma_f=[0.006, 0.150], chi=[1.0, 0.0])

n, h = 40, 5.0
mesh = Mesh(hx=np.full(n, h), hy=np.full(n, h),
            mat_map=np.zeros((n, n), int), materials=[fuel],
            bc=dict(xmin="reflective", ymin="reflective",
                    xmax="vacuum", ymax="vacuum"))

sol = solve(mesh, verbose=True)
print(sol.k_eff, sol.flux.shape)   # flux is (G, ny, nx)
```

Group 0 is the fastest group. Lengths in cm, cross sections in 1/cm.

### Where to take it next

Roughly in order of payoff:

1. **Real benchmarks.** The 2-D IAEA PWR benchmark and BIBLIS are the standard
   next targets — heterogeneous, published reference values, and they will
   expose mesh-convergence issues this smooth example doesn't.
2. **Assembly discontinuity factors / nodal methods.** Finite-difference
   diffusion needs a fine mesh to be accurate. Real core simulators use nodal
   expansion (NEM) or analytic nodal methods to get node-averaged accuracy on
   one node per assembly. This is the single biggest step toward a practical
   code.
3. **Depletion.** Couple to a Bateman-equation solver (CRAM is the modern
   choice) to burn the fuel and update cross sections over cycle life.
4. **Thermal-hydraulic feedback.** Density and Doppler feedback on the cross
   sections, Picard-iterated with the neutronics. This is where reactor
   physics gets genuinely interesting and where diffusion codes earn their keep.
5. **Beyond diffusion.** Diffusion fails near strong absorbers and voids. SN
   discrete ordinates or method of characteristics is the next rung. The
   eigenvalue machinery you have here carries over unchanged — only the
   transport operator changes.

### Cross-section data

The example data here is made up for testing. For real work you need real
multigroup constants:

- **OpenMC** (open source, Python API) — can generate multigroup libraries
  from continuous-energy ENDF/B data, and is the easiest independent check on
  your k values.
- **NJOY** — the standard tool for processing ENDF/B evaluations into
  multigroup form.
- **Published benchmark specs** (IAEA, BIBLIS, ANL) ship their own cross
  sections, which is why they're the right place to start.

Budget real time for the data pipeline. On most projects like this it is
larger than the solver.
