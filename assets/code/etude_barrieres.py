"""
etude_barrieres.py -- Les barrieres de fission, de la goutte liquide a la
double bosse des actinides.

Cinq questions, dans l'ordre :
  1. Le calcul des formes est-il juste ?        (sphere, spheroides exacts,
                                                 petites deformations)
  2. Que donne la goutte seule ?                (Bohr-Wheeler, point selle,
                                                 MS-LD contre LSD)
  3. La goutte prefere-t-elle un partage symetrique ?
  4. Que changent les couches ?                 (Nilsson + Strutinsky + BCS
                                                 le long des spheroides)
  5. Ou s'arrete ce modele ?                    (spheroides contre vallee
                                                 avec col)

Lancer :  python etude_barrieres.py   ->  barrieres_fission.png   (~1 min)

Valeurs de reference des barrieres : RIPL-3 (R. Capote et al., Nucl. Data
Sheets 110 (2009) 3107), telles que reprises dans la table 1 de W. Ryssens
et al., Eur. Phys. J. A 59 (2023) 96 ; energies des isomeres de fission :
M. Samyn et al., Phys. Rev. C 70 (2004) 044309, via la meme table.
"""

import numpy as np
from barrieres import (Forme, Carte, fonctions_forme, energies_sphere,
                       energie_deformation, fissilite, barriere_reduite,
                       bohr_wheeler, penetrabilite)
from frdm import rapport_axes, surface_coulomb
from nilsson import spectre, hbar_omega
from strutinsky import Lissage, correction_appariement

# (Z, N, nom) ; references : (barriere primaire E_I, secondaire E_II, isomere)
NOYAUX = [(92, 142, "U-234"), (92, 144, "U-236"), (92, 146, "U-238"),
          (94, 146, "Pu-240"), (94, 148, "Pu-242"), (96, 148, "Cm-244")]
REFERENCE = {"U-234": (5.50, 4.80, None), "U-236": (5.67, 5.00, 2.30),
             "U-238": (6.30, 5.50, 2.60), "Pu-240": (6.05, 5.45, 2.25),
             "Pu-242": (5.85, 5.05, None), "Cm-244": (6.18, 5.10, 1.04)}

N_MAX = 20                      # couches d'oscillateur (convergence : etape 4)
DELTAS = np.round(np.arange(-0.10, 0.6001, 0.025), 4)


def titre(t):
    print()
    print("=" * 76)
    print(t)
    print("=" * 76)


# ==========================================================================
def verifications():
    titre("1. LE CALCUL DES FORMES EST-IL JUSTE ?")
    Bs, Bk, Bc = fonctions_forme(Forme(1.0))
    print(f"\n  sphere : B_s = {Bs:.6f}, B_k = {Bk:.6f}, B_c = {Bc:.6f}   (attendu : 1)")
    print("\n  spheroides (h = -(c-1)/4) contre formules exactes (frdm.py) :")
    print(f"  {'c':>5} {'B_s calc':>10} {'B_s exact':>10} {'B_c calc':>10} {'B_c exact':>10}")
    for c in (1.2, 1.5, 2.0):
        q = c ** 1.5
        d = 3 * (q * q - 1) / (2 + 4 * q * q)
        bs, _, bc = fonctions_forme(Forme(c, -(c - 1) / 4))
        es, ec = surface_coulomb(d)
        print(f"  {c:5.2f} {bs:10.6f} {es:10.6f} {bc:10.6f} {ec:10.6f}")
    print("\n  petites deformations, R = R0 (1 + a2 P2) avec a2 = c - 1 :")
    print("  Bohr et Wheeler : B_s - 1 = (2/5) a2^2,  B_c - 1 = -(1/5) a2^2")
    a2 = 0.02
    bs, _, bc = fonctions_forme(Forme(1 + a2, -a2 / 4))
    print(f"  a2 = {a2} : (B_s - 1)/a2^2 = {(bs - 1) / a2**2:.3f},"
          f"  (B_c - 1)/a2^2 = {(bc - 1) / a2**2:.3f}")
    print("\n  convergence de la quadrature (c = 1,6, h = 0) :")
    for n in (48, 96, 160):
        bs, bk, bc = fonctions_forme(Forme(1.6, 0.0), n_z=n, n_psi=n)
        print(f"    n = {n:3d} : B_s = {bs:.7f}  B_k = {bk:.7f}  B_c = {bc:.7f}")


# ==========================================================================
def goutte(carte):
    titre("2. LA GOUTTE SEULE : BOHR-WHEELER, POINT SELLE, MS-LD CONTRE LSD")
    print("""
  Barriere reduite xi(x) = E_f / E_s d'une goutte sans courbure, calculee
  sur la carte (c, h), contre le developpement de Bohr et Wheeler
  (98/135)(1 - x)^3, exact quand x -> 1 :
""")
    print(f"  {'x':>5} {'xi calcule':>12} {'Bohr-Wheeler':>13} {'rapport':>8} {'c selle':>8}")
    xs = (0.75, 0.80, 0.85, 0.90, 0.95)
    xis = []
    for x in xs:
        xi, c = barriere_reduite(carte, x)
        xis.append(xi)
        print(f"  {x:5.2f} {xi:12.5f} {bohr_wheeler(x):13.5f} {xi / bohr_wheeler(x):8.3f} {c:8.3f}")
    print("""
  En dessous de x ~ 0,7, le point selle s'approche de la scission : la vallee
  "minimum en h" saute alors vers la branche des fragments separes, et le
  maximum le long de c ne donne plus qu'une borne INFERIEURE de la barriere.
  On s'en tient donc aux x >= 0,75 (actinides).""")

    print("""
  Barrieres macroscopiques avec les deux jeux de parametres publies
  (Pomorski et Dudek 2003). B_f = max sur c du min sur h de E_def(c, h).
""")
    print(f"  {'noyau':<8} {'x MS-LD':>8} {'B_f MS-LD':>10} {'x LSD':>7} {'B_f LSD':>8} {'c selle':>8} {'h selle':>8}")
    lignes = {}
    for Z, N, nom in NOYAUX + [(98, 154, "Cf-252")]:
        bm = carte.barriere(Z, N, "MS-LD")
        bl = carte.barriere(Z, N, "LSD")
        lignes[nom] = (bm[0], bl[0])
        print(f"  {nom:<8} {fissilite(Z, N, 'MS-LD'):8.3f} {bm[0]:10.2f} "
              f"{fissilite(Z, N, 'LSD'):7.3f} {bl[0]:8.2f} {bl[1]:8.3f} {bl[2]:8.2f}")
    return np.array(xs), np.array(xis), lignes


# ==========================================================================
def asymetrie(carte):
    titre("3. LA GOUTTE PREFERE-T-ELLE UN PARTAGE SYMETRIQUE ?")
    Z, N = 92, 144
    _, cs, hs = carte.barriere(Z, N, "LSD")
    c = round(cs / 0.025) * 0.025
    print(f"\n  U-236, LSD, forme du point selle (c = {c:.3f}, h = {hs:.2f}), asymetrie alpha :")
    alphas = np.array([0.0, 0.05, 0.10, 0.15, 0.20])
    E = []
    for a in alphas:
        fo = Forme(c, hs, a)
        E.append(energie_deformation(Z, N, fonctions_forme(fo), "LSD"))
        print(f"    alpha = {a:4.2f} : E_def = {E[-1]:6.2f} MeV")
    print("\n  L'energie croit avec alpha : au point selle, la goutte est STABLE")
    print("  contre l'asymetrie. Elle predit une fission symetrique.")
    return c, hs, alphas, np.array(E)


# ==========================================================================
def micro(Z, N, gamma=1.2):
    """dE_couches et dE_couches + dE_appariement (MeV) le long de DELTAS."""
    A = Z + N
    sh, tot = [], []
    for d in DELTAS:
        s = p = 0.0
        for t, n in (("proton", Z), ("neutron", N)):
            L = Lissage(spectre(d, t, N_MAX), gamma)
            hw = hbar_omega(Z, N, t)
            s += L.correction_couches(n) * hw
            p += correction_appariement(L, n, hw, A)
        sh.append(s)
        tot.append(s + p)
    return np.array(sh), np.array(tot)


def couches():
    titre("4. QUE CHANGENT LES COUCHES ? LA DOUBLE BOSSE")
    print("""
  E(delta) = E_goutte LSD(spheroide delta) + dE_couches + dE_appariement,
  niveaux de Nilsson (nilsson.py), correction de Strutinsky et BCS
  (strutinsky.py). Spheroides seulement, comme dans frdm.py.
""")
    print("  Convergence de la base (U-236, dE_couches + dE_app., MeV) :")
    print(f"  {'N_max':>6} " + " ".join(f"{d:6.2f}" for d in (0.2, 0.4, 0.5, 0.55, 0.6)))
    for nm in (14, 18, 20):
        v = []
        for d in (0.2, 0.4, 0.5, 0.55, 0.6):
            v.append(sum(Lissage(spectre(d, t, nm), 1.2).correction_couches(n) * hbar_omega(92, 144, t)
                         + correction_appariement(Lissage(spectre(d, t, nm), 1.2), n,
                                                  hbar_omega(92, 144, t), 236)
                         for t, n in (("proton", 92), ("neutron", 144))))
        print(f"  {nm:6d} " + " ".join(f"{x:6.2f}" for x in v))
    print("  -> converge a ~0,1 MeV jusqu'a delta = 0,55 ; on s'arrete a 0,60.")

    Bsph = []
    for d in DELTAS:
        c = rapport_axes(d) ** (2 / 3)
        Bsph.append(fonctions_forme(Forme(c, -(c - 1) / 4)))

    print(f"\n  {'noyau':<8} {'delta_fond':>10} {'B_A':>6} {'d_A':>6} {'E_iso':>6} {'d_iso':>6}"
          f"   {'ref. E_I':>8} {'E_II':>6} {'E_iso':>6}")
    courbes = {}
    for Z, N, nom in NOYAUX:
        sh, tot = micro(Z, N)
        mac = np.array([energie_deformation(Z, N, b, "LSD") for b in Bsph])
        E = mac + tot
        i0 = int(np.argmin(np.where(DELTAS < 0.35, E, np.inf)))
        m2 = (DELTAS > 0.38) & (DELTAS < 0.6)
        i2 = int(np.where(m2)[0][np.argmin(E[m2])])
        iA = i0 + int(np.argmax(E[i0:i2 + 1]))
        egs = E[i0]
        ref = REFERENCE[nom]
        iso = f"{ref[2]:6.2f}" if ref[2] else "     -"
        print(f"  {nom:<8} {DELTAS[i0]:10.3f} {E[iA] - egs:6.2f} {DELTAS[iA]:6.3f} "
              f"{E[i2] - egs:6.2f} {DELTAS[i2]:6.3f}   {ref[0]:8.2f} {ref[1]:6.2f} {iso}")
        courbes[nom] = dict(mac=mac, sh=sh, tot=tot, E=E - egs, i0=i0, iA=iA, i2=i2)
    print("""
  B_A : barriere interne ; E_iso : fond du second puits (isomere de fission),
  tous deux mesures depuis l'etat fondamental du modele. References : RIPL-3
  (barriere primaire E_I = la plus haute, secondaire E_II = la plus basse).""")
    return Bsph, courbes


# ==========================================================================
def limites(carte, Bsph):
    titre("5. OU S'ARRETE CE MODELE ? SPHEROIDES CONTRE VALLEE AVEC COL")
    Z, N = 92, 144
    Ev, hv = carte.vallee(Z, N, "LSD")
    Qv = np.array([Forme(c, h).quadrupole() if np.isfinite(e) else np.nan
                   for c, h, e in zip(carte.cs, hv, Ev)])
    print("""
  Energie de la goutte (LSD, U-236) a moment quadrupolaire Q2 egal :
  spheroide contre meilleure forme (c, h) de la vallee de fission.
""")
    print(f"  {'delta':>6} {'Q2':>6} {'E spheroide':>12} {'E vallee':>9} {'ecart':>7}")
    res = []
    for d in (0.3, 0.4, 0.5, 0.6):
        i = int(np.argmin(np.abs(DELTAS - d)))
        c = rapport_axes(d) ** (2 / 3)
        q = Forme(c, -(c - 1) / 4).quadrupole()
        es = energie_deformation(Z, N, Bsph[i], "LSD")
        ok = np.isfinite(Qv) & (carte.cs < 1.9)
        ev = float(np.interp(q, Qv[ok], Ev[ok]))
        res.append((d, q, es, ev))
        print(f"  {d:6.2f} {q:6.3f} {es:12.2f} {ev:9.2f} {es - ev:7.2f}")
    print("""
  Au-dela de delta ~ 0,45, un ellipsoide coute bien plus cher qu'une forme
  qui peut se creuser : la barriere externe et la profondeur du second puits
  ne sont pas calculables ici.""")

    titre("PENETRABILITE D'UNE BARRIERE PARABOLIQUE (Hill et Wheeler 1953)")
    print("\n  T = 1 / (1 + exp(2 pi (B_f - E) / hbar_omega)) :")
    print(f"  {'E - B_f (MeV)':>14} " + " ".join(f"{x:7.1f}" for x in (-2, -1, -0.5, 0, 0.5, 1)))
    for hw in (0.5, 1.0):
        T = penetrabilite(np.array([-2, -1, -0.5, 0, 0.5, 1.0]), 0.0, hw)
        print(f"  hbar_omega={hw:3.1f} " + " ".join(f"{t:7.1e}" for t in T))
    return res


# ==========================================================================
def figure(carte, xs, xis, asym, Bsph, courbes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Z, N = 92, 144
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
    Ev, hv = carte.vallee(Z, N, "LSD")
    Bf, cs, hs = carte.barriere(Z, N, "LSD")

    # (a) formes le long de la vallee
    a = ax[0, 0]
    choix = [1.0, 1.25, round(cs / 0.025) * 0.025, 1.75, 1.95]
    for k, c in enumerate(choix):
        i = int(np.argmin(np.abs(carte.cs - c)))
        fo = Forme(carte.cs[i], hv[i])
        z = np.linspace(-fo.c, fo.c, 400)
        r = np.sqrt(np.maximum(fo.f(z), 0))
        dx = 2.6 * k
        a.fill_between(z, -r - dx, r - dx, color=plt.cm.viridis(k / 5), alpha=0.8)
        etiq = f"c = {fo.c:.2f}" + ("  (point selle)" if k == 2 else "")
        a.text(2.15, -dx, etiq, va="center", fontsize=9)
    a.set_xlim(-2.2, 3.6)
    a.set_ylim(-2.6 * 4 - 1.3, 1.3)
    a.set_aspect("equal")
    a.axis("off")
    a.set_title("(a) U-236 : formes le long de la vallee de fission (LSD)")

    # (b) carte E(c, h)
    a = ax[0, 1]
    E = carte.energie(Z, N, "LSD")
    niveaux = np.arange(-6, 14.1, 1.0)
    cf = a.contourf(carte.cs, carte.hs, E.T, levels=niveaux, cmap="RdYlBu_r", extend="both")
    a.contour(carte.cs, carte.hs, E.T, levels=niveaux, colors="k", linewidths=0.3)
    av = carte.cs <= 1.85
    a.plot(carte.cs[av], hv[av], "k--", lw=1.2, label="vallee (min. en h)")
    a.plot([cs], [hs], "w*", ms=14, mec="k", label=f"point selle, {Bf:.2f} MeV")
    c_sph = np.linspace(1, 2.2, 50)
    a.plot(c_sph, -(c_sph - 1) / 4, color="purple", lw=1, label="spheroides")
    a.set_ylim(carte.hs[0], carte.hs[-1])
    a.set_xlabel("elongation c")
    a.set_ylabel("parametre de col h")
    a.legend(fontsize=8, loc="upper left")
    fig.colorbar(cf, ax=a, label="E_def (MeV)")
    a.set_title("(b) U-236 : surface d'energie de la goutte (LSD)")

    # (c) profil de la vallee
    a = ax[0, 2]
    av = carte.cs <= 1.85
    for jeu, st in (("MS-LD", "--"), ("LSD", "-")):
        for (z, n, nom, col) in ((92, 144, "U-236", "C0"), (98, 154, "Cf-252", "C3")):
            ev, _ = carte.vallee(z, n, jeu)
            a.plot(carte.cs[av], ev[av], st, color=col, label=f"{nom}, {jeu}")
    a.axhline(0, color="gray", lw=0.5)
    a.set_ylim(-3, 9)
    a.set_xlabel("elongation c")
    a.set_ylabel("E_def le long de la vallee (MeV)")
    a.legend(fontsize=8)
    a.set_title("(c) Barriere de la goutte : deux jeux de parametres")

    # (d) Bohr-Wheeler
    a = ax[1, 0]
    xx = np.linspace(0.7, 1.0, 100)
    a.semilogy(xx, bohr_wheeler(xx), "k-", label="Bohr-Wheeler (98/135)(1-x)^3")
    a.semilogy(xs, xis, "o", color="C1", label="calcul (c, h)")
    a.set_xlabel("parametre de fissilite x")
    a.set_ylabel("E_f / E_s")
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    a.set_title("(d) Barriere reduite de la goutte")

    # (e) correction microscopique
    a = ax[1, 1]
    cu = courbes["U-236"]
    a.plot(DELTAS, cu["sh"], "C0--", label="couches (Strutinsky)")
    a.plot(DELTAS, cu["tot"], "C0-", label="couches + appariement (BCS)")
    a.plot(DELTAS, cu["mac"], "C3:", label="goutte LSD (spheroides)")
    a.axhline(0, color="gray", lw=0.5)
    a.set_xlabel("deformation delta")
    a.set_ylabel("MeV")
    a.legend(fontsize=8)
    a.set_title("(e) U-236 : les deux morceaux de l'energie")

    # (f) energie totale
    a = ax[1, 2]
    for nom, col in (("U-236", "C0"), ("U-238", "C2"), ("Pu-240", "C3")):
        c_ = courbes[nom]
        a.plot(DELTAS, c_["E"], color=col, label=nom)
        a.plot(DELTAS[c_["iA"]], c_["E"][c_["iA"]], "^", color=col)
        a.plot(DELTAS[c_["i2"]], c_["E"][c_["i2"]], "v", color=col)
    a.axvspan(0.45, 0.6, color="gray", alpha=0.15)
    a.text(0.525, 8.5, "spheroides\npeu fiables", ha="center", fontsize=8)
    a.set_ylim(-1, 10)
    a.axhline(0, color="gray", lw=0.5)
    a.set_xlabel("deformation delta")
    a.set_ylabel("E - E(fondamental) (MeV)")
    a.legend(fontsize=8, loc="upper left")
    a.set_title("(f) Double bosse : barriere interne et isomere")

    fig.tight_layout()
    fig.savefig("barrieres_fission.png", dpi=110)
    print("\n  --> barrieres_fission.png ecrit")


# ==========================================================================
if __name__ == "__main__":
    verifications()
    carte = Carte(np.round(np.arange(1.0, 2.501, 0.025), 4),
                  np.round(np.arange(-0.30, 0.501, 0.02), 4))
    xs, xis, _ = goutte(carte)
    asym = asymetrie(carte)
    Bsph, courbes = couches()
    limites(carte, Bsph)
    figure(carte, xs, xis, asym, Bsph, courbes)
