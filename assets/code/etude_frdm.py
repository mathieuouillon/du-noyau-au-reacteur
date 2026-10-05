"""
etude_frdm.py -- Le modele macroscopique-microscopique, piece par piece,
puis face aux masses mesurees.

Prerequis : micro.npz (python precalcul_micro.py, ~10 min, une seule fois).
Lancer :    python etude_frdm.py    ->  modele_frdm.png
"""

import numpy as np
from donnees_liaison import table_liaison
from modeles_masse import ModeleLineaire, colonnes_M2, colonnes_M3
from nilsson import spectre, hbar_omega, delta_vers_beta2
from strutinsky import Lissage
from frdm import ModeleMacMic, GRILLE_DELTA, TERMES


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


def titre(t):
    print("\n" + "=" * 76 + "\n" + t + "\n" + "=" * 76)


def charger():
    t = table_liaison(mesurees_seulement=True)
    t = t[(t["Z"] >= 8) & (t["N"] >= 8)]
    return t, t["Z"].astype(float), t["N"].astype(float), t["B_MeV"]


# ==========================================================================
def etape1_nilsson():
    titre("1. MODELE DE NILSSON : LES NOMBRES MAGIQUES EMERGENT")
    print("""
  On diagonalise le hamiltonien de Nilsson a deformation nulle et on mesure
  l'ecart entre le dernier niveau occupe et le suivant, pour chaque nombre
  de nucleons. Un grand ecart = une couche fermee = un nombre magique.
  Rien n'est ajuste ici : kappa et mu viennent de la litterature.
""")
    for typ in ("neutron", "proton"):
        e = spectre(0.0, typ)
        g = np.diff(e)
        print(f"  {typ}s, ecart apres n nucleons (en hbar*w0) :")
        ligne = "   "
        for n in (8, 20, 28, 40, 50, 82, 114, 126):
            ligne += f" {n}:{g[n // 2 - 1]:.2f}"
        print(ligne)
    print("""
  Neutrons : 8, 20, 28, 50, 82, 126 ressortent ; 40 (magique pour un
  oscillateur pur) est affaibli par le spin-orbite. C'est le resultat de
  Mayer et Jensen (1949).

  Protons : pas d'ecart a 126, mais un grand ecart a 114. Ce n'est pas une
  erreur : c'est la prediction historique du modele de Nilsson (Nilsson et
  al. 1969), a l'origine de l'idee d'un "ilot de stabilite" des elements
  superlourds autour de Z = 114. Sans consequence ici (Z <= 100).
""")


# ==========================================================================
def etape2_plateau():
    titre("2. STRUTINSKY : LA CONDITION DE PLATEAU")
    print("""
  La correction de couches ne doit pratiquement pas dependre de la largeur
  de lissage gamma. C'est le test de validite de la methode.
""")
    Z, N = 82, 126
    hn, hp = hbar_omega(Z, N, "neutron"), hbar_omega(Z, N, "proton")
    en, ep = spectre(0.0, "neutron"), spectre(0.0, "proton")
    gam, tot = [], []
    print(f"  {'gamma':>8} {'dE neutrons':>13} {'dE protons':>12} {'total':>9}")
    for g in (0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5):
        dn = Lissage(en, g).correction_couches(N) * hn
        dp = Lissage(ep, g).correction_couches(Z) * hp
        gam.append(g)
        tot.append(dn + dp)
        print(f"  {g:>6.1f}hw {dn:>13.2f} {dp:>12.2f} {dn + dp:>9.2f} MeV")
    t = np.array(tot)[1:]
    print(f"""
  Entre 1,0 et 1,5 hbar*w0, la correction du plomb 208 varie de
  {t.max() - t.min():.2f} MeV seulement : le plateau existe. On prend gamma = 1,2.
""")
    return gam, tot


# ==========================================================================
def etape3_courbes(mm):
    titre("3. L'ENERGIE EN FONCTION DE LA FORME")
    print("""
  Pour chaque noyau : E(delta) = E_def_macro(delta) + E_micro(delta).
  La goutte seule prefere toujours la sphere (sa surface coute plus que ce
  que Coulomb rend, tant que x < 1). Les couches peuvent renverser cela.
""")
    courbes = {}
    print(f"  {'noyau':<8} {'macro(0.25)':>12} {'micro(0)':>10} {'micro(0.25)':>12} "
          f"{'delta*':>8}")
    i0 = list(GRILLE_DELTA).index(0.0)
    i25 = list(GRILLE_DELTA).index(0.25)
    for nom, z, n in (("Pb-208", 82, 126), ("Sm-152", 62, 90),
                      ("Er-166", 68, 98), ("U-238", 92, 146)):
        tot, mac = mm.energie_forme(z, n)
        mic = tot - mac
        d, _ = mm.forme(z, n)
        courbes[nom] = (mac, mic, tot)
        print(f"  {nom:<8} {mac[i25]:>12.2f} {mic[i0]:>10.2f} {mic[i25]:>12.2f} {d:>8.3f}")
    print("""
  Pb-208 : la sphere est un puits microscopique profond (-11 MeV) : il
  reste spherique. Er-166, U-238 : la sphere est PENALISEE de 4,5 a 6 MeV
  par les couches, alors qu'une forme allongee en gagne ; ce gain depasse le
  cout macroscopique et le noyau se deforme spontanement.

  C'est l'explication microscopique de la deformation des noyaux : ce
  n'est pas la goutte qui veut s'allonger, ce sont les nucleons de valence
  qui trouvent, dans une forme allongee, des niveaux plus bas.
""")
    return courbes


# ==========================================================================
def etape4_epreuves(t, Z, N, B):
    titre("4. FACE AUX MASSES MESUREES : TROIS EPREUVES")
    rng = np.random.default_rng(1)
    plis = rng.integers(0, 5, len(t))
    bord = np.zeros(len(t), bool)
    for z in np.unique(t["Z"]):
        idx = np.where(t["Z"] == z)[0]
        o = idx[np.argsort(N[idx])]
        bord[o[:2]] = True
        bord[o[-2:]] = True
    c0 = ModeleLineaire("M2", colonnes_M2).ajuster(Z, N, B).coefs

    def fabriques():
        return [
            ("M2  goutte etendue", 8,
             lambda z, n, b: ModeleLineaire("M2", colonnes_M2).ajuster(z, n, b)),
            ("M3  + couches comptees", 11,
             lambda z, n, b: ModeleLineaire("M3", colonnes_M3).ajuster(z, n, b)),
            ("MM  mac-mic, couches seules", 8,
             lambda z, n, b: ModeleMacMic(appariement=False).ajuster(z, n, b, c0)),
            ("MM  mac-mic + appariement BCS", 8,
             lambda z, n, b: ModeleMacMic(appariement=True).ajuster(z, n, b, c0)),
        ]

    print(f"\n  {len(t)} masses mesurees AME2020 (Z, N >= 8). RMS sur B en MeV.\n")
    print(f"  {'modele':<32} {'param.':>6} {'ajust.':>8} {'interp.':>8} {'extrap.':>8}")
    res = []
    for nom, p, fab in fabriques():
        m = fab(Z, N, B)
        r_aj = rms(B - m.predire(Z, N))
        err = []
        for k in range(5):
            a, s = plis != k, plis == k
            mk = fab(Z[a], N[a], B[a])
            err.append(B[s] - mk.predire(Z[s], N[s]))
        r_cv = rms(np.concatenate(err))
        me = fab(Z[~bord], N[~bord], B[~bord])
        r_ex = rms(B[bord] - me.predire(Z[bord], N[bord]))
        print(f"  {nom:<32} {p:>6} {r_aj:>8.3f} {r_cv:>8.3f} {r_ex:>8.3f}")
        res.append((nom, r_aj, r_cv, r_ex))
    print("""
  Les "parametres" comptent uniquement ce qui est AJUSTE sur les masses.
  Les deux modeles mac-mic n'ajustent que les 8 coefficients de la goutte :
  toute la physique des couches vient de kappa, mu et Delta~ = 12/sqrt(A),
  pris dans la litterature. Ils battent pourtant M3 en ajustement et en
  interpolation, alors que M3 ajuste ses 3 coefficients de couches
  directement sur les donnees.

  En extrapolation aux bords, le bilan est mitige : la version avec
  appariement fait mieux que M3, la version sans appariement moins bien.
  Les bords contiennent beaucoup de noyaux legers et de noyaux riches en
  protons, ou les deux faiblesses identifiees en section 5 (noyaux legers,
  couches fermees) pesent le plus.
""")
    return res, c0, bord


# ==========================================================================
def etape5_regions(Z, N, B, mm, mm0):
    titre("5. OU LE MODELE REUSSIT, OU IL ECHOUE")
    A = Z + N
    m2 = ModeleLineaire("M2", colonnes_M2).ajuster(Z, N, B)
    m3 = ModeleLineaire("M3", colonnes_M3).ajuster(Z, N, B)
    R = {"M2": B - m2.predire(Z, N), "M3": B - m3.predire(Z, N),
         "MM": B - mm.predire(Z, N)}
    print(f"\n  {'region':<12} {'noyaux':>7} {'M2':>7} {'M3':>7} {'MM':>7}")
    for lab, msk in (("A < 40", A < 40), ("40-100", (A >= 40) & (A < 100)),
                     ("100-150", (A >= 100) & (A < 150)),
                     ("150-200", (A >= 150) & (A < 200)), ("A >= 200", A >= 200)):
        print(f"  {lab:<12} {msk.sum():>7} " +
              " ".join(f"{rms(R[k][msk]):>7.2f}" for k in ("M2", "M3", "MM")))
    i = np.where((Z == 82) & (N == 126))[0][0]
    r_avec = R["MM"][i]
    r_sans = (B - mm0.predire(Z, N))[i]
    print(f"""
  * TERRES RARES (150-200) : meilleur modele, nettement. C'est la region
    des noyaux fortement deformes, ou le comptage de valence de M3 ne sait
    pas representer l'energie gagnee en se deformant.

  * NOYAUX LOURDS (A >= 200) : moins bon que M3. Le cas d'ecole est le
    plomb 208 : erreur {r_avec:+.2f} MeV avec appariement, {r_sans:+.2f} MeV sans.

    Diagnostic : dans le plomb 208, la couche fermee laisse un tel ecart
    entre niveaux que l'equation du gap BCS n'a plus de solution non nulle.
    Le gap s'EFFONDRE a zero, et la correction penalise le noyau de toute
    l'energie d'appariement moyenne (~1,3 MeV par type de nucleon). Or un
    vrai noyau magique garde des correlations d'appariement : BCS, qui ne
    conserve pas le nombre de particules, les rate. C'est precisement le
    defaut que corrige la methode de Lipkin-Nogami (Lipkin 1960, Nogami
    1964), retenue par FRDM.

  * NOYAUX LEGERS (A < 40) : tous les modeles souffrent ; quelques nucleons
    ne font pas une goutte, et les parametres de Nilsson N = 0, 1 sont
    prolonges sans verification.
""")
    print("  Le coefficient d'appariement MACROSCOPIQUE apres ajustement :")
    print(f"    sans appariement microscopique : {mm0.coefs[4]:6.2f} MeV")
    print(f"    avec appariement microscopique : {mm.coefs[4]:6.2f} MeV")
    print("""
  Avec le BCS et le blocage du nucleon celibataire, le terme empirique
  delta/sqrt(A) devient inutile : l'effet pair-impair des masses est
  reproduit A PARTIR DES NIVEAUX. Le blocage, c'est la raison physique de
  cet effet : un nucleon non apparie occupe seul son niveau et prive
  l'appariement de cet etat.
""")


# ==========================================================================
def etape6_deformations(Z, N, mm):
    titre("6. LES DEFORMATIONS PREDITES")
    d = mm.deformations(Z, N)
    du = mm.forme(92, 146)[0]
    print(f"""
  Fraction des noyaux mesures predits deformes (|delta| > 0,05) : {np.mean(np.abs(d) > 0.05):.0%}
    dont aplatis (oblate) : {np.mean(d < -0.05):.0%} ; allonges (prolate) : {np.mean(d > 0.05):.0%}

  Les noyaux allonges dominent largement, comme dans la nature. C'est un
  resultat connu et non trivial : la preference prolate vient de la
  structure des orbitales de Nilsson pres du niveau de Fermi.

  Uranium 238 : delta* = {du:.3f}, soit beta2 ~ {float(delta_vers_beta2(du)):.3f}.
  Mesure : beta2 = 0,286, tire de B(E2) = 12,09 e^2 b^2 (compilation de
  S. Raman et al., At. Data Nucl. Data Tables 78 (2001) 1).

  Le modele SOUS-ESTIME la deformation de l'uranium. Trois raisons :
    - une seule variable de forme (delta) ; FRDM minimise dans un espace
      a quatre dimensions (eps2, eps3, eps4, eps6), et la deformation
      hexadecapolaire eps4 ajoute au moment quadrupolaire mesure ;
    - le beta2 "experimental" de Raman est deduit de B(E2) avec un rayon
      R0 = 1,2 A^(1/3) fm et une surface abrupte : c'est une convention, pas
      la meme grandeur que notre delta ;
    - la conversion delta -> beta2 est au premier ordre seulement.
  L'accord qualitatif (spherique / deforme, prolate / oblate) est bon ;
  l'accord quantitatif sur beta2 ne l'est pas, et il ne faut pas le
  surinterpreter dans un sens ou dans l'autre.
""")
    return d


# ==========================================================================
def etape7_distance(t, Z, N, B, c0, K=8):
    titre("7. EXTRAPOLATION : LA PHYSIQUE TIENT-ELLE LA DISTANCE ?")
    test = np.zeros(len(t), bool)
    dist = np.zeros(len(t), int)
    for z in np.unique(t["Z"]):
        idx = np.where(t["Z"] == z)[0]
        if len(idx) < K + 6:
            continue
        o = idx[np.argsort(N[idx])]
        test[o[-K:]] = True
        dist[o[-K:]] = (N[o[-K:]] - N[o[-K - 1]]).astype(int)
    modeles = {
        "M3": ModeleLineaire("M3", colonnes_M3).ajuster(Z[~test], N[~test], B[~test]),
        "MM": ModeleMacMic().ajuster(Z[~test], N[~test], B[~test], c0),
    }
    courbes = {}
    print(f"\n  On retire les {K} isotopes les plus riches en neutrons de chaque element")
    print(f"  ({test.sum()} noyaux), on ajuste sur le reste, on les predit.\n")
    print(f"  {'d':>3} {'M3':>8} {'MM':>8}")
    for k, m in modeles.items():
        e = B[test] - m.predire(Z[test], N[test])
        courbes[k] = [(dd, rms(e[dist[test] == dd])) for dd in range(1, K + 1)
                      if (dist[test] == dd).sum() >= 5]
    for (dd, a), (_, b) in zip(courbes["M3"], courbes["MM"]):
        print(f"  {dd:>3} {a:>8.3f} {b:>8.3f}")
    a1, a8 = courbes["M3"][0][1], courbes["M3"][-1][1]
    b1, b8 = courbes["MM"][0][1], courbes["MM"][-1][1]
    print(f"""
  M3 : x{a8 / a1:.2f} de d=1 a d=8.   MM : x{b8 / b1:.2f}.
  (A comparer a la correction statistique M4 de MODELES_MASSE.md : x3,5.)
  Le modele mac-mic est meilleur a toutes les distances et se degrade un
  peu moins vite. L'ecart reste modeste : ses 8 coefficients de goutte
  sont ajustes eux aussi. Mais ses couches viennent de niveaux CALCULES et
  non d'un ajustement : cette partie-la n'a rien appris de la region
  connue, et n'a donc rien a oublier en s'en eloignant.
""")
    return courbes


# ==========================================================================
def etape8_fissilite(mm):
    titre("8. RETOUR A LA FISSION : LE PARAMETRE DE FISSILITE")
    c = mm.coefs
    print("""
  Bohr et Wheeler (1939) : la goutte spherique est stable tant que
  x = E_coulomb / (2 E_surface) < 1. Avec les coefficients AJUSTES ici :
""")
    print(f"  {'noyau':<9} {'E_surface':>10} {'E_coulomb':>10} {'x':>7}")
    for nom, z, n in (("Pb-208", 82, 126), ("U-236", 92, 144), ("Pu-240", 94, 146),
                      ("Fm-256", 100, 156)):
        A, I = z + n, n - z
        Es = c[1] * A ** (2 / 3) - c[5] * I ** 2 / A ** (4 / 3)
        Ec = c[2] * z * (z - 1) / A ** (1 / 3) - c[6] * z ** (4 / 3) / A ** (1 / 3)
        print(f"  {nom:<9} {Es:>10.0f} {Ec:>10.0f} {Ec / (2 * Es):>7.3f}")
    print("""
  L'uranium 236 est a x ~ 0,7 : la goutte est encore stable, mais la
  barriere qui la retient n'est que de quelques MeV. C'est pour cela que
  les 6,5 MeV apportes par la capture d'un neutron (voir FISSION_URANIUM.md)
  suffisent a la faire fissionner.

  Calculer la barriere elle-meme demande des formes tres allongees, avec un
  col (necking) : hors de portee des spheroides de ce code. C'est l'un des
  usages principaux de FRDM (barrieres a double bosse des actinides).
""")


# ==========================================================================
def graphiques(t, Z, N, B, plateau, courbes, res, dist, d_pred):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(14, 15.5))
    gs = fig.add_gridspec(3, 2, hspace=0.38, wspace=0.24)

    # (a) diagramme de Nilsson, neutrons
    ax = fig.add_subplot(gs[0, 0])
    ds = np.linspace(-0.3, 0.45, 61)
    E = np.array([spectre(d, "neutron", 10) for d in ds])
    for k in range(E.shape[1]):
        if 5.0 < E[0, k] < 8.6 or 5.0 < E[-1, k] < 8.6:
            ax.plot(ds, E[:, k], lw=0.7, color="#1f4e79")
    e0 = spectre(0.0, "neutron", 10)
    for n in (50, 82, 126):
        y = 0.5 * (e0[n // 2 - 1] + e0[n // 2])
        ax.text(-0.29, y, str(n), color="crimson", fontsize=9, va="center",
                fontweight="bold")
    ax.axvline(0, color="k", lw=0.6, ls=":")
    ax.set_ylim(5.0, 8.6)
    ax.set_xlabel("deformation delta  (< 0 aplati, > 0 allonge)")
    ax.set_ylabel("energie  [hbar w0]")
    ax.set_title("(a) Diagramme de Nilsson, neutrons\n(les couches 50, 82, 126 se brisent "
                 "a la deformation)")

    # (b) plateau
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(plateau[0], plateau[1], "o-", color="#1f4e79", lw=2)
    ax.axvspan(1.0, 1.5, color="seagreen", alpha=0.12)
    t = np.array(plateau[1])[1:]
    ax.text(1.25, t.max() + 0.12, f"plateau : variation de {t.max() - t.min():.2f} MeV\n"
            "sur une correction de ~14 MeV\n(axe tres dilate)", ha="center",
            fontsize=8.5, color="seagreen")
    ax.set_xlabel("largeur de lissage gamma  [hbar w0]")
    ax.set_ylabel("correction de couches  [MeV]")
    ax.set_title("(b) Strutinsky : condition de plateau\n(plomb 208, spherique)")
    ax.grid(alpha=0.3)

    # (c) E(delta)
    ax = fig.add_subplot(gs[1, 0])
    couleurs = {"Pb-208": "#555555", "Sm-152": "#2a9d8f", "Er-166": "#e76f51",
                "U-238": "#1f4e79"}
    for nom, (mac, mic, tot) in courbes.items():
        ax.plot(GRILLE_DELTA, tot - tot[list(GRILLE_DELTA).index(0.0)], "-",
                color=couleurs[nom], lw=2.2, label=nom)
        ax.plot(GRILLE_DELTA, mac, ":", color=couleurs[nom], lw=1.2)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("deformation delta")
    ax.set_ylabel("E(delta) - E(0)  [MeV]")
    ax.set_title("(c) Energie de deformation\n(pointilles : goutte seule ; "
                 "traits pleins : goutte + couches + appariement)")
    ax.set_ylim(-8, 10)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)

    # (d) carte des deformations
    ax = fig.add_subplot(gs[1, 1])
    sc = ax.scatter(N, Z, c=d_pred, s=5, cmap="RdBu_r", vmin=-0.3, vmax=0.3)
    for m in (28, 50, 82, 126):
        ax.axvline(m, color="k", ls=":", lw=0.6)
    for m in (28, 50, 82):
        ax.axhline(m, color="k", ls=":", lw=0.6)
    ax.set_xlabel("N")
    ax.set_ylabel("Z")
    ax.set_title("(d) Deformation predite delta*\n(rouge allonge, bleu aplati, "
                 "blanc spherique)")
    fig.colorbar(sc, ax=ax, label="delta*")

    # (e) epreuves
    ax = fig.add_subplot(gs[2, 0])
    noms = ["M2", "M3", "MM\ncouches", "MM\n+ appar."]
    x = np.arange(len(res))
    for i, (lab, col) in enumerate((("ajustement", "#9db4cf"),
                                    ("interpolation", "#1f4e79"),
                                    ("extrapolation", "#e76f51"))):
        ax.bar(x + (i - 1) * 0.27, [r[i + 1] for r in res], 0.27, label=lab, color=col)
    ax.axhline(0.669, color="seagreen", ls="--", lw=1)
    ax.text(-0.4, 0.75, "FRDM(1992) : 0,669 MeV", color="seagreen", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels(noms)
    ax.set_ylabel("ecart RMS sur B  [MeV]")
    ax.set_title("(e) Trois epreuves (masses mesurees AME2020)")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, axis="y")

    # (f) distance
    ax = fig.add_subplot(gs[2, 1])
    for k, col in (("M3", "#9db4cf"), ("MM", "#1f4e79")):
        ax.plot([c[0] for c in dist[k]], [c[1] for c in dist[k]], "o-", color=col,
                lw=2, label={"M3": "M3 couches comptees",
                             "MM": "mac-mic (couches calculees)"}[k])
    ax.set_xlabel("neutrons au-dela du dernier isotope connu")
    ax.set_ylabel("ecart RMS sur B  [MeV]")
    ax.set_title("(f) Extrapolation vers les noyaux riches en neutrons")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)

    fig.text(0.5, 0.005, "Donnees : AME2020 (masses mesurees). Nilsson : Bengtsson & "
             "Ragnarsson 1985 ; Strutinsky 1967 ; Brack et al. 1972.",
             ha="center", fontsize=8.5, color="#444444")
    fig.savefig("modele_frdm.png", dpi=120, bbox_inches="tight")
    print("\n  --> modele_frdm.png ecrit")


# ==========================================================================
def main():
    t, Z, N, B = charger()
    etape1_nilsson()
    plateau = etape2_plateau()
    res, c0, _ = etape4_epreuves(t, Z, N, B)
    mm = ModeleMacMic(appariement=True).ajuster(Z, N, B, c0)
    mm0 = ModeleMacMic(appariement=False).ajuster(Z, N, B, c0)
    titre("COEFFICIENTS MACROSCOPIQUES AJUSTES (MeV)")
    for nom, v in zip(TERMES, mm.coefs):
        print(f"  {nom:<16} {v:>9.3f}")
    courbes = etape3_courbes(mm)
    etape5_regions(Z, N, B, mm, mm0)
    d_pred = etape6_deformations(Z, N, mm)
    dist = etape7_distance(t, Z, N, B, c0)
    etape8_fissilite(mm)
    graphiques(t, Z, N, B, plateau, courbes, res, dist, d_pred)


if __name__ == "__main__":
    main()
