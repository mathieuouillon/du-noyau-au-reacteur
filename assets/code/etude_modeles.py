"""
etude_modeles.py -- Comparer les modeles de masse, et surtout les juger sur
ce qu'ils PREDISENT, pas sur ce qu'ils reproduisent.

Trois epreuves de difficulte croissante :
  1. ajustement    : erreur sur les noyaux qui ont servi a l'ajustement
  2. interpolation : validation croisee a 5 plis (noyaux tires au hasard)
  3. extrapolation : on retire les noyaux les plus exotiques de chaque
                     element, et on les predit

Lancer :  python etude_modeles.py   ->  modeles_masse.png
"""

import numpy as np
from donnees_liaison import table_liaison
from modeles_masse import (tous_les_modeles, ModeleLineaire, ModeleCorrige,
                           colonnes_M1, colonnes_M2, colonnes_M3, valence)
from fission import (A_VOLUME, A_SURFACE, A_COULOMB, A_ASYMETRIE,
                     A_APPARIEMENT)


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


def titre(texte):
    print()
    print("=" * 76)
    print(texte)
    print("=" * 76)


def charger():
    t = table_liaison(mesurees_seulement=True)
    t = t[(t["Z"] >= 8) & (t["N"] >= 8)]
    return t, t["Z"].astype(float), t["N"].astype(float), t["B_MeV"]


def masque_bords(t, k=2):
    """Les k isotopes les plus riches ET les k plus pauvres en neutrons de
    chaque element : la frontiere de la carte connue."""
    bord = np.zeros(len(t), bool)
    for z in np.unique(t["Z"]):
        idx = np.where(t["Z"] == z)[0]
        o = idx[np.argsort(t["N"][idx])]
        bord[o[:k]] = True
        bord[o[-k:]] = True
    return bord


# ==========================================================================
def epreuves(t, Z, N, B):
    titre("1. TROIS EPREUVES : AJUSTER, INTERPOLER, EXTRAPOLER")
    print(f"""
  {len(t)} masses MESUREES d'AME2020, Z >= 8 et N >= 8. Metrique : ecart
  RMS sur l'energie de liaison TOTALE B, en MeV (standard du domaine).
""")
    rng = np.random.default_rng(1)
    plis = rng.integers(0, 5, len(t))
    bord = masque_bords(t)

    resultats = []
    print(f"  {'modele':<32} {'param.':>7} {'ajust.':>8} {'interp.':>8} "
          f"{'extrap.':>8} {'extrap/interp':>14}")
    for m in tous_les_modeles():
        m.ajuster(Z, N, B)
        r_aj = rms(B - m.predire(Z, N))
        p = m.nombre_parametres()

        err = []
        for k in range(5):
            a, s = plis != k, plis == k
            m.ajuster(Z[a], N[a], B[a])
            err.append(B[s] - m.predire(Z[s], N[s]))
        r_cv = rms(np.concatenate(err))

        m.ajuster(Z[~bord], N[~bord], B[~bord])
        r_ex = rms(B[bord] - m.predire(Z[bord], N[bord]))

        etiq = f"~{len(t)}" if p > 100 else str(p)
        print(f"  {m.nom:<32} {etiq:>7} {r_aj:>8.3f} {r_cv:>8.3f} "
              f"{r_ex:>8.3f} {r_ex / r_cv:>13.2f}x")
        resultats.append((m.nom, r_aj, r_cv, r_ex))

    print(f"""
  ({bord.sum()} noyaux de bord : les 2 plus riches et les 2 plus pauvres en
  neutrons de chaque element, retires de l'ajustement puis predits.)

  Trois lecons :

  * Pour M0-M3, ajustement et interpolation sont quasi identiques. Avec 5
    a 11 parametres pour {len(t)} noyaux, aucun risque de sur-apprentissage :
    le modele ne peut pas "memoriser" les donnees.

  * Chaque piece de physique paie : du manuel au modele avec couches,
    l'erreur est divisee par 3. Le terme de couches a lui seul (M2 -> M3)
    divise l'erreur par 2 avec 3 parametres de plus.

  * M4 ecrase tout en interpolation, mais c'est lui qui se degrade le plus
    en extrapolation. Ce rapport extrap/interp est le chiffre a surveiller
    quand on choisit un modele pour PREDIRE.
""")
    return resultats


# ==========================================================================
def coefficients(t, Z, N, B):
    titre("2. CE QUE DISENT LES COEFFICIENTS")
    m1 = ModeleLineaire("M1", colonnes_M1).ajuster(Z, N, B)
    manuel = [A_VOLUME, A_SURFACE, A_COULOMB, A_ASYMETRIE, A_APPARIEMENT]
    print("\n  M1 : la goutte liquide reajustee sur AME2020\n")
    print(f"  {'terme':<14} {'manuel':>10} {'reajuste':>10}")
    for nom, a, b in zip(m1.termes, manuel, m1.coefs):
        print(f"  {nom:<14} {a:>10.3f} {b:>10.3f}")
    print("""
  Les coefficients bougent peu : les valeurs de manuel etaient deja
  ajustees, sur des tables plus anciennes et avec Z^2 au lieu de Z(Z-1).
  Le gain de M0 a M1 vient de la, pas d'une physique nouvelle.
""")
    m3 = ModeleLineaire("M3", colonnes_M3).ajuster(Z, N, B)
    print("  M3 : tous les termes, en MeV\n")
    print(f"  {'terme':<16} {'coefficient':>12}   lecture")
    lectures = {
        "volume": "energie de liaison d'un nucleon au coeur du noyau",
        "surface": "cout des nucleons en surface",
        "Coulomb": "repulsion entre protons",
        "asymetrie": "cout de l'exces de neutrons (Pauli)",
        "appariement": "bonus des noyaux pair-pair",
        "sym. surface": "l'asymetrie coute MOINS en surface (peau de n)",
        "echange coul.": "correction quantique a Coulomb",
        "Wigner": "surliaison des noyaux N ~ Z",
        "couches S": "penalite d'eloignement des couches fermees",
        "couches S^2": "courbure",
        "couches p-n": "interaction p-n de valence (deformation)",
    }
    for nom, c in zip(m3.termes, m3.coefs):
        print(f"  {nom:<16} {c:>12.4f}   {lectures[nom]}")
    J = m3.coefs[m3.termes.index("asymetrie")]
    Q = m3.coefs[m3.termes.index("sym. surface")]
    print(f"""
  SYMETRIE. Le volume coute J = {J:.1f} MeV, mais la surface en rend
  Q = {Q:.1f} MeV : le coefficient effectif d'un noyau fini est
      a_sym(A) = J - Q * A^(-1/3)
  soit {J - Q * 40 ** (-1 / 3):.1f} MeV pour le calcium 40, {J - Q * 208 ** (-1 / 3):.1f} MeV pour le plomb 208.
  M1 n'a qu'un coefficient ({m1.coefs[3]:.1f} MeV) : une moyenne effective, proche
  de la valeur des noyaux lourds qui dominent l'echantillon. Physiquement, les neutrons en exces se logent en surface, ou ils
  coutent moins cher : c'est la PEAU DE NEUTRONS, mesuree sur le plomb
  208 (experience PREX) et directement liee a la physique des etoiles a
  neutrons.
""")
    def couches(z, n):
        xp, xn = valence(np.array([z])), valence(np.array([n]))
        S = xp + xn
        return float(sum(m3.coefs[m3.termes.index(k)] * v for k, v in
                         (("couches S", S), ("couches S^2", S ** 2),
                          ("couches p-n", xp * xn)))[0])
    print("  COUCHES. Valeur du terme de couches sur quelques noyaux :\n")
    print(f"  {'noyau':<10} {'Z':>4} {'N':>4} {'x_p':>6} {'x_n':>6} {'terme [MeV]':>12}")
    for nom, z, n in (("Pb-208", 82, 126), ("Sn-132", 50, 82), ("Ca-48", 20, 28),
                      ("Sm-152", 62, 90), ("Er-166", 68, 98), ("U-238", 92, 146)):
        xp, xn = float(valence(np.array([z]))[0]), float(valence(np.array([n]))[0])
        print(f"  {nom:<10} {z:>4} {n:>4} {xp:>6.2f} {xn:>6.2f} {couches(z, n):>12.2f}")
    print("""
  Les noyaux doublement magiques (Pb-208, Sn-132, Ca-48) sont a zero par
  construction ; ceux de milieu de couche (Er-166) perdent plus de 10 MeV.
  C'est l'ordre de grandeur des corrections de couches des grands modeles
  macroscopiques-microscopiques. Le terme p-n est POSITIF : l'interaction
  entre protons et neutrons de valence rend de la liaison, et c'est elle
  qui pousse les noyaux de milieu de couche a se deformer.
""")
    return m3


# ==========================================================================
def incertitudes(Z, N, B):
    titre("2 bis. CE QUE LES MASSES CONTRAIGNENT MAL")
    print("""
  Moindres carres B = X a + residu. Erreurs FORMELLES des coefficients :
      cov(a) = s^2 (X^T X)^-1,   s^2 = somme(residus^2) / (n - p)
  Elles supposent des residus independants ; ici les residus sont des
  erreurs de MODELE, correlees sur la carte : ces erreurs sont donc des
  bornes inferieures. Les correlations, elles, disent quelles combinaisons
  de termes les masses ne savent pas separer.
""")
    paires = (("volume", "surface"), ("surface", "Coulomb"), ("volume", "Coulomb"),
              ("asymetrie", "sym. surface"), ("Coulomb", "echange coul."))
    for nom, fab in (("M1", colonnes_M1), ("M2", colonnes_M2), ("M3", colonnes_M3)):
        m = ModeleLineaire(nom, fab).ajuster(Z, N, B)
        X = m._X(Z, N)
        n, p = X.shape
        r = B - X @ m.coefs
        s2 = r @ r / (n - p)
        cov = s2 * np.linalg.inv(X.T @ X)
        err = np.sqrt(np.diag(cov))
        corr = cov / np.outer(err, err)
        cond = np.linalg.cond(X / np.linalg.norm(X, axis=0))
        T = m.termes
        print(f"  {nom} : s = {np.sqrt(s2):.3f} MeV, conditionnement de X (colonnes normees) = {cond:.0f}")
        for k, c, e in zip(T, m.coefs, err):
            print(f"      {k:<14} {c:9.4f} +- {e:.4f}")
        for a, b in paires:
            if a in T and b in T:
                print(f"      correlation({a}, {b}) = {corr[T.index(a), T.index(b)]:+.3f}")
        # fissilite de l'U-236 et son erreur formelle (propagation lineaire)
        c = dict(zip(T, m.coefs))
        z, nn = 92, 144
        A, I = z + nn, nn - z
        Es = c["surface"] * A ** (2 / 3) - c.get("sym. surface", 0) * I ** 2 / A ** (4 / 3)
        Ec = c["Coulomb"] * z * (z - 1) / A ** (1 / 3) - c.get("echange coul.", 0) * z ** (4 / 3) / A ** (1 / 3)
        x = Ec / (2 * Es)
        g = np.zeros(p)
        g[T.index("Coulomb")] = z * (z - 1) / A ** (1 / 3) / (2 * Es)
        g[T.index("surface")] = -x / Es * A ** (2 / 3)
        if "echange coul." in T:
            g[T.index("echange coul.")] = -z ** (4 / 3) / A ** (1 / 3) / (2 * Es)
        if "sym. surface" in T:
            g[T.index("sym. surface")] = x / Es * I ** 2 / A ** (4 / 3)
        print(f"      fissilite de l'U-236 : x = {x:.3f} +- {np.sqrt(g @ cov @ g):.3f} (formelle)\n")
        if nom == "M1":
            r0 = 0.6 * 1.44 / c["Coulomb"]
            sig = c["surface"] / (4 * np.pi * r0 ** 2)
            hb2m, rho0 = 20.736, 0.16
            kF = (1.5 * np.pi ** 2 * rho0) ** (1 / 3)
            eF = hb2m * kF ** 2
            print(f"      a_C = (3/5) e^2 / r0       ->  r0 = {r0:.3f} fm")
            print(f"      a_S = 4 pi r0^2 sigma      ->  tension de surface sigma = {sig:.2f} MeV/fm^2")
            print(f"      gaz de Fermi, rho0 = {rho0} fm^-3 : k_F = {kF:.3f} fm^-1, e_F = {eF:.1f} MeV,")
            print(f"      part cinetique de l'asymetrie e_F/3 = {eF / 3:.1f} MeV (ajuste : {c['asymetrie']:.1f})")
            ex = 0.75 * (3 / (2 * np.pi)) ** (2 / 3) * 1.44 / r0
            print(f"      echange coulombien de Slater (3/4)(3/2pi)^(2/3) e^2/r0 = {ex:.2f} MeV\n")
    print("""  Surface et Coulomb sont correles a ~0,95 : les masses fixent leur
  COMBINAISON, pas chacun separement. Or la barriere de fission depend de
  leur RAPPORT, x = E_c / 2 E_s (lecon 6). D'un modele a l'autre, x(U-236)
  bouge bien plus que son erreur formelle : c'est une erreur de modele.
""")


# ==========================================================================
def distance(t, Z, N, B, K=8):
    titre("3. COMMENT L'ERREUR CROIT AVEC LA DISTANCE AU CONNU")
    print(f"""
  Experience : pour chaque element ayant au moins {K + 6} isotopes mesures,
  on retire les {K} plus riches en neutrons. On ajuste sur le reste, puis on
  predit les noyaux retires. On note d = nombre de neutrons au-dela du
  dernier isotope connu de l'element. C'est exactement la situation reelle :
  predire la masse de noyaux de plus en plus exotiques, vers la limite
  d'existence des noyaux (drip line) ou se joue le processus r.
""")
    test = np.zeros(len(t), bool)
    dist = np.zeros(len(t), int)
    for z in np.unique(t["Z"]):
        idx = np.where(t["Z"] == z)[0]
        if len(idx) < K + 6:
            continue
        o = idx[np.argsort(N[idx])]
        test[o[-K:]] = True
        dernier = N[o[-K - 1]]
        dist[o[-K:]] = (N[o[-K:]] - dernier).astype(int)

    modeles = [ModeleLineaire("M2", colonnes_M2), ModeleLineaire("M3", colonnes_M3),
               ModeleCorrige()]
    courbes = {}
    for m in modeles:
        m.ajuster(Z[~test], N[~test], B[~test])
        e = B[test] - m.predire(Z[test], N[test])
        courbes[m.nom[:2]] = [(d, rms(e[dist[test] == d]), int((dist[test] == d).sum()))
                              for d in range(1, K + 1) if (dist[test] == d).sum() >= 5]
    print(f"  {test.sum()} noyaux predits, sur {len(np.unique(t['Z'][test]))} elements\n")
    print(f"  {'d':>3} {'noyaux':>7} {'M2':>8} {'M3':>8} {'M4':>8}")
    for i, (d, _, n) in enumerate(courbes["M3"]):
        print(f"  {d:>3} {n:>7} {courbes['M2'][i][1]:>8.3f} "
              f"{courbes['M3'][i][1]:>8.3f} {courbes['M4'][i][1]:>8.3f}")
    print("""
  M4 est imbattable a un ou deux pas du connu (4 fois meilleur que M3),
  puis son erreur croit bien plus vite : x3,5 de d=1 a d=8, contre x1,4
  pour M3. L'avantage fond de 4x a 1,6x sur la plage testee.

  C'est mecanique : sa correction est une moyenne ponderee des erreurs des
  voisins, avec un poids en exp(-d^2 / 2 l^2). Plus on s'eloigne, moins il
  y a de voisins, plus la correction s'efface et plus il ne reste que la
  physique de M3. (Elle ne s'efface pas tout a fait a d=8 : les elements
  voisins, Z-1 et Z+1, fournissent encore un peu d'information en
  diagonale sur la carte.)

  LA LECON GENERALE. Un modele statistique ne sait rien au-dela de ses
  donnees ; seule la physique extrapole. C'est pourquoi les tables de
  masse utilisees en astrophysique reposent sur des modeles physiques
  (FRDM, HFB, WS4, Duflo-Zuker) et que les corrections par apprentissage
  sont reservees au voisinage des noyaux mesures.
""")
    return courbes


# ==========================================================================
def litterature(resultats):
    titre("4. OU SE SITUENT CES MODELES ?")
    print("""
  Ecart RMS sur les masses d'AME2020 de modeles publies (valeurs tirees de
  la litterature recente, voir MODELES_MASSE.md pour les references) :

    Gogny D1M (HFB, interaction de portee finie)          0,81 MeV
    Duflo-Zuker DZ10 (10 parametres)                      0,56 MeV
    HFB-27 (Skyrme, Goriely et al.)                       0,52 MeV
    WS4 (Weizsacker-Skyrme, Wang et al.)                  ~0,3 MeV
    reseaux de neurones sur le residu                     ~0,2 MeV
""")
    print("  Nos modeles, en interpolation :\n")
    for nom, _, r_cv, _ in resultats:
        print(f"    {nom:<44} {r_cv:.2f} MeV")
    print("""
  M3, 11 parametres et 30 lignes de code, fait 1,2 MeV : deux fois l'erreur
  de Duflo-Zuker, avec une physique bien plus fruste. Pour descendre
  encore, il faut ce que M3 n'a pas : la DEFORMATION des noyaux (la
  plupart des noyaux de milieu de couche sont en ballon de rugby, et
  l'energie de deformation se chiffre en MeV) et un vrai calcul des
  niveaux individuels au lieu d'un simple comptage de valence.
""")


# ==========================================================================
def graphiques(t, Z, N, B, resultats, courbes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(14, 13))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.15, 1], hspace=0.42, wspace=0.22)

    # (a) les trois epreuves
    ax = fig.add_subplot(gs[0, 0])
    noms = [r[0][:2] for r in resultats]
    x = np.arange(len(noms))
    for i, (lab, col) in enumerate((("ajustement", "#9db4cf"),
                                    ("interpolation", "#1f4e79"),
                                    ("extrapolation", "#e76f51"))):
        ax.bar(x + (i - 1) * 0.27, [r[i + 1] for r in resultats], 0.27,
               label=lab, color=col)
    ax.set_xticks(x)
    ax.set_xticklabels(noms)
    ax.set_ylabel("ecart RMS sur B  [MeV]")
    ax.set_title("Trois epreuves")
    ax.axhline(0.3, color="seagreen", ls="--", lw=1)
    ax.text(-0.45, 0.42, "WS4 ~0,3 MeV (meilleur modele global)", color="seagreen", fontsize=8.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, axis="y")

    # (b) erreur vs distance
    ax = fig.add_subplot(gs[0, 1])
    for cle, col, lab in (("M2", "#999999", "M2 goutte etendue"),
                          ("M3", "#1f4e79", "M3 + couches"),
                          ("M4", "#e76f51", "M4 + regression a noyau")):
        d = [c[0] for c in courbes[cle]]
        r = [c[1] for c in courbes[cle]]
        ax.plot(d, r, "o-", color=col, lw=2, label=lab)
    ax.set_xlabel("neutrons au-dela du dernier isotope connu")
    ax.set_ylabel("ecart RMS sur B  [MeV]")
    ax.set_title("Extrapolation vers les noyaux riches en neutrons")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)

    # (c) residus en fonction de N, M2 / M3 / M4
    modeles = [ModeleLineaire("M2", colonnes_M2), ModeleLineaire("M3", colonnes_M3),
               ModeleCorrige()]
    titres = ("M2 goutte etendue : les couches dominent le residu",
              "M3 + couches : les grandes bosses ont disparu",
              "M4 + regression a noyau (ajustement)")
    sub = gs[1, :].subgridspec(1, 3, wspace=0.08)
    axes = [fig.add_subplot(sub[0, i]) for i in range(3)]
    residus = []
    for ax, m, tt in zip(axes, modeles, titres):
        m.ajuster(Z, N, B)
        r = B - m.predire(Z, N)
        residus.append((m, r))
        ax.scatter(N, r, c=Z, s=3, cmap="viridis")
        for nm in (28, 50, 82, 126):
            ax.axvline(nm, color="crimson", ls=":", lw=1)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylim(-9, 9)
        ax.set_title(f"{tt}\nRMS = {rms(r):.2f} MeV", fontsize=9.5)
        ax.set_xlabel("N")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("mesure - modele  [MeV]")
    for ax in axes[1:]:
        ax.set_yticklabels([])

    # (d) carte du residu de M3
    ax = fig.add_subplot(gs[2, 0])
    r3 = residus[1][1]
    sc = ax.scatter(N, Z, c=r3, s=5, cmap="RdBu_r", vmin=-4, vmax=4)
    for nm in (28, 50, 82, 126):
        ax.axvline(nm, color="k", ls=":", lw=0.6)
    for zm in (28, 50, 82):
        ax.axhline(zm, color="k", ls=":", lw=0.6)
    ax.set_xlabel("N")
    ax.set_ylabel("Z")
    ax.set_title("Residu de M3 sur la carte : ce qui reste a expliquer")
    fig.colorbar(sc, ax=ax, label="MeV")

    # (e) carte de la correction apprise par M4
    ax = fig.add_subplot(gs[2, 1])
    m4 = residus[2][0]
    corr = m4.predire(Z, N) - m4.base.predire(Z, N)
    sc = ax.scatter(N, Z, c=corr, s=5, cmap="RdBu_r", vmin=-4, vmax=4)
    for nm in (28, 50, 82, 126):
        ax.axvline(nm, color="k", ls=":", lw=0.6)
    for zm in (28, 50, 82):
        ax.axhline(zm, color="k", ls=":", lw=0.6)
    ax.set_xlabel("N")
    ax.set_title(f"Correction apprise par M4 (l = {m4.l:g} nucleons)")
    fig.colorbar(sc, ax=ax, label="MeV")

    fig.text(0.5, 0.005, "Donnees : AME2020, masses mesurees seulement "
             f"({len(t)} noyaux, Z,N >= 8) -- M. Wang et al., Chin. Phys. C 45 (2021) 030003",
             ha="center", fontsize=8.5, color="#444444")
    fig.savefig("modeles_masse.png", dpi=125, bbox_inches="tight")
    print("\n  --> modeles_masse.png ecrit")
    return m4


# ==========================================================================
def main():
    t, Z, N, B = charger()
    resultats = epreuves(t, Z, N, B)
    coefficients(t, Z, N, B)
    incertitudes(Z, N, B)
    courbes = distance(t, Z, N, B)
    litterature(resultats)
    graphiques(t, Z, N, B, resultats, courbes)


if __name__ == "__main__":
    main()
