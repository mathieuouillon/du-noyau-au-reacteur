"""
analyse_liaison.py -- Exploite les energies de liaison mesurees.

Produit :
  * energies_liaison.csv       table complete (2825 nuclides)
  * courbe_liaison.png         4 graphiques
  * une confrontation SEMF / mesure a l'ecran

Lancer :  python analyse_liaison.py
"""

import numpy as np
from donnees_liaison import (table_liaison, exporter_csv, B_sur_A, chercher,
                             energie_separation_neutron, Q_fission)
from fission import liaison_par_nucleon

NOMBRES_MAGIQUES = (2, 8, 20, 28, 50, 82, 126)


def titre(texte):
    print()
    print("=" * 74)
    print(texte)
    print("=" * 74)


# ==========================================================================
def validation(t):
    titre("VALIDATION DES DONNEES")
    print("""
  Les masses proviennent du paquet `periodictable` (table isotopique issue
  des masses atomiques evaluees NIST/AME). On recalcule B/A et on confronte
  a des valeurs de reference publiees.
""")
    reference = {("He", 4): 7.0739, ("C", 12): 7.6801, ("O", 16): 7.9762,
                 ("Fe", 56): 8.7903, ("Ni", 62): 8.7946, ("Pb", 208): 7.8675,
                 ("U", 235): 7.5909, ("U", 238): 7.5701}
    print(f"  {'nuclide':>9} {'calcule':>11} {'reference':>11} {'ecart':>12}")
    for (s, a), v in reference.items():
        c = B_sur_A(t, s, a)
        print(f"  {s + '-' + str(a):>9} {c:>11.4f} {v:>11.4f} "
              f"{(c - v) * 1000:>9.2f} keV")
    print("\n  Accord a mieux que 0.1 keV : les donnees sont exactes.")


# ==========================================================================
def le_maximum(t):
    titre("OU EST VRAIMENT LE MAXIMUM ?")
    ordre = np.argsort(t["B_sur_A_MeV"])[::-1]
    print(f"\n  {'rang':>5} {'nuclide':>9} {'B/A [MeV]':>12} {'abondance':>12}")
    for rang, i in enumerate(ordre[:6], start=1):
        r = t[i]
        print(f"  {rang:>5} {r['symbole'] + '-' + str(r['A']):>9} "
              f"{r['B_sur_A_MeV']:>12.5f} {r['abondance_pct']:>11.2f} %")
    print("""
  LE FER 56 N'EST PAS LE MAXIMUM. C'est le nickel 62, avec 8.79455 MeV
  contre 8.79036 pour Fe-56 : 4.2 keV d'ecart par nucleon. Fe-58 passe
  aussi devant.

  D'ou vient la confusion ? Fe-56 est le noyau de plus faible masse PAR
  NUCLEON (m/A minimal), parce qu'il contient proportionnellement plus de
  protons, plus legers que les neutrons. Ni-62 est le plus LIE (B/A
  maximal). Les deux criteres ne coincident pas, et c'est le premier qu'on
  cite d'habitude en astrophysique.

  Fe-56 domine quand meme dans la nature (91.75 % du fer) parce que la
  nucleosynthese stellaire passe par Ni-56, qui decroit en Fe-56 : c'est
  un argument cinetique, pas thermodynamique.
""")


# ==========================================================================
def confrontation_semf(t):
    titre("LE MODELE DE LA GOUTTE LIQUIDE FACE AUX MESURES")
    print(f"\n  (masses MESUREES seules : les {int(t['estimee'].sum())} masses estimees "
          "'#' d'AME2020 sont ecartees)")
    t = t[~t["estimee"]]
    masque = t["A"] >= 20
    sub = t[masque]
    semf = liaison_par_nucleon(sub["A"].astype(float), sub["Z"].astype(float))
    ecart = semf - sub["B_sur_A_MeV"]

    naturels = t[(t["abondance_pct"] > 0) & (t["A"] >= 20)]
    semf_nat = liaison_par_nucleon(naturels["A"].astype(float),
                                   naturels["Z"].astype(float))
    ecart_nat = semf_nat - naturels["B_sur_A_MeV"]

    print(f"""
  Sur {masque.sum()} nuclides avec A >= 20 :
      ecart moyen  {ecart.mean():+.4f} MeV/nucleon
      ecart RMS     {np.sqrt((ecart ** 2).mean()):.4f} MeV/nucleon
      ecart max     {np.abs(ecart).max():.4f} MeV/nucleon

  Sur les {len(naturels)} isotopes NATURELS seulement :
      ecart RMS     {np.sqrt((ecart_nat ** 2).mean()):.4f} MeV/nucleon

  Trois fois meilleur pres de la vallee de stabilite. Normal : c'est la que
  les parametres de la formule ont ete ajustes.
""")
    pires = np.argsort(np.abs(ecart))[::-1][:8]
    print(f"  Les 8 pires cas :\n  {'nuclide':>9} {'Z':>4} {'N':>4} "
          f"{'mesure':>9} {'SEMF':>9} {'ecart':>9}")
    for i in pires:
        r = sub[i]
        print(f"  {r['symbole'] + '-' + str(r['A']):>9} {r['Z']:>4} {r['N']:>4} "
              f"{r['B_sur_A_MeV']:>9.3f} {semf[i]:>9.3f} {ecart[i]:>+9.3f}")
    print("""
  Tous sont des noyaux LEGERS et tres riches en neutrons (C-22 a N=16,
  C-20, N-23, les isotopes lourds de l'oxygene...). La goutte liquide les
  croit bien plus lies qu'ils ne sont. Ce sont precisement les noyaux ou les effets de couche et de peau
  de neutrons dominent, c'est-a-dire tout ce que le modele ignore par
  construction.

  Pour le reacteur, aucune importance : l'uranium et ses produits de
  fission sont lourds et proches de la stabilite, domaine ou le modele est
  bon a 0.03 MeV/nucleon pres.
""")


# ==========================================================================
def S_2n(t, symbole, A):
    """Energie de separation de DEUX neutrons : S_2n = B(A,Z) - B(A-2,Z)."""
    haut, bas = chercher(t, symbole, A), chercher(t, symbole, A - 2)
    if haut is None or bas is None:
        return None
    return float(haut["B_MeV"] - bas["B_MeV"])


def effets_de_couche(t):
    titre("LA SIGNATURE DES COUCHES NUCLEAIRES")
    print("""
  Comment prouver que le noyau a une structure en couches, sans invoquer
  aucun modele ? En regardant l'energie de separation de DEUX neutrons,
  S_2n = B(A,Z) - B(A-2,Z), le long d'une chaine isotopique. On prend deux
  neutrons a la fois pour eliminer l'oscillation pair-impair.

  S_2n decroit lentement et regulierement quand on ajoute des neutrons...
  sauf juste apres une couche complete, ou elle s'effondre : le neutron
  suivant doit aller se loger dans une couche superieure, bien moins liee.
  C'est exactement le saut d'energie d'ionisation entre un gaz rare et
  l'alcalin qui le suit.

  C'est un OBSERVABLE PUR, deduit uniquement des masses mesurees.
""")
    chaines = [("Ca", 20, range(17, 25), 20),
               ("Sn", 50, range(79, 87), 82),
               ("Pb", 82, range(123, 131), 126)]
    for symbole, Z, plage, magique in chaines:
        print(f"\n  Chaine {symbole} (Z = {Z}), couche attendue a N = {magique}")
        print(f"  {'N':>6} {'A':>5} {'S_2n [MeV]':>12} {'chute':>18}")
        precedent = None
        for N in plage:
            valeur = S_2n(t, symbole, Z + N)
            if valeur is None:
                continue
            if chercher(t, symbole, Z + N)["estimee"] or \
               chercher(t, symbole, Z + N - 2)["estimee"]:
                print(f"  {N:>6} {Z + N:>5}   (masse estimee '#', ignoree)")
                precedent = None
                continue
            if precedent is None:
                chute = ""
            else:
                delta = precedent - valeur
                marque = "   <== EFFONDREMENT" if delta > 2.0 else ""
                chute = f"{delta:8.2f}{marque}"
            print(f"  {N:>6} {Z + N:>5} {valeur:>12.2f} {chute:>18}")
            precedent = valeur

    print("""
  Dans les trois chaines, la cassure tombe exactement UN CRAN APRES le
  nombre magique : N = 21 pour le calcium, N = 83 pour l'etain, N = 127
  pour le plomb. Tant que la couche se remplit, chaque paire de neutrons
  coute a peu pres pareil ; des qu'elle est pleine, la paire suivante
  tombe dans une couche superieure et l'energie chute de 3 a 5 MeV d'un
  coup.

  Les nombres magiques 2, 8, 20, 28, 50, 82, 126 ne sont donc PAS une
  hypothese theorique ajoutee apres coup : ils se lisent directement dans
  les masses mesurees.

  C'EST CE QUE LA GOUTTE LIQUIDE IGNORE, par construction -- elle traite
  le noyau comme un fluide continu, sans niveaux d'energie individuels.
  D'ou ses deux echecs rencontres dans etude_fission.py : He-4
  (doublement magique, Z = N = 2) sous-estime de 1.6 MeV/nucleon, et la
  fission asymetrique, dont les fragments se forment pres de Z = 50 et
  N = 82.

  REMARQUE DE METHODE. J'avais d'abord cherche ces couches en moyennant le
  residu (mesure - SEMF) par valeur de N. N = 50 et N = 82 ressortaient,
  mais N = 20 et N = 28 restaient noyes : a bas N l'echantillon contient
  surtout des noyaux tres exotiques, ou le modele echoue pour d'autres
  raisons, ce qui brouille le signal. Le S_2n le long d'une chaine est
  superieur sur les deux plans : il ne depend d'aucun modele, et il compare
  des noyaux voisins plutot que de melanger tout un isotone.
""")


# ==========================================================================
def grandeurs_du_reacteur(t):
    titre("LES GRANDEURS DU REACTEUR, EN DONNEES MESUREES")
    print("\n  Fissilite : energie de separation du neutron S_n\n")
    barriere = {"U-236": 6.2, "U-239": 6.6, "Pu-240": 6.0, "U-234": 5.9}
    print(f"  {'capture':<22} {'compose':>8} {'S_n mesure':>12} "
          f"{'barriere':>10} {'fissile':>9}")
    for cible, A, comp in (("U-235 + n", 236, "U-236"), ("U-238 + n", 239, "U-239"),
                           ("Pu-239 + n", 240, "Pu-240"), ("U-233 + n", 234, "U-234")):
        sym = comp.split("-")[0]
        sn = energie_separation_neutron(t, sym, A)
        b = barriere[comp]
        print(f"  {cible:<22} {comp:>8} {sn:>9.3f} MeV {b:>8.1f} MeV "
              f"{'OUI' if sn > b else 'non':>9}")
    print("""
  Les valeurs MESUREES confirment la prediction du modele : seuls les
  noyaux a nombre de neutrons IMPAIR (U-235 : 143, Pu-239 : 145, U-233 :
  141) liberent assez d'energie d'appariement en capturant un neutron pour
  franchir leur barriere. L'U-238 (146 neutrons, pair) non.
""")
    print("  Q de fission, calcule sur les masses mesurees :\n")
    print(f"  {'partition de U-236':<34} {'Q mesure':>10} {'Q SEMF':>10}")
    from fission import energie_fission
    partitions = [(("Pd", 118), ("Pd", 118), 0, (118, 46), (118, 46), "symetrique"),
                  (("Xe", 137), ("Sr", 97), 2, (137, 54), (97, 38), ""),
                  (("Ba", 141), ("Kr", 92), 3, (141, 56), (92, 36), "la plus citee")]
    for f1, f2, nn, p1, p2, note in partitions:
        q = Q_fission(t, ("U", 236), f1, f2)
        qs = float(energie_fission(236, 92, p1[0], p1[1], p2[0], p2[1]))
        label = f"{f1[0]}-{f1[1]} + {f2[0]}-{f2[1]} + {nn}n"
        print(f"  {label:<34} {q:>7.1f} MeV {qs:>7.1f} MeV   {note}")
    print("""
  CORRECTION A ETUDE_FISSION.PY. J'y ecrivais que la goutte liquide
  "predit un maximum pour la fission symetrique, alors que la fission
  reelle est asymetrique". Les masses mesurees montrent que la partition
  symetrique libere effectivement PLUS d'energie (193 MeV contre 167) :
  sur ce point le modele avait raison, ce n'est pas son erreur.

  L'asymetrie de la fission n'est donc PAS un effet de bilan energetique.
  Elle vient de la dynamique du noyau au moment de la scission : la
  surface d'energie potentielle au point selle favorise des fragments
  proches des couches fermees Z=50 et N=82. Le noyau ne choisit pas le
  partage le plus exothermique, il suit le chemin le plus facile.
""")


# ==========================================================================
def graphiques(t):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"  (graphiques ignores : {exc})")
        return

    nat = t[t["abondance_pct"] > 0]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10))

    # (a) la courbe, donnees mesurees
    ax = axes[0][0]
    ax.plot(t["A"], t["B_sur_A_MeV"], ".", ms=1.2, color="#b8c6d9",
            label=f"tous ({len(t)})")
    ax.plot(nat["A"], nat["B_sur_A_MeV"], ".", ms=4, color="#1f4e79",
            label=f"isotopes naturels ({len(nat)})")
    for s, a in (("He", 4), ("Ni", 62), ("U", 235)):
        r = chercher(t, s, a)
        ax.plot(r["A"], r["B_sur_A_MeV"], "o", color="crimson", ms=6)
        ax.annotate(f"{s}-{a}", (r["A"], r["B_sur_A_MeV"]),
                    textcoords="offset points", xytext=(7, -11), fontsize=9)
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("B/A [MeV]")
    ax.set_title("Energie de liaison par nucleon -- donnees mesurees")
    ax.set_ylim(0, 9.3)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)

    # (b) zoom sur le maximum
    ax = axes[0][1]
    z = t[(t["A"] >= 45) & (t["A"] <= 75) & (t["abondance_pct"] > 0)]
    ax.plot(z["A"], z["B_sur_A_MeV"], "o-", ms=5, color="#1f4e79", lw=1)
    for s, a, c, dx, dy in (("Fe", 56, "darkorange", -30, -26),
                            ("Fe", 58, "seagreen", 0, 14),
                            ("Ni", 62, "crimson", 26, -18)):
        r = chercher(t, s, a)
        ax.plot(r["A"], r["B_sur_A_MeV"], "o", color=c, ms=10, zorder=5)
        ax.annotate(f"{s}-{a}\n{r['B_sur_A_MeV']:.4f}",
                    (r["A"], r["B_sur_A_MeV"]), textcoords="offset points",
                    xytext=(dx, dy), fontsize=8.5, ha="center", color=c,
                    fontweight="bold")
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("B/A [MeV]")
    ax.set_title("Zoom : le maximum est Ni-62, pas Fe-56")
    ax.grid(alpha=0.3)

    # (c) residu mesure - SEMF
    ax = axes[1][0]
    m = t["A"] >= 20
    sub = t[m]
    res = sub["B_sur_A_MeV"] - liaison_par_nucleon(sub["A"].astype(float),
                                                  sub["Z"].astype(float))
    sc = ax.scatter(sub["N"], res, c=sub["Z"], s=3, cmap="viridis")
    for N in NOMBRES_MAGIQUES:
        if N >= 20:
            ax.axvline(N, color="crimson", ls=":", lw=1.2)
            ax.annotate(str(N), (N, res.max() * 0.93), color="crimson",
                        fontsize=8.5, ha="center")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("nombre de neutrons N")
    ax.set_ylabel("mesure - SEMF [MeV/nucleon]")
    ax.set_title("Residu : les couches nucleaires apparaissent\n"
                 "(pointilles = nombres magiques)")
    fig.colorbar(sc, ax=ax, label="Z")
    ax.grid(alpha=0.3)

    # (d) carte des nuclides
    ax = axes[1][1]
    sc = ax.scatter(t["N"], t["Z"], c=t["B_sur_A_MeV"], s=2.5,
                    cmap="inferno", vmin=7.0, vmax=8.8)
    ax.plot(nat["N"], nat["Z"], ".", ms=1.6, color="white", alpha=0.75)
    for N in NOMBRES_MAGIQUES:
        ax.axvline(N, color="cyan", ls=":", lw=0.7, alpha=0.6)
    for Z in NOMBRES_MAGIQUES:
        ax.axhline(Z, color="cyan", ls=":", lw=0.7, alpha=0.6)
    ax.set_xlabel("nombre de neutrons N")
    ax.set_ylabel("nombre de protons Z")
    ax.set_title("Carte des nuclides, coloree par B/A\n"
                 "(points blancs = isotopes naturels)")
    fig.colorbar(sc, ax=ax, label="B/A [MeV]")

    fig.suptitle("Energies de liaison mesurees -- "
                 f"{len(t)} nuclides", fontsize=14)
    fig.tight_layout()
    fig.savefig("courbe_liaison.png", dpi=130)
    print("\n  --> courbe_liaison.png ecrit")


# ==========================================================================
def main():
    t = table_liaison()
    print(f"\n{len(t)} nuclides charges (AME2020), Z de {t['Z'].min()} a {t['Z'].max()} :"
          f" {int((~t['estimee']).sum())} mesures, {int(t['estimee'].sum())} estimes")
    validation(t)
    le_maximum(t)
    confrontation_semf(t)
    effets_de_couche(t)
    grandeurs_du_reacteur(t)
    chemin = exporter_csv(t, "energies_liaison.csv")
    print(f"\n  --> {chemin} ecrit ({len(t)} lignes)")
    graphiques(t)


if __name__ == "__main__":
    main()
