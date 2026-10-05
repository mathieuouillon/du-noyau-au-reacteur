"""
trace_termes.py -- Chaque terme des modeles de masse, montre separement.

On ajuste M3 (le plus complet des modeles lineaires, modeles_masse.py) sur les
masses mesurees d'AME2020, puis on regarde la contribution de chacun de ses
11 termes a l'energie de liaison :

  le long de la vallee de stabilite, une figure par panneau :
      termes_goutte.png   la goutte : volume, surface, Coulomb, asymetrie,
                          en MeV par nucleon ;
      termes_m2.png       les corrections de M2 et l'appariement, en MeV ;
      termes_couches.png  les trois termes de couches de M3, en MeV ;
      termes_gain.png     ce que rapporte chaque terme : l'ecart RMS quand on
                          les ajoute un a un.
  * termes_carte.png
      six termes sur la carte (N, Z) : ou chacun agit.

Convention : contribution a B (B > 0 pour un noyau lie). Positif = le terme
LIE davantage le noyau ; negatif = il coute de la liaison.

Lancer :  python trace_termes.py   ->  termes_goutte.png, termes_m2.png,
                                       termes_couches.png, termes_gain.png,
                                       termes_carte.png
"""

import numpy as np
from donnees_liaison import table_liaison
from modeles_masse import ModeleLineaire, colonnes_M3

# Palette de reference (voir la page Code) : surface, encres, series
SURFACE = "#fcfcfb"
ENCRE = "#0b0b0b"
ENCRE_2 = "#52514e"
DISCRET = "#898781"
GRILLE = "#e1e0d9"
AXE = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
BLEUS = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ROUGE, GRIS_MILIEU = "#e34948", "#f0efec"

NOMS = {
    "volume": "volume", "surface": "surface", "Coulomb": "Coulomb",
    "asymetrie": "asymétrie", "appariement": "appariement",
    "sym. surface": "symétrie de surface", "echange coul.": "échange coulombien",
    "Wigner": "Wigner", "couches S": "couches, a₁·S", "couches S^2": "couches, a₂·S²",
    "couches p-n": "couches, a₃·xₚxₙ",
}
SYMBOLES = {8: "O", 20: "Ca", 28: "Ni", 40: "Zr", 50: "Sn", 56: "Ba", 82: "Pb"}


def nombre(x):
    """1 decimale, virgule francaise, vrai signe moins, pas de -0,0."""
    if abs(x) < 0.05:
        return "0"
    return f"{x:.1f}".replace(".", ",").replace("-", "−")


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


def charger():
    t = table_liaison(mesurees_seulement=True)
    t = t[(t["Z"] >= 8) & (t["N"] >= 8)]
    return t["Z"].astype(float), t["N"].astype(float), t["B_MeV"]


def contributions(Z, N, modele):
    """{terme: contribution a B en MeV}, pour chaque noyau."""
    cols = colonnes_M3(Z, N)
    return {k: modele.coefs[i] * cols[k] for i, k in enumerate(modele.termes)}


def vallee(Z, N, B):
    """Pour chaque A, l'isobare le plus lie parmi les noyaux mesures."""
    A = Z + N
    idx = []
    for a in np.unique(A):
        i = np.where(A == a)[0]
        idx.append(i[np.argmax(B[i])])
    return np.array(idx)


def ablation(Z, N, B, termes):
    """Ecart RMS quand on ajoute les termes de M3 un a un, dans l'ordre."""
    cols = colonnes_M3(Z, N)
    res = []
    for k in range(1, len(termes) + 1):
        X = np.column_stack([cols[t] for t in termes[:k]])
        a, *_ = np.linalg.lstsq(X, B, rcond=None)
        res.append(rms(B - X @ a))
    return res


# ==========================================================================
def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXE)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=DISCRET, labelcolor=ENCRE_2, labelsize=11, width=0.8)
    ax.grid(True, color=GRILLE, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(ENCRE_2)
    ax.yaxis.label.set_color(ENCRE_2)
    ax.xaxis.label.set_fontsize(11.5)
    ax.yaxis.label.set_fontsize(11.5)


def titre(ax, texte, sous_titre=None):
    ax.set_title(texte, loc="left", fontsize=13.5, color=ENCRE, fontweight="semibold",
                 pad=24 if sous_titre else 8)
    if sous_titre:
        ax.text(0, 1.015, sous_titre, transform=ax.transAxes, fontsize=10.5,
                color=ENCRE_2, va="bottom")


def ecarter(ys, ecart):
    """Positions d'etiquettes : dans l'ordre des valeurs, espacees d'au moins
    `ecart` (unites de donnees), le plus pres possible des valeurs."""
    ordre = np.argsort(ys)
    pos = np.array(ys, float)[ordre]
    for _ in range(200):
        bouge = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < ecart:
                milieu = 0.5 * (pos[i] + pos[i - 1])
                pos[i - 1], pos[i] = milieu - ecart / 2, milieu + ecart / 2
                bouge = True
        if not bouge:
            break
    out = np.empty(len(ys))
    out[ordre] = pos
    return out


def etiquettes(ax, xfin, yfin, noms, couleurs, ecart):
    """Etiquette directe en bout de courbe : trait de rappel fin, segment de
    couleur (l'identite), texte en encre (jamais dans la couleur de la serie)."""
    yl = ecarter(yfin, ecart)
    x0 = xfin
    dx = 0.04 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    for y0, y1, nom, col in zip(yfin, yl, noms, couleurs):
        ax.plot([x0, x0 + dx], [y0, y1], color=DISCRET, lw=0.6, clip_on=False)
        ax.plot([x0 + dx, x0 + 2.2 * dx], [y1, y1], color=col, lw=2.4,
                solid_capstyle="round", clip_on=False)
        ax.text(x0 + 2.6 * dx, y1, nom, va="center", fontsize=10.5, color=ENCRE)


# ==========================================================================
def figure_vallee(Z, N, B, m3, contrib):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter
    import periodictable

    iv = vallee(Z, N, B)
    A = (Z + N)[iv]
    o = np.argsort(A)
    iv = iv[o]
    A, Zv, Nv = A[o], Z[iv], N[iv]
    c = {k: v[iv] for k, v in contrib.items()}
    xlim = (A.min() - 2, A.max() + 2)

    figures = {}

    def panneau(fichier, gauche=0.09, droite=0.80):
        """Une figure par panneau : chacune s'affiche en pleine largeur, a cote
        du paragraphe qui la commente."""
        fig, ax = plt.subplots(figsize=(9, 5.6), facecolor=SURFACE)
        fig.subplots_adjust(left=gauche, right=droite, top=0.83, bottom=0.16)
        figures[fichier] = fig
        return ax

    # (a) la goutte, par nucleon
    ax = panneau("termes_goutte.png", gauche=0.10, droite=0.80)
    style(ax)
    cles = ("volume", "surface", "Coulomb", "asymetrie")
    for i, k in enumerate(cles):
        ax.plot(A, c[k] / A, color=SERIES[i], lw=2, solid_joinstyle="round")
    tot = sum(c[k] for k in c) / A
    ax.plot(A, tot, color=ENCRE_2, lw=1.2)
    ax.text(A[len(A) // 3], tot[len(A) // 3] + 0.9, "B/A total (M3)", fontsize=10.5, color=ENCRE_2)
    ax.axhline(0, color=AXE, lw=0.8)
    ax.set_xlim(*xlim)
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("contribution à B/A (MeV par nucléon)")
    etiquettes(ax, xlim[1], [c[k][-1] / A[-1] for k in cles], [NOMS[k] for k in cles],
               SERIES, ecart=1.1)
    titre(ax, "(a) La goutte : quatre termes qui se disputent B/A",
          "Le volume lie (+15,5 MeV par nucléon) ; surface, Coulomb et asymétrie coûtent.")

    # (b) corrections de M2 + appariement, en MeV
    ax = panneau("termes_m2.png", gauche=0.10, droite=0.70)
    style(ax)
    cles = ("echange coul.", "sym. surface", "Wigner", "appariement")
    for i, k in enumerate(cles[:3]):
        ax.plot(A, c[k], color=SERIES[i], lw=2, solid_joinstyle="round")
    ax.plot(A, c["appariement"], "o", ms=3.4, color=SERIES[3], mec=SURFACE, mew=0.6)
    ax.axhline(0, color=AXE, lw=0.8)
    ax.set_xlim(*xlim)
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("contribution à B (MeV)")
    fins = [c[k][-1] for k in cles[:3]] + [float(np.max(c["appariement"][-15:]))]
    etiquettes(ax, xlim[1], fins, [NOMS[k] + (" (points)" if k == "appariement" else "")
                                   for k in cles], SERIES, ecart=6.5)
    titre(ax, "(b) Les corrections de M2, et l'appariement",
          "Échange et symétrie de surface rendent de la liaison ; Wigner et appariement restent petits.")

    # (c) les couches
    ax = panneau("termes_couches.png", gauche=0.10, droite=0.82)
    style(ax)
    tot_c = c["couches S"] + c["couches S^2"] + c["couches p-n"]
    series = [("a₁·S", c["couches S"]), ("a₂·S²", c["couches S^2"]),
              ("a₃·xₚxₙ", c["couches p-n"]), ("somme", tot_c)]
    for i, (nom, y) in enumerate(series):
        ax.plot(A, y, color=SERIES[i], lw=2 if i < 3 else 2.6, solid_joinstyle="round")
    ax.axhline(0, color=AXE, lw=0.8)
    bas = min(float(y.min()) for _, y in series)
    haut = max(float(y.max()) for _, y in series)
    ax.set_ylim(bas - 2, haut + 6)
    ax.set_xlim(*xlim)
    # les noyaux ou la somme remonte pres de zero : doublement magiques
    for i in range(1, len(A) - 1):
        if tot_c[i] > -1.5 and tot_c[i] >= tot_c[i - 1] and tot_c[i] >= tot_c[i + 1] and A[i] > 10:
            nom = f"{int(A[i])}{periodictable.elements[int(Zv[i])].symbol}"
            ax.annotate(nom, (A[i], tot_c[i]), xytext=(A[i], haut + 3.5), ha="center",
                        fontsize=10, color=ENCRE_2,
                        arrowprops=dict(arrowstyle="-", color=DISCRET, lw=0.6))
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("contribution à B (MeV)")
    etiquettes(ax, xlim[1], [y[-1] for _, y in series], [n for n, _ in series], SERIES,
               ecart=0.07 * (haut - bas + 8))
    titre(ax, "(c) Les couches de M3 : nulles aux noyaux doublement magiques",
          "a₁·S pénalise l'éloignement des couches ; a₂·S² et a₃·xₚxₙ en rendent une part.")

    # (d) ablation
    ax = panneau("termes_gain.png", gauche=0.31, droite=0.95)
    style(ax)
    termes = list(m3.termes)
    r = ablation(Z, N, B, termes)
    y = np.arange(len(termes))[::-1]
    ax.plot(r, y, color=SERIES[0], lw=2)
    ax.plot(r, y, "o", ms=8, color=SERIES[0], mec=SURFACE, mew=2)
    # valeur en bas a droite du point : le segment qui monte vers le terme
    # precedent part a droite, celui qui descend vers le suivant part a gauche
    for xi, yi in zip(r, y):
        ax.text(xi * 1.07, yi - 0.34, f"{xi:.2f}".replace(".", ","), va="center",
                fontsize=10.5, color=ENCRE)
    ax.set_yticks(y)
    ax.set_yticklabels(["+ " + NOMS[k] if i else NOMS[k] + " seul" for i, k in enumerate(termes)])
    ax.set_xscale("log")
    ax.set_xlim(0.8, max(r) * 2.2)
    ticks = [1, 2, 5, 10, 20, 50]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel("écart RMS sur B (MeV, échelle logarithmique)")
    ax.grid(True, axis="x", color=GRILLE, lw=0.8)
    ax.grid(False, axis="y")
    titre(ax, "(d) Ce que rapporte chaque terme",
          "Écart RMS sur les 2367 masses mesurées, en ajoutant les termes un à un.")

    for fichier, fig in figures.items():
        if fichier == "termes_gain.png":
            note = ("Ajustements successifs sur les 2367 masses mesurées d'AME2020 (Z, N ≥ 8), "
                    "un terme de plus à chaque ligne.")
        else:
            note = ("Coefficients de M3 ajustés sur AME2020, masses mesurées seulement (Z, N ≥ 8) ; "
                    "vallée = isobare le plus lié pour chaque A.")
        fig.text(0.02, 0.02, note, fontsize=8.5, color=DISCRET)
        fig.savefig(fichier, dpi=120, facecolor=SURFACE)
        plt.close(fig)
    return r


# ==========================================================================
def figure_carte(Z, N, contrib):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    seq = LinearSegmentedColormap.from_list("bleus", BLEUS[::-1])          # fonce = negatif
    seq_pos = LinearSegmentedColormap.from_list("bleus_pos", BLEUS)         # fonce = positif
    div = LinearSegmentedColormap.from_list("div", [SERIES[0], GRIS_MILIEU, ROUGE])

    panneaux = [
        ("Coulomb", "Coulomb : croît avec Z", seq, None),
        ("asymetrie", "Asymétrie : loin de N = Z", seq, None),
        ("sym. surface", "Symétrie de surface : rend une part", seq_pos, None),
        ("Wigner", "Wigner : nul en N = Z, une pointe de liaison", seq, None),
        ("appariement", "Appariement : un damier pair/impair", div, "sym"),
        ("couches", "Couches : s'effacent aux nombres magiques", seq, None),
    ]
    c = dict(contrib)
    c["couches"] = contrib["couches S"] + contrib["couches S^2"] + contrib["couches p-n"]

    fig, axes = plt.subplots(3, 2, figsize=(11, 15), facecolor=SURFACE)
    plt.subplots_adjust(left=0.08, right=0.91, top=0.86, bottom=0.05, wspace=0.42, hspace=0.45)
    for ax, (cle, texte, cmap, mode) in zip(axes.flat, panneaux):
        style(ax)
        ax.grid(False)
        v = c[cle]
        if mode == "sym":
            m = float(np.max(np.abs(v)))
            norm = TwoSlopeNorm(vcenter=0.0, vmin=-m, vmax=m)
        else:
            norm = None
        sc = ax.scatter(N, Z, c=v, cmap=cmap, norm=norm, s=9, marker="s", linewidths=0)
        for mg in (8, 20, 28, 50, 82, 126):
            ax.axvline(mg, color=GRILLE, lw=0.8, zorder=0)
            if mg <= 100:
                ax.axhline(mg, color=GRILLE, lw=0.8, zorder=0)
        ax.plot([0, 100], [0, 100], color=AXE, lw=0.8, zorder=0)
        ax.text(92, 98, "N = Z", fontsize=10, color=DISCRET, ha="right")
        ax.set_xlim(0, 162)
        ax.set_ylim(0, 105)
        ax.set_xlabel("N")
        ax.set_ylabel("Z")
        cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=10.5, colors=DISCRET, labelcolor=ENCRE_2)
        cb.set_label("contribution à B (MeV)", color=ENCRE_2, fontsize=11)
        titre(ax, texte, f"de {nombre(v.min())} à {nombre(v.max())} MeV")

    fig.text(0.08, 0.965, "Où agit chaque terme : les contributions de M3",
             fontsize=16, color=ENCRE, fontweight="semibold")
    fig.text(0.08, 0.945, "Un carré par noyau mesuré. Lignes fines : nombres magiques.\n"
             "Plus foncé = contribution plus grande en valeur absolue.",
             fontsize=11.5, color=ENCRE_2, va="top")
    fig.text(0.08, 0.012, "Données : AME2020, masses mesurées seulement (Z, N ≥ 8) -- "
             "M. Wang et al., Chin. Phys. C 45 (2021) 030003", fontsize=8.5, color=DISCRET)
    fig.savefig("termes_carte.png", dpi=110, facecolor=SURFACE)
    plt.close(fig)


# ==========================================================================
def tableau(Z, N, B, m3, contrib):
    """La table de valeurs (equivalent texte des figures)."""
    print("\n  Contribution de chaque terme de M3 a B, en MeV :\n")
    noyaux = (("O-16", 8, 8), ("Ca-48", 20, 28), ("Sn-120", 50, 70),
              ("Pb-208", 82, 126), ("U-238", 92, 146))
    print(f"  {'terme':<22}" + "".join(f"{n:>10}" for n, _, _ in noyaux))
    idx = []
    for _, z, n in noyaux:
        idx.append(int(np.where((Z == z) & (N == n))[0][0]))
    for k in m3.termes:
        print(f"  {NOMS[k]:<22}" + "".join(f"{contrib[k][i]:10.1f}" for i in idx))
    print(f"  {'total M3':<22}" + "".join(f"{sum(contrib[k][i] for k in m3.termes):10.1f}" for i in idx))
    print(f"  {'mesure AME2020':<22}" + "".join(f"{B[i]:10.1f}" for i in idx))


if __name__ == "__main__":
    Z, N, B = charger()
    m3 = ModeleLineaire("M3", colonnes_M3).ajuster(Z, N, B)
    contrib = contributions(Z, N, m3)
    r = figure_vallee(Z, N, B, m3, contrib)
    print("  Ecart RMS en ajoutant les termes un a un (MeV) :")
    for k, x in zip(m3.termes, r):
        print(f"    + {NOMS[k]:<22} {x:8.3f}")
    figure_carte(Z, N, contrib)
    tableau(Z, N, B, m3, contrib)
    print("\n  --> termes_goutte.png, termes_m2.png, termes_couches.png, "
          "termes_gain.png, termes_carte.png ecrits")
