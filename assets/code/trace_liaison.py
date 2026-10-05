"""
trace_liaison.py -- B/A : donnees mesurees brutes face au modele de la goutte
liquide, et construction du modele terme par terme.

Lancer :  python trace_liaison.py   ->  energie_liaison_modeles.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from donnees_liaison import table_liaison, chercher
from fission import (liaison_par_nucleon, A_VOLUME, A_SURFACE, A_COULOMB,
                     A_ASYMETRIE)

MAGIQUES = (20, 28, 50, 82, 126)


def z_optimal(A):
    """Z qui maximise B a A fixe, selon la goutte liquide (sans appariement).

    dB/dZ = 0  ->  Z0 = A / (2 + (a_C / (2 a_A)) * A^(2/3))
    C'est la vallee de stabilite PREDITE par le modele.
    """
    return A / (2.0 + (A_COULOMB / (2.0 * A_ASYMETRIE)) * A ** (2.0 / 3.0))


def termes_par_nucleon(A, Z):
    """Contributions de chaque terme a B/A (sans appariement), en MeV."""
    vol = np.full_like(A, A_VOLUME, dtype=float)
    surf = -A_SURFACE * A ** (-1.0 / 3.0)
    coul = -A_COULOMB * Z ** 2 / A ** (4.0 / 3.0)
    asym = -A_ASYMETRIE * (A - 2 * Z) ** 2 / A ** 2
    return vol, surf, coul, asym


def main():
    t = table_liaison()
    nat = t[t["abondance_pct"] > 0]

    A_lisse = np.linspace(4, 260, 600)
    Z_lisse = z_optimal(A_lisse)
    vol, surf, coul, asym = termes_par_nucleon(A_lisse, Z_lisse)
    semf_lisse = vol + surf + coul + asym

    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.32, wspace=0.22)

    # ------------------------------------------------------------------
    # (a) donnees brutes + modele
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    mes, est = t[~t["estimee"]], t[t["estimee"]]
    ax.plot(mes["A"], mes["B_sur_A_MeV"], ".", ms=2, color="#b4c2d3",
            label=f"AME2020, masses mesurees ({len(mes)})", zorder=1)
    ax.plot(est["A"], est["B_sur_A_MeV"], "x", ms=2.5, mew=0.6,
            color="#e9a96b",
            label=f"AME2020, masses estimees '#' ({len(est)})", zorder=1)
    ax.plot(nat["A"], nat["B_sur_A_MeV"], "o", ms=3.2, color="#1f4e79",
            label=f"isotopes naturels ({len(nat)}, tous mesures)", zorder=3)
    ax.plot(A_lisse, semf_lisse, "-", color="crimson", lw=2.2,
            label="goutte liquide (Bethe-Weizsacker), vallee predite", zorder=4)

    notes = [("He", 4, (12, -4), "He-4 : doublement magique\nmodele trop bas de 1,6 MeV"),
             ("Ni", 62, (25, -60), "Ni-62 : maximum mesure\n8,7946 MeV"),
             ("U", 235, (-40, -55), "U-235 : 7,5909 MeV")]
    for s, a, decal, texte in notes:
        r = chercher(t, s, a)
        ax.plot(r["A"], r["B_sur_A_MeV"], "o", ms=9, mfc="none",
                mec="black", mew=1.6, zorder=5)
        ax.annotate(texte, (r["A"], r["B_sur_A_MeV"]),
                    textcoords="offset points", xytext=decal, fontsize=9,
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.8))

    # fleche fission / fusion
    pd = chercher(t, "Pd", 118)
    ax.plot(pd["A"], pd["B_sur_A_MeV"], "o", ms=9, mfc="none",
            mec="darkorange", mew=1.8, zorder=5)
    ax.annotate("", xy=(118, pd["B_sur_A_MeV"] - 0.08), xytext=(236, 7.50),
                arrowprops=dict(arrowstyle="->", color="darkorange", lw=2.2,
                                connectionstyle="arc3,rad=-0.35"))
    ax.text(178, 6.05, "FISSION\nU-236 -> 2 x Pd-118\n+0,8 MeV/nucleon x 236\n= 193 MeV",
            color="darkorange", fontsize=9.5, ha="center", fontweight="bold")
    ax.annotate("", xy=(4, 7.07), xytext=(2, 1.11),
                arrowprops=dict(arrowstyle="->", color="seagreen", lw=2))
    ax.text(14, 3.0, "FUSION\nD + D -> He-4\n+6 MeV/nucleon", color="seagreen",
            fontsize=9.5, fontweight="bold")

    ax.set_xlim(0, 262)
    ax.set_ylim(0, 9.4)
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("B/A  [MeV par nucleon]")
    ax.set_title("Energie de liaison par nucleon : mesures brutes et modele "
                 "de la goutte liquide", fontsize=12.5)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)

    # ------------------------------------------------------------------
    # (b) construction du modele terme par terme
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    c1 = vol
    c2 = c1 + surf
    c3 = c2 + coul
    c4 = c3 + asym
    ax.plot(A_lisse, c1, color="#555555", lw=1.6, ls="--",
            label=f"volume seul  (+{A_VOLUME})")
    ax.plot(A_lisse, c2, color="#2a9d8f", lw=1.8, label="- surface")
    ax.plot(A_lisse, c3, color="#e76f51", lw=1.8, label="- Coulomb")
    ax.plot(A_lisse, c4, color="crimson", lw=2.4, label="- asymetrie = modele")
    ax.fill_between(A_lisse, c1, c2, color="#2a9d8f", alpha=0.15)
    ax.fill_between(A_lisse, c2, c3, color="#e76f51", alpha=0.18)
    ax.fill_between(A_lisse, c3, c4, color="crimson", alpha=0.12)
    ax.plot(nat["A"], nat["B_sur_A_MeV"], ".", ms=2.5, color="#1f4e79",
            label="mesures (naturels)")
    ax.text(150, 13.3, "surface : penalise\nles PETITS noyaux",
            color="#2a9d8f", fontsize=8.5)
    ax.text(150, 10.2, "Coulomb : penalise\nles GROS noyaux", color="#e76f51",
            fontsize=8.5)
    ax.set_xlim(0, 262)
    ax.set_ylim(4, 16.5)
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("contribution cumulee a B/A  [MeV]")
    ax.set_title("Le modele, terme par terme\n(le maximum nait du duel "
                 "surface / Coulomb)")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)

    # ------------------------------------------------------------------
    # (c) residu mesure - modele
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    m = (t["A"] >= 16) & (~t["estimee"])     # residu : masses mesurees seulement
    sub = t[m]
    res = sub["B_sur_A_MeV"] - liaison_par_nucleon(sub["A"].astype(float),
                                                  sub["Z"].astype(float))
    sc = ax.scatter(sub["N"], res * 1000, c=sub["Z"], s=3, cmap="viridis")
    for N in MAGIQUES:
        ax.axvline(N, color="crimson", ls=":", lw=1.2)
        ax.text(N + (-1.5 if N == 20 else 1.5 if N == 28 else 0), 330, str(N),
                color="crimson", fontsize=9,
                ha="right" if N == 20 else "left" if N == 28 else "center")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylim(-150, 360)
    ax.set_xlabel("nombre de neutrons N")
    ax.set_ylabel("mesure - modele  [keV par nucleon]")
    ax.set_title("Ce que le modele ignore : les couches\n"
                 f"(masses mesurees seules, {len(sub)} nuclides)")
    fig.colorbar(sc, ax=ax, label="Z")
    ax.grid(alpha=0.3)

    fig.text(0.5, 0.005, "Donnees : AME2020 -- M. Wang et al., Chinese Physics C 45 (2021) "
             "030003, via le paquet Python periodictable 2.1.0 (massround.mas20.txt)",
             ha="center", fontsize=8.5, color="#444444")
    fig.savefig("energie_liaison_modeles.png", dpi=130, bbox_inches="tight")
    print("energie_liaison_modeles.png ecrit")

    # quelques chiffres pour verifier les annotations
    for s, a in (("He", 4), ("Ni", 62), ("U", 235), ("U", 236), ("Pd", 118)):
        r = chercher(t, s, a)
        print(f"{s}-{a}: mesure {r['B_sur_A_MeV']:.4f}  modele "
              f"{float(liaison_par_nucleon(a, r['Z'])):.4f}")
    r = chercher(t, "H", 2)
    print(f"H-2: {r['B_sur_A_MeV']:.4f}")
    i = int(np.argmax(semf_lisse))
    print(f"maximum du modele lisse : A = {A_lisse[i]:.0f}, B/A = {semf_lisse[i]:.3f}")


if __name__ == "__main__":
    main()
