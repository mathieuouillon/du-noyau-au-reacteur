"""
galerie_orbitales.py -- Tracer les orbitales des nucleons comme on trace
celles de l'electron, puis la forme du noyau entier.

Produit : orbitales_nucleaires.png  et  orbitales_3d.json (pour la page 3D).
Lancer  : python galerie_orbitales.py
"""

import json
import numpy as np
from orbitales import (etats, densite_orbitale, densite_noyau,
                       longueur_oscillateur, rayon_rms)
from nilsson import hbar_omega

R_MAX = 11.0
rho1 = np.linspace(0, R_MAX, 221)
z1 = np.linspace(-R_MAX, R_MAX, 441)
RHO, Z = np.meshgrid(rho1, z1)


def plein_plan(d):
    """Densite de revolution -> coupe complete (x, z) par symetrie."""
    return np.hstack([d[:, :0:-1], d])


def choisir(E, parent, om_x2):
    cands = [e for e in E if e["parent"] == parent and e["om_x2"] == om_x2]
    return max(cands, key=lambda e: e["poids_parent"]) if cands else None


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    export = {"grille": {"rho": rho1[::4].round(3).tolist(),
                         "z": z1[::4].round(3).tolist()}, "objets": []}

    def exporter(nom, groupe, texte, d):
        export["objets"].append({"nom": nom, "groupe": groupe, "texte": texte,
                                 "densite": (d[::4, ::4] / d.max()).round(4).tolist()})

    fig = plt.figure(figsize=(16, 11.2))
    gs = fig.add_gridspec(3, 6, hspace=0.22, wspace=0.08, height_ratios=[1, 1, 1.55])
    ext = [-R_MAX, R_MAX, -R_MAX, R_MAX]

    def tracer(ax, d, titre, cmap="magma"):
        ax.imshow(plein_plan(d), origin="lower", extent=ext, cmap=cmap, aspect="equal")
        ax.set_title(titre, fontsize=9.5)
        ax.set_xticks([])
        ax.set_yticks([])

    # ---------------------------------------------------------------- spheriques
    Zn, Nn = 82, 126
    b0 = longueur_oscillateur(Zn, Nn, "neutron", 0.0)
    E0 = etats(0.0, "neutron")
    choix = [("1s1/2", 1), ("1p3/2", 1), ("1p3/2", 3),
             ("1d5/2", 1), ("1d5/2", 5), ("2s1/2", 1)]
    print("Ligne 1 : orbitales spheriques (noyau de masse 208)")
    for k, (p, o) in enumerate(choix):
        e = choisir(E0, p, o)
        d = densite_orbitale(e, RHO, Z, b0)
        titre = f"{p}, $\\Omega$ = {o}/2"
        tracer(fig.add_subplot(gs[0, k]), d, titre)
        exporter(f"{p} Ω={o}/2", "Orbitales sphériques", titre, d)
        print(f"  {p:7s} Omega={o}/2  rms = {rayon_rms(d, RHO, Z)[0]:.2f} fm")

    # ---------------------------------------------------------------- 1i13/2 deforme
    delta = 0.25
    Zu, Nu = 92, 146
    bd = longueur_oscillateur(Zu, Nu, "neutron", delta)
    hw = hbar_omega(Zu, Nu, "neutron")
    Ed = etats(delta, "neutron")
    famille = [choisir(Ed, "1i13/2", o) for o in (1, 5, 9, 13)]
    e_ref = min(e["energie"] for e in famille)
    print(f"\nLigne 2 : famille 1i13/2 a delta = {delta} (neutrons, noyau de masse 238)")
    for k, (e, o) in enumerate(zip(famille, (1, 5, 9, 13))):
        d = densite_orbitale(e, RHO, Z, bd)
        dE = (e["energie"] - e_ref) * hw
        titre = f"1i13/2, $\\Omega$ = {o}/2\n+{dE:.1f} MeV"
        ax = fig.add_subplot(gs[1, k])
        tracer(ax, d, titre, cmap="viridis")
        exporter(f"1i13/2 Ω={o}/2 (δ=0,25)", "Noyau allongé : famille 1i13/2",
                 f"1i13/2, Ω = {o}/2, +{dE:.1f} MeV", d)
        print(f"  Omega={o:2d}/2  energie relative = +{dE:5.2f} MeV  "
              f"(poids du parent spherique {e['poids_parent']:.0%})")
    ax = fig.add_subplot(gs[1, 4:])
    ax.axis("off")
    ax.text(0.02, 0.5,
            "Les mêmes quatorze états sphériques 1i13/2,\n"
            "dans un noyau allongé selon l'axe vertical.\n\n"
            "Ω petit : le moment angulaire est presque\n"
            "perpendiculaire à l'axe ; l'orbite passe par\n"
            "les pôles, le long du noyau, là où il y a\n"
            "de la matière. Énergie la plus basse.\n\n"
            "Ω grand : le moment angulaire est aligné\n"
            "sur l'axe ; l'orbite tourne dans le plan\n"
            "équatorial, là où le noyau est mince.\n"
            "Énergie la plus haute.\n\n"
            "C'est l'éventail du diagramme de Nilsson.",
            fontsize=10, va="center")

    # ---------------------------------------------------------------- noyaux entiers
    print("\nLigne 3 : densite totale de nucleons")
    for k, (nom, z, n, dl) in enumerate((("Pb-208", 82, 126, 0.0),
                                         ("Er-166", 68, 98, 0.241),
                                         ("U-238", 92, 146, 0.193))):
        D = densite_noyau(z, n, dl, RHO, Z)
        tot = D["proton"] + D["neutron"]
        ax = fig.add_subplot(gs[2, 2 * k:2 * k + 2])
        im = ax.imshow(plein_plan(tot), origin="lower", extent=ext, cmap="inferno",
                       aspect="equal", vmin=0, vmax=0.22)
        ax.contour(np.linspace(-R_MAX, R_MAX, 2 * len(rho1) - 1), z1, plein_plan(tot),
                   levels=[0.08], colors="cyan", linewidths=1)
        w = 2 * np.pi * RHO * tot
        rz = np.sqrt(np.sum(w * Z ** 2) / np.sum(w))
        rr = np.sqrt(np.sum(w * RHO ** 2) / (2 * np.sum(w)))
        ax.set_title(f"{nom}, $\\delta$ = {dl:.2f}\n"
                     f"rapport des rayons z / x = {rz / rr:.2f}", fontsize=10)
        ax.set_xlabel("x  [fm]")
        if k == 0:
            ax.set_ylabel("z (axe de symétrie)  [fm]")
        fig.colorbar(im, ax=ax, fraction=0.046, label="nucléons / fm³")
        exporter(nom, "Noyaux entiers", f"{nom}, δ = {dl:.2f}", tot)
        rp = rayon_rms(D["proton"], RHO, Z)[0]
        rn = rayon_rms(D["neutron"], RHO, Z)[0]
        print(f"  {nom:7s} rms protons {rp:.2f} fm, neutrons {rn:.2f} fm, "
              f"z/x = {rz / rr:.3f}, densite max {tot.max():.3f} fm^-3")

    fig.suptitle("Orbitales nucléaires (modèle de Nilsson) : densité de probabilité "
                 "dans un plan contenant l'axe de symétrie", fontsize=13)
    fig.text(0.5, 0.075, "Contour cyan : 0,08 nucléon/fm³, la moitié de la densité de "
             "la matière nucléaire. Chaque image fait 22 fm de côté.",
             ha="center", fontsize=9.5, color="#333333")
    fig.savefig("orbitales_nucleaires.png", dpi=115, bbox_inches="tight")
    with open("orbitales_3d.json", "w") as f:
        json.dump(export, f, separators=(",", ":"))
    print("\n--> orbitales_nucleaires.png, orbitales_3d.json")


if __name__ == "__main__":
    main()
