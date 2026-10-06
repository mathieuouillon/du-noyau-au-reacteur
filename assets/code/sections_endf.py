"""
sections_endf.py -- Les vraies sections efficaces de l'U-235 et de l'U-238,
reconstruites a partir de l'evaluation ENDF/B-VIII.0.

  1. Telecharge, depuis le NNDC (Brookhaven), les deux seuls fichiers
     n-092_U_235.endf et n-092_U_238.endf de l'archive ENDF/B-VIII.0
     (lecture partielle de l'archive zip : ~16 Mo transferes au lieu de 296).
  2. Reconstruit les resonances resolues (Reich-Moore, reconstruction.py)
     et les elargit par effet Doppler a 293,6 K.
  3. Valide : valeurs thermiques annoncees par l'evaluation, facteur de
     Westcott, integrales de resonance (comparees a JENDL-4.0).
  4. Trace : vue d'ensemble, foret de resonances, effet Doppler.

Dependances : numpy, matplotlib, endf (lecture ENDF-6), remotezip.
Lancer :  python sections_endf.py   (~1 min, plus le telechargement initial)
    ->  endf_vue_ensemble.png, endf_resonances.png, endf_doppler.png

Donnees : D.A. Brown et al., ENDF/B-VIII.0, Nucl. Data Sheets 148 (2018) 1.
"""

import os
import warnings
import numpy as np
from reconstruction import Noyau, KB
from sections_efficaces import breit_wigner_capture, U238_RES
from trace_termes import (SURFACE, ENCRE, ENCRE_2, DISCRET, GRILLE, AXE, SERIES,
                          style, titre, etiquettes)

warnings.filterwarnings("ignore", message="MF=")
URL = "https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_neutrons.zip"
DOSSIER = "donnees_endf"
FICHIERS = {"U-235": "n-092_U_235.endf", "U-238": "n-092_U_238.endf"}
T_AMB = 293.6            # K : temperature de reference des sections "thermiques"
E0 = KB * T_AMB          # 0,0253 eV


def telecharger():
    """Extrait les deux fichiers de l'archive distante, sans la telecharger
    en entier (requetes HTTP partielles)."""
    os.makedirs(DOSSIER, exist_ok=True)
    manquants = [f for f in FICHIERS.values() if not os.path.exists(os.path.join(DOSSIER, f))]
    if not manquants:
        return
    from remotezip import RemoteZip
    with RemoteZip(URL) as z:
        for info in z.infolist():
            nom = os.path.basename(info.filename)
            if nom in manquants:
                print(f"  telechargement de {nom} ({info.compress_size / 1e6:.1f} Mo compresses)")
                with z.open(info) as src, open(os.path.join(DOSSIER, nom), "wb") as dst:
                    dst.write(src.read())


def titre_section(t):
    print()
    print("=" * 76)
    print(t)
    print("=" * 76)


def valeur(E, s, e):
    return float(np.interp(e, E, s))


def westcott(E, s, T=T_AMB):
    """g = <sigma v> / (sigma(E0) v0) pour une maxwellienne a T (E0 = kT)."""
    kT = KB * T
    m = E < 40 * kT
    e, sig = E[m], s[m]
    M = 2 / np.sqrt(np.pi) * kT ** -1.5 * np.sqrt(e) * np.exp(-e / kT)
    return float(np.trapezoid(sig * np.sqrt(e / E0) * M, e) / valeur(E, s, E0))


def integrale_resonance(E, s, emin=0.5, emax=2e7):
    m = (E >= emin) & (E <= emax)
    return float(np.trapezoid(s[m] / E[m], E[m]))


# ==========================================================================
def calculer():
    donnees = {}
    for nom, fichier in FICHIERS.items():
        n = Noyau(os.path.join(DOSSIER, fichier))
        Ef, s0 = n.zero_kelvin(fine=True)          # integrales, Doppler
        E, _ = n.zero_kelvin(fine=False)           # affichage
        sT = n.temperature(Ef, s0, T_AMB, E_sortie=E)
        s0_aff = {k: np.interp(E, Ef, v) for k, v in s0.items()}
        donnees[nom] = dict(noyau=n, Ef=Ef, s0f=s0, E=E, s0=s0_aff, sT=sT)
        print(f"  {nom} : {len(Ef)} points a 0 K, {len(E)} points a {T_AMB} K, "
              f"resonances resolues jusqu'a {n.rm.EH:g} eV")
    return donnees


def validations(d):
    titre_section("1. LA RECONSTRUCTION EST-ELLE JUSTE ?")
    u5, u8 = d["U-235"], d["U-238"]
    print("""
  Valeurs a 0,0253 eV et 0 K annoncees dans l'en-tete du fichier U-235
  (MF=1, MT=451), contre la reconstruction :
""")
    print(f"    fission U-235 : evaluation 586.7870 b, reconstruit {valeur(u5['Ef'], u5['s0f']['fission'], E0):.4f} b")
    print(f"    capture U-235 : evaluation  99.3909 b, reconstruit {valeur(u5['Ef'], u5['s0f']['capture'], E0):.4f} b")
    print(f"\n  A {T_AMB} K (effet Doppler), contre JENDL-4.0 a 300 K (autre evaluation) :")
    for nom, reac, j in (("U-235", "fission", 585.1), ("U-235", "capture", 98.71),
                         ("U-238", "capture", 2.683), ("U-238", "elastique", 9.300)):
        x = d[nom]
        print(f"    {reac:<9} {nom} : {valeur(x['E'], x['sT'][reac], E0):9.3f} b   (JENDL-4.0 : {j} b)")
    print("\n  Facteur de Westcott g (maxwellienne a 293,6 K) :")
    for nom, reac, j in (("U-235", "fission", 0.977), ("U-235", "capture", 0.991),
                         ("U-238", "capture", 1.003)):
        x = d[nom]
        print(f"    {reac:<9} {nom} : g = {westcott(x['E'], x['sT'][reac]):.4f}   (JENDL-4.0 : {j})")
    print("\n  Integrales de resonance, de 0,5 eV a 20 MeV, a 0 K (JENDL-4.0 : autre")
    print("  evaluation, bornes d'integration non precisees dans la table) :")
    for nom, reac, j in (("U-235", "fission", 274.4), ("U-235", "capture", 139.0),
                         ("U-238", "capture", 275.6)):
        x = d[nom]
        print(f"    {reac:<9} {nom} : {integrale_resonance(x['Ef'], x['s0f'][reac]):8.2f} b   (JENDL-4.0 : {j} b)")


def lecture(d):
    titre_section("2. CE QUE DISENT LES COURBES")
    u5, u8 = d["U-235"], d["U-238"]
    f5 = lambda e: valeur(u5["E"], u5["sT"]["fission"], e)
    f8 = lambda e: valeur(u8["E"], u8["sT"]["fission"], e)
    print(f"""
  Fission de l'U-235 : {f5(E0):.1f} b a 0,0253 eV, {f5(1e6):.3f} b a 1 MeV,
  {f5(2e6):.3f} b a 2 MeV : rapport thermique / 2 MeV = {f5(E0) / f5(2e6):.0f}.
  Fission de l'U-238 :""")
    for e in (0.5e6, 1e6, 1.2e6, 1.5e6, 2e6, 5e6):
        print(f"    {e / 1e6:4.1f} MeV : {f8(e):.4g} b")
    Es = u8["E"][(u8["E"] > 1e5)]
    sf = np.interp(Es, u8["E"], u8["sT"]["fission"])
    seuil = Es[np.argmax(sf > 0.5 * sf[Es > 3e6].mean())]
    print(f"  la fission de l'U-238 atteint la moitie de son plateau (2-5 MeV) vers {seuil / 1e6:.2f} MeV")
    # foret de resonances
    for nom, x, reac in (("U-235", u5, "fission"), ("U-238", u8, "capture")):
        sec = np.concatenate([np.asarray(s["ER"]) for s in x["noyau"].rm.sections])
        print(f"  {nom} : {np.sum((sec > 0) & (sec < x['noyau'].rm.EH))} resonances resolues "
              f"(au-dessus de 0) jusqu'a {x['noyau'].rm.EH:g} eV")
    return seuil


def effet_doppler(d):
    titre_section("3. L'EFFET DOPPLER SUR LA RESONANCE A 6,67 eV")
    u8 = d["U-238"]
    n = u8["noyau"]
    Ef, s0 = u8["Ef"], u8["s0f"]
    E = np.linspace(6.0, 7.4, 4001)                 # grille d'affichage du zoom
    courbes = {0: np.interp(E, Ef, s0["capture"])}
    for T in (T_AMB, 900.0, 1800.0):
        courbes[T] = n.temperature(Ef, s0, T, reactions=("capture",), E_sortie=E,
                                   emin=5.5, emax=8.0)["capture"]
    u8["E_zoom"] = E
    m = (E >= 6.0) & (E <= 7.4)
    print(f"\n  {'T (K)':>8} {'pic (b)':>10} {'largeur a mi-hauteur (eV)':>27} {'aire 6,0-7,4 eV (b eV)':>24}")
    res = {}
    for T, s in courbes.items():
        e, c = E[m], s[m]
        pic = c.max()
        au_dessus = e[c >= pic / 2]
        fwhm = au_dessus.max() - au_dessus.min()
        aire = np.trapezoid(c, e)
        res[T] = (pic, fwhm, aire)
        print(f"  {T:>8.1f} {pic:>10.0f} {fwhm:>27.4f} {aire:>24.1f}")
    print("""
  Le pic s'abaisse et s'elargit quand la temperature monte, mais l'aire
  reste a peu pres la meme. Dans un combustible, le flux est creuse au
  centre de la resonance (autoprotection) : une resonance plus large et
  moins haute absorbe donc PLUS de neutrons. C'est le coefficient Doppler,
  negatif et instantane, premiere barriere de surete d'un REP.""")
    return courbes, res


# ==========================================================================
def figures(d, seuil, courbes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

    def fmt_e(v, _):
        if v < 1e-3:
            return f"10$^{{{int(round(np.log10(v)))}}}$ eV"
        if v >= 1e6:
            return f"{v / 1e6:g} MeV"
        if v >= 1e3:
            return f"{v / 1e3:g} keV"
        return f"{v:g} eV".replace(".", ",")

    def puissance(v, _):
        e = int(round(np.log10(v)))
        if abs(np.log10(v) - e) > 1e-6:
            return ""
        if -1 <= e <= 3:
            return f"{v:g}".replace(".", ",")
        return f"10$^{{{e}}}$"

    def figure(gauche=0.11, droite=0.78):
        fig, ax = plt.subplots(figsize=(9, 5.6), facecolor=SURFACE)
        fig.subplots_adjust(left=gauche, right=droite, top=0.83, bottom=0.15)
        return fig, ax

    def note(fig, texte):
        fig.text(0.02, 0.02, texte, fontsize=8.5, color=DISCRET)

    u5, u8 = d["U-235"], d["U-238"]

    # --- 1. vue d'ensemble
    fig, ax = figure()
    style(ax)
    ax.axvspan(1e-3, 0.5, color=GRILLE, alpha=0.45, lw=0)
    ax.axvspan(1e5, 1e7, color=GRILLE, alpha=0.45, lw=0)
    ax.text(0.022, 2e-5, "thermique", ha="center", fontsize=10, color=ENCRE_2)
    ax.text(1e6, 2e-5, "rapide", ha="center", fontsize=10, color=ENCRE_2)
    ax.text(40, 2e-5, "résonances", ha="center", fontsize=10, color=ENCRE_2)
    series = [(u5["E"], u5["sT"]["fission"], "fission U-235"),
              (u8["E"], u8["sT"]["capture"], "capture U-238"),
              (u8["E"], u8["sT"]["fission"], "fission U-238")]
    for i, (E, s, _) in enumerate(series):
        m = s > 0
        ax.plot(E[m], s[m], color=SERIES[i], lw=0.9 if i < 2 else 1.4)
    e = np.logspace(-5, 0, 50)
    ax.plot(e, 3 * valeur(u5["E"], u5["sT"]["fission"], E0) * np.sqrt(E0 / e), color=DISCRET,
            lw=0.8, zorder=1)
    ax.text(3e-4, 3 * valeur(u5["E"], u5["sT"]["fission"], E0) * np.sqrt(E0 / 3e-4) * 1.6,
            "pente 1/v", fontsize=9.5, color=ENCRE_2, rotation=-14)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-5, 2e7)
    ax.set_ylim(1e-5, 1e5)
    ax.xaxis.set_major_locator(FixedLocator([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5, 1e7]))
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_e))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(FixedLocator([10.0 ** k for k in range(-5, 6)]))
    ax.yaxis.set_major_formatter(FuncFormatter(puissance))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("énergie du neutron")
    ax.set_ylabel("section efficace (barns)")
    etiquettes(ax, None, [valeur(E, s, 1.9e7) for E, s, _ in series], [n for *_, n in series],
               SERIES[:3], ecart=0.45)
    titre(ax, "Les vraies sections efficaces de l'uranium (ENDF/B-VIII.0)",
          "Fission de l'U-235 : 586 b en thermique, ~1 b en rapide. L'U-238 ne fissionne qu'au-dessus d'environ 1 MeV.")
    note(fig, "Résonances résolues : reconstruction Reich-Moore et Doppler à 293,6 K (reconstruction.py). "
              "Au-delà : sections tabulées.")
    fig.savefig("endf_vue_ensemble.png", dpi=120, facecolor=SURFACE)
    plt.close(fig)

    # --- 2. foret de resonances + modele a une resonance
    fig, ax = figure(droite=0.95)
    style(ax)
    for i, (E, s) in enumerate(((u5["E"], u5["sT"]["fission"]), (u8["E"], u8["sT"]["capture"]))):
        m = (E > 0.01) & (E < 100)
        ax.plot(E[m], s[m], color=SERIES[i], lw=1.1)
    e = np.logspace(-2, 2, 4000)
    ax.plot(e, breit_wigner_capture(e, **U238_RES), color=ENCRE_2, lw=1)
    ax.text(0.012, 1.5e3, "fission de l'U-235", fontsize=10.5, color=ENCRE)
    ax.text(0.012, 4.0, "capture de l'U-238", fontsize=10.5, color=ENCRE)
    ax.annotate("modèle à une seule résonance\n(Breit-Wigner, leçon 2)", (0.3, breit_wigner_capture(0.3, **U238_RES)),
                xytext=(0.05, 0.03), fontsize=9.5, color=ENCRE_2,
                arrowprops=dict(arrowstyle="-", color=DISCRET, lw=0.6))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0.01, 1e5)
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_e))
    ax.yaxis.set_major_locator(FixedLocator([10.0 ** k for k in range(-2, 6)]))
    ax.yaxis.set_major_formatter(FuncFormatter(puissance))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("énergie du neutron")
    ax.set_ylabel("section efficace à 293,6 K (barns)")
    titre(ax, "Une forêt de résonances",
          "Au-delà de 1 eV, chaque pic est un état du noyau composé ; le modèle à une résonance n'en décrit qu'un.")
    note(fig, "ENDF/B-VIII.0 reconstruit (reconstruction.py). Modèle : résonance de 6,67 eV seule, "
              "paramètres de Mughabghab, sans effet Doppler.")
    fig.savefig("endf_resonances.png", dpi=120, facecolor=SURFACE)
    plt.close(fig)

    # --- 3. Doppler
    fig, ax = figure(droite=0.78)
    style(ax)
    E = u8["E_zoom"]
    m = (E > 6.2) & (E < 7.15)
    noms = {0: "0 K", T_AMB: "293,6 K", 900.0: "900 K", 1800.0: "1 800 K"}
    ordre = [0, T_AMB, 900.0, 1800.0]
    for i, T in enumerate(ordre):
        ax.plot(E[m], courbes[T][m] / 1e3, color=SERIES[i], lw=2)
    ax.set_xlim(6.2, 7.15)
    ax.set_ylim(0, 24)
    ax.set_xlabel("énergie du neutron (eV)")
    ax.set_ylabel("capture de l'U-238 (milliers de barns)")
    pics = [float(courbes[T][m].max()) / 1e3 for T in ordre]
    etiquettes(ax, None, pics, [f"{noms[T]} : {p:.1f} kb".replace(".", ",") for T, p in zip(ordre, pics)],
               SERIES[:4], ecart=1.6)
    titre(ax, "L'effet Doppler sur la résonance de 6,67 eV",
          "Quand la température monte, le pic s'abaisse et s'élargit ; son aire ne change presque pas.")
    note(fig, "ENDF/B-VIII.0, élargissement en gaz libre (noyau exact). 900 et 1 800 K : "
              "ordres de grandeur d'un combustible en fonctionnement.")
    fig.savefig("endf_doppler.png", dpi=120, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    telecharger()
    d = calculer()
    validations(d)
    seuil = lecture(d)
    courbes, _ = effet_doppler(d)
    figures(d, seuil, courbes)
    print("\n  --> endf_vue_ensemble.png, endf_resonances.png, endf_doppler.png")
