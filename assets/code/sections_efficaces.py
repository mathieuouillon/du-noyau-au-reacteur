"""
sections_efficaces.py -- Pourquoi des neutrons LENTS : noyau compose, taille
quantique du neutron, resonances et loi en 1/v ; puis le partage des
fragments et l'origine de leur energie cinetique.

Donnees :
  * masses : AME2020 (donnees_liaison.py) ;
  * sections efficaces evaluees : JENDL-4.0 a 300 K (K. Shibata et al.,
    J. Nucl. Sci. Technol. 48 (2011) 1), lues dans les tables publiees par
    la JAEA (wwwndc.jaea.go.jp, « Tab80 ») et recopiees ci-dessous ;
  * resonance de l'U-238 a 6,67 eV : S.F. Mughabghab, Atlas of Neutron
    Resonances (Elsevier, 6e ed., 2018) ;
  * barrieres de reference : RIPL-3 (Capote et al. 2009), via la table 1 de
    Ryssens et al., Eur. Phys. J. A 59 (2023) 96.

Lancer :  python sections_efficaces.py
    ->  sections_taille.png, fission_seuil.png, fission_partage.png
    (les vraies courbes, reconstruites depuis ENDF/B-VIII.0 : sections_endf.py)
"""

import numpy as np
from donnees_liaison import table_liaison, energie_separation_neutron, chercher
from trace_termes import (SURFACE, ENCRE, ENCRE_2, DISCRET, GRILLE, AXE, SERIES,
                          style, titre, etiquettes)

HBARC = 197.327      # MeV fm
MN = 939.565         # masse du neutron, MeV
E2 = 1.44            # e^2, MeV fm
BARN = 100.0         # 1 b = 100 fm^2
E_TH = 0.0253        # eV, neutron "thermique" conventionnel (2200 m/s)

# JENDL-4.0 (300 K) : fission et capture a 0,0253 eV (b), moyenne de la
# fission sur le spectre des neutrons de fission (b), nu prompt et retarde
JENDL = {
    ("Th", 232): dict(sf=53.71e-6, sg=7.338, sf_rapide=79.94e-3, nup=1.851, nud=4.900e-2),
    ("U", 233): dict(sf=531.3, sg=45.26, sf_rapide=1.908, sg_rapide=58.96e-3, nup=2.478, nud=6.700e-3),
    ("U", 235): dict(sf=585.1, sg=98.71, sf_rapide=1.218, sg_rapide=87.41e-3, nup=2.420, nud=1.585e-2),
    ("U", 238): dict(sf=16.80e-6, sg=2.683, sf_rapide=306.4e-3, nup=2.281, nud=4.634e-2),
    ("Pu", 239): dict(sf=747.4, sg=271.5, sf_rapide=1.802, sg_rapide=52.61e-3, nup=2.872, nud=6.220e-3),
    ("Pu", 240): dict(sf=36.21e-3, sg=289.3, sf_rapide=1.328, nup=2.825, nud=9.110e-3),
    ("Pu", 241): dict(sf=1012.0, sg=363.1, sf_rapide=1.626, sg_rapide=71.97e-3, nup=2.931, nud=1.600e-2),
}
G_FISSION = {("U", 235): 0.977, ("Pu", 239): 1.058, ("U", 233): 0.997}

# Resonance s de l'U-238 (Mughabghab) : E0, largeur neutron (a E0), largeur gamma
U238_RES = dict(E0=6.674, Gn=1.48e-3, Gg=23.0e-3, g=1.0, A=238)

# Barrieres primaires recommandees (RIPL-3) pour les noyaux composes pair-pair
BARRIERE = {("U", 234): 5.50, ("U", 236): 5.67, ("Pu", 240): 6.05, ("Pu", 242): 5.85}

ZSYM = {90: "Th", 92: "U", 94: "Pu"}


def titre_section(t):
    print()
    print("=" * 76)
    print(t)
    print("=" * 76)


# ==========================================================================
# 1. Le noyau compose
# ==========================================================================
def noyau_compose(t):
    titre_section("1. LE NOYAU COMPOSE : CE QU'APPORTE LE NEUTRON")
    print("""
  Le neutron absorbe forme un noyau compose (Bohr 1936), excite de
      E* = S_n(noyau compose) + E_n * A/(A+1)
  Pour un neutron thermique, E_n = 0,025 eV : tout vient de S_n.
""")
    print(f"  {'cible':<8} {'N cible':>7} {'compose':>8} {'S_n (MeV)':>10} {'barriere':>9} "
          f"{'marge':>7} {'sigma_f thermique':>18}")
    lignes = []
    for (s, a), d in JENDL.items():
        z = {"Th": 90, "U": 92, "Pu": 94}[s]
        sn = energie_separation_neutron(t, s, a + 1)
        b = BARRIERE.get((s, a + 1))
        marge = f"{sn - b:+.2f}" if b else "  -"
        bt = f"{b:.2f}" if b else "  -"
        sf = d["sf"]
        sft = f"{sf:.1f} b" if sf >= 1 else (f"{sf * 1e3:.2f} mb" if sf >= 1e-3 else f"{sf * 1e6:.2f} ub")
        print(f"  {s}-{a:<5} {a - z:>7} {s}-{a + 1:<5} {sn:>10.3f} {bt:>9} {marge:>7} {sft:>18}")
        lignes.append((s, a, a - z, sn, sf))
    print("""
  Cible a N IMPAIR -> compose a N pair -> S_n grand (6,3 a 6,8 MeV) -> fissile.
  Cible a N PAIR   -> compose a N impair -> S_n petit (4,8 a 5,2 MeV) -> non.
""")
    print("  La signature de l'appariement dans les masses : S_n le long des uranium")
    print(f"  {'noyau':<8} {'N':>4} {'S_n (MeV)':>10}")
    for a in range(233, 241):
        sn = energie_separation_neutron(t, "U", a)
        print(f"  U-{a:<6} {a - 92:>4} {sn:>10.3f}  {'(N pair)' if (a - 92) % 2 == 0 else ''}")
    s_pair = [energie_separation_neutron(t, "U", a) for a in (234, 236, 238)]
    s_imp = [energie_separation_neutron(t, "U", a) for a in (235, 237, 239)]
    print(f"\n  moyenne N pair {np.mean(s_pair):.3f} MeV, N impair {np.mean(s_imp):.3f} MeV, "
          f"ecart {np.mean(s_pair) - np.mean(s_imp):.2f} MeV")
    print(f"  (gap d'appariement empirique 12/sqrt(A) = {12 / np.sqrt(236):.2f} MeV ; l'ecart pair-impair de S_n vaut ~2 Delta)")
    return lignes


# ==========================================================================
# 2. La taille quantique du neutron
# ==========================================================================
def lambda_reduite(E_eV, A):
    """Longueur d'onde reduite du neutron dans le centre de masse, en fm."""
    p = np.sqrt(2 * MN * np.asarray(E_eV) * 1e-6)          # MeV/c
    return HBARC / p * (A + 1) / A


def taille():
    titre_section("2. LA TAILLE APPARENTE D'UN NEUTRON LENT")
    A = 235
    R = 1.2 * A ** (1 / 3)
    print(f"\n  Rayon de l'U-235 : R = 1,2 A^(1/3) = {R:.2f} fm ; aire geometrique pi R^2 = "
          f"{np.pi * R * R / BARN:.2f} b\n")
    print(f"  {'E neutron':>12} {'lambda reduite':>16} {'pi lambda^2':>14} {'pi lambda^2 / pi R^2':>22}")
    for E, nom in ((E_TH, "0,0253 eV"), (1.0, "1 eV"), (1e3, "1 keV"), (1e6, "1 MeV"), (2e6, "2 MeV")):
        lam = lambda_reduite(E, A)
        aire = np.pi * lam ** 2 / BARN
        print(f"  {nom:>12} {lam:>13.1f} fm {aire:>12.4g} b {aire / (np.pi * R * R / BARN):>20.4g}")


# ==========================================================================
# 3. Breit-Wigner et loi en 1/v
# ==========================================================================
def breit_wigner_capture(E, E0, Gn, Gg, g, A):
    """Capture (n,gamma) d'une resonance s isolee, en barns (formule a un
    niveau de Breit et Wigner, 1936). La largeur neutron croit comme
    sqrt(E) : Gn(E) = Gn(E0) sqrt(E/E0)."""
    E = np.asarray(E, float)
    Gn_E = Gn * np.sqrt(E / E0)
    G = Gn_E + Gg
    lam = lambda_reduite(E, A)
    return np.pi * lam ** 2 / BARN * g * Gn_E * Gg / ((E - E0) ** 2 + G ** 2 / 4)


def resonances():
    titre_section("3. RESONANCES ET LOI EN 1/v")
    r = U238_RES
    pic = breit_wigner_capture(r["E0"], **r)
    th = breit_wigner_capture(E_TH, **r)
    jendl = JENDL[("U", 238)]["sg"]
    print(f"""
  Resonance s de l'U-238 : E0 = {r['E0']} eV, Gamma_n = {r['Gn'] * 1e3:.2f} meV,
  Gamma_gamma = {r['Gg'] * 1e3:.1f} meV (Mughabghab).
    capture au pic          : {pic:,.0f} b
    capture a 0,0253 eV     : {th:.3f} b   (cette seule resonance)
    JENDL-4.0 a 0,0253 eV   : {jendl:.3f} b
    -> la resonance a 6,67 eV explique {100 * th / jendl:.0f} % de la capture thermique ;
       le reste vient des autres resonances (et d'un niveau lie, sous le seuil).

  Loin sous la resonance (E << E0), (E - E0)^2 ~ E0^2 est constant,
  Gamma_n ~ sqrt(E) et pi lambda^2 ~ 1/E : sigma ~ sqrt(E)/E = 1/sqrt(E) ~ 1/v.""".replace(",", " "))
    print("\n  Verification sur le modele : sigma(E) * sqrt(E) doit etre constant")
    for E in (0.001, 0.01, 0.0253, 0.1, 0.5):
        print(f"    E = {E:<7} eV : sigma = {breit_wigner_capture(E, **r):8.3f} b, "
              f"sigma*sqrt(E) = {breit_wigner_capture(E, **r) * np.sqrt(E):.4f}")
    print("""
  Facteur de Westcott g (JENDL-4.0) : moyenne maxwellienne / sigma(0,0253 eV),
  qui vaut 1 pour une loi en 1/v exacte :""")
    for (s, a), g in G_FISSION.items():
        print(f"    fission {s}-{a} : g = {g:.3f}")
    print("""  Pu-239 s'ecarte nettement de 1/v (g > 1) : sa grande resonance de fission
  a 0,3 eV est deja dans la queue du spectre thermique.""")


def lent_rapide():
    titre_section("4. LENT CONTRE RAPIDE")
    print(f"\n  {'cible':<8} {'sigma_f thermique':>18} {'sigma_f rapide':>15} {'rapport':>10} "
          f"{'nu':>6} {'eta thermique':>14} {'eta rapide':>11}")
    for (s, a), d in JENDL.items():
        nu = d["nup"] + d["nud"]
        eta = nu * d["sf"] / (d["sf"] + d["sg"])
        er = (f"{nu * d['sf_rapide'] / (d['sf_rapide'] + d['sg_rapide']):>11.3f}"
              if "sg_rapide" in d else f"{'-':>11}")
        sft = f"{d['sf']:.1f} b" if d["sf"] >= 1 else (f"{d['sf'] * 1e3:.2f} mb" if d["sf"] >= 1e-3 else f"{d['sf'] * 1e6:.2f} ub")
        print(f"  {s}-{a:<5} {sft:>18} {d['sf_rapide']:>13.3f} b {d['sf'] / d['sf_rapide']:>10.3g} "
              f"{nu:>6.3f} {eta:>14.3f} {er}")
    print("""
  sigma_f rapide : moyenne sur le spectre des neutrons de fission de l'U-235
  (energie moyenne ~2 MeV), JENDL-4.0.
  eta = nu sigma_f / (sigma_f + sigma_gamma) : neutrons emis par neutron absorbe.
  eta rapide : sections moyennees sur le spectre de fission, avec le nu
  THERMIQUE (nu croit avec l'energie du neutron : c'est une borne basse).""")


# ==========================================================================
# 4 bis. Franchir la barriere
# ==========================================================================
def franchissement(t):
    titre_section("4 bis. FRANCHIR LA BARRIERE : EFFET TUNNEL ET SEUIL")
    print("""
  Rapport sigma_f / sigma_gamma a 0,0253 eV (JENDL-4.0) : la probabilite
  relative, pour le noyau compose, de fissionner plutot que d'emettre un gamma.
""")
    for (s, a), d in JENDL.items():
        sn = energie_separation_neutron(t, s, a + 1)
        print(f"    {s}-{a:<5} S_n = {sn:.3f} MeV   sigma_f/sigma_gamma = {d['sf'] / d['sg']:.3g}")
    print("""
  Penetrabilite de Hill-Wheeler T = 1/(1 + exp(2 pi (B - E*)/hbar_omega)) pour
  le noyau compose U-239 (cible U-238, neutron thermique, E* = S_n = 4,806 MeV).
  Sa barriere n'est pas dans notre table : on balaie la bande des actinides
  pair-pair voisins (RIPL-3, 5,5 a 6,3 MeV).
""")
    sn = energie_separation_neutron(t, "U", 239)
    print(f"  {'B (MeV)':>8} " + " ".join(f"{'hw=' + str(h):>11}" for h in (0.5, 0.75, 1.0))
          + f"   {'E_n pour T = 1/2':>17}")
    for B in (5.5, 6.0, 6.3):
        Ts = [1 / (1 + np.exp(2 * np.pi * (B - sn) / h)) for h in (0.5, 0.75, 1.0)]
        En = (B - sn) * 239 / 238
        print(f"  {B:>8.1f} " + " ".join(f"{x:>11.2e}" for x in Ts) + f"   {En:>13.2f} MeV")
    print("""
  La penetrabilite va de 1e-8 a 1e-2 selon B et hbar_omega : on ne peut pas
  en tirer une valeur precise. Mais le rapport mesure (6e-6) a l'ordre de
  grandeur d'un passage par effet tunnel sous une barriere d'environ 6 MeV.
  Pour la franchir par-dessus, il faut un neutron de 0,7 a 1,5 MeV : c'est
  le seuil observe, ~1 MeV.""")


# ==========================================================================
# 5. Le partage des fragments
# ==========================================================================
def partage(t, Zc=92, Ac=236):
    titre_section("5. LE PARTAGE : ENERGIE DISPONIBLE ET REPULSION COULOMBIENNE")
    Bc = float(chercher(t, ZSYM[Zc], Ac)["B_MeV"])
    res = []
    for AL in range(70, Ac // 2 + 1):
        AH = Ac - AL
        meilleur = None
        for ZL in range(25, Zc // 2 + 1 + 6):
            ZH = Zc - ZL
            m = (t["Z"] == ZL) & (t["A"] == AL)
            h = (t["Z"] == ZH) & (t["A"] == AH)
            if m.any() and h.any() and not t["estimee"][m][0] and not t["estimee"][h][0]:
                Q = float(t["B_MeV"][m][0] + t["B_MeV"][h][0] - Bc)
                if meilleur is None or Q > meilleur[0]:
                    meilleur = (Q, ZL, ZH)
        if meilleur:
            Q, ZL, ZH = meilleur
            r0 = 1.2
            EC = ZL * ZH * E2 / (r0 * (AL ** (1 / 3) + AH ** (1 / 3)))
            res.append((AL, AH, ZL, ZH, Q, EC))
    res = np.array(res)
    print(f"\n  U-236 -> deux fragments (sans neutron), masses mesurees AME2020 ;")
    print("  pour chaque A leger, la charge qui maximise Q.\n")
    print(f"  {'A leger':>8} {'A lourd':>8} {'Z':>8} {'Q (MeV)':>9} {'Coulomb au contact (MeV)':>26}")
    for row in res:
        if int(row[0]) in (80, 90, 95, 100, 104, 110, 118):
            print(f"  {int(row[0]):>8} {int(row[1]):>8} {int(row[2]):>3}/{int(row[3]):<4} {row[4]:>9.1f} {row[5]:>26.1f}")
    i = int(np.argmax(res[:, 4]))
    k = int(np.argmin(np.abs(res[:, 0] - Ac / 2)))
    print(f"\n  Q maximal pour A leger = {int(res[i, 0])}, A lourd = {int(res[i, 1])} "
          f"(Z = {int(res[i, 2])}/{int(res[i, 3])}) : {res[i, 4]:.1f} MeV")
    print(f"  partage symetrique {int(res[k, 0])}/{int(res[k, 1])} : {res[k, 4]:.1f} MeV")
    print("  -> le maximum n'est pas au partage symetrique mais au fragment lourd Sn-132,")
    print("     doublement magique (Z = 50, N = 82) : les couches se lisent deja dans Q.")
    TKE = 168.0
    j = int(np.argmin(np.abs(res[:, 0] - 95)))
    ZL, ZH, AL, AH = res[j, 2], res[j, 3], res[j, 0], res[j, 1]
    d = ZL * ZH * E2 / TKE
    contact = 1.2 * (AL ** (1 / 3) + AH ** (1 / 3))
    print(f"""
  Energie cinetique des fragments ~ {TKE:.0f} MeV (bilan de fission.py) : c'est
  la repulsion coulombienne des deux fragments, liberee a la scission.
  Pour le partage A = {int(AL)}/{int(AH)} (Z = {int(ZL)}/{int(ZH)}) :
    spheres au contact (r0 = 1,2 fm) : distance des centres {contact:.1f} fm,
                                       Coulomb {res[j, 5]:.0f} MeV
    distance qui donne {TKE:.0f} MeV      : {d:.1f} fm
  -> a la scission, les centres sont ~{d - contact:.0f} fm plus loin qu'au contact :
     les fragments sont ALLONGES, relies par un col. Leur energie de
     deformation, rendue apres la rupture, chauffe les fragments, qui
     l'evacuent en neutrons prompts et en gammas.""")
    return res


# ==========================================================================
# Figures
# ==========================================================================
def figures(lignes, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter

    def figure(gauche=0.11, droite=0.80):
        fig, ax = plt.subplots(figsize=(9, 5.6), facecolor=SURFACE)
        fig.subplots_adjust(left=gauche, right=droite, top=0.83, bottom=0.15)
        return fig, ax

    def note(fig, texte):
        fig.text(0.02, 0.02, texte, fontsize=8.5, color=DISCRET)

    def fmt_e(v, _):
        if v >= 1e6:
            return f"{v / 1e6:g} MeV"
        if v >= 1e3:
            return f"{v / 1e3:g} keV"
        return f"{v:g} eV".replace(".", ",")

    from matplotlib.ticker import FuncFormatter, FixedLocator

    def puissance(v, _):
        if v <= 0:
            return ""
        e = int(round(np.log10(v)))
        if abs(np.log10(v) - e) > 1e-6:
            return ""
        if -1 <= e <= 3:
            return f"{v:g}".replace(".", ",")
        return f"10$^{{{e}}}$"

    # --- 1. taille apparente
    fig, ax = figure(droite=0.74)
    style(ax)
    E = np.logspace(-3, 7, 400)
    A = 235
    aire = np.pi * lambda_reduite(E, A) ** 2 / BARN
    geo = np.pi * (1.2 * A ** (1 / 3)) ** 2 / BARN
    ax.plot(E, aire, color=SERIES[0], lw=2)
    ax.axhline(geo, color=SERIES[1], lw=2)
    sf = JENDL[("U", 235)]
    ax.plot([E_TH, 2e6], [sf["sf"], sf["sf_rapide"]], "s", ms=8, color=SERIES[2],
            mec=SURFACE, mew=2, zorder=6)
    for e, lab, dxy, ha in ((E_TH, "neutron thermique (0,0253 eV)", (12, 4), "left"),
                            (2e6, "neutron de fission (2 MeV)", (-12, -18), "right")):
        a = np.pi * lambda_reduite(e, A) ** 2 / BARN
        ax.plot(e, a, "o", ms=8, color=SERIES[0], mec=SURFACE, mew=2, zorder=5)
        txt = (f"{a:.2g}".replace("e+0", " × 10^").replace(".", ",") + " b") if a > 1e4 else f"{a:.2f} b".replace(".", ",")
        txt = txt.replace(" × 10^7", " × 10⁷")
        ax.annotate(f"{lab} : π ƛ² = {txt}", (e, a), xytext=dxy, textcoords="offset points",
                    fontsize=10, color=ENCRE, ha=ha)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 1e7)
    ax.set_ylim(0.05, 3e8)
    ax.xaxis.set_major_locator(FixedLocator([1e-3, 1e-1, 1e1, 1e3, 1e5, 1e7]))
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_e))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]))
    ax.yaxis.set_major_formatter(FuncFormatter(puissance))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("énergie du neutron")
    ax.set_ylabel("aire (barns)")
    etiquettes(ax, None, [float(aire[-1]), geo, float(sf["sf_rapide"])],
               ["π ƛ² : taille quantique", "π R² : taille de l'U-235", "σ_f mesuré de l'U-235"],
               SERIES[:3], ecart=0.55)
    titre(ax, "Un neutron lent « voit » un noyau des millions de fois plus gros",
          "La longueur d'onde réduite ƛ = ħ/p croît quand le neutron ralentit ; le noyau, lui, ne change pas.")
    note(fig, "ƛ dans le centre de masse, U-235 ; R = 1,2 A^(1/3) fm. Carrés : fission de l'U-235 à 0,0253 eV "
              "et moyennée sur le spectre de fission (JENDL-4.0).")
    fig.savefig("sections_taille.png", dpi=120, facecolor=SURFACE)
    plt.close(fig)

    # --- 3. seuil de fissilite
    fig, ax = figure(gauche=0.13, droite=0.95)
    style(ax)
    ax.axvspan(5.5, 6.3, color=GRILLE, alpha=0.6, lw=0)
    ax.text(5.9, 4e-4, "barrières\nRIPL-3\n5,5 à 6,3 MeV", ha="center", fontsize=10, color=ENCRE_2)
    decalages = {("U", 235): (10, -12, "left"), ("Pu", 239): (-10, 10, "right"),
                 ("U", 238): (-10, 0, "right"), ("Pu", 241): (-10, 8, "right")}
    for s, a, n, sn, sfv in lignes:
        col = SERIES[0] if n % 2 else SERIES[1]
        ax.plot(sn, sfv, "o", ms=10, color=col, mec=SURFACE, mew=2, zorder=5)
        dx, dy, ha = decalages.get((s, a), (10, 0, "left"))
        ax.annotate(f"{s}-{a}", (sn, sfv), xytext=(dx, dy), textcoords="offset points",
                    fontsize=10.5, color=ENCRE, va="center", ha=ha)
    ax.set_yscale("log")
    ax.set_xlim(4.5, 7.1)
    ax.set_ylim(3e-6, 2e4)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax.yaxis.set_major_formatter(FuncFormatter(
        lambda v, _: (f"{v * 1e6:g} μb" if v < 1e-3 else (f"{v * 1e3:g} mb" if v < 1 else f"{v:g} b"))))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("énergie apportée par un neutron lent : S_n du noyau composé (MeV)")
    ax.set_ylabel("fission par neutron thermique")
    ax.plot([], [], "o", color=SERIES[0], ms=9, label="cible à N impair")
    ax.plot([], [], "o", color=SERIES[1], ms=9, label="cible à N pair")
    ax.legend(loc="upper left", frameon=False, fontsize=10.5, labelcolor=ENCRE_2)
    titre(ax, "Fissile ou non : un saut de sept ordres de grandeur",
          "Au-dessus de la barrière, des centaines de barns ; au-dessous, des microbarns à quelques millibarns.")
    note(fig, "Abscisse : AME2020. Ordonnée : JENDL-4.0, fission à 0,0253 eV. Bande : barrières primaires "
              "RIPL-3 des actinides pair-pair (U-234 à U-238).")
    fig.savefig("fission_seuil.png", dpi=120, facecolor=SURFACE)
    plt.close(fig)

    # --- 4. partage
    fig, ax = figure(droite=0.72)
    style(ax)
    AL, Q, EC = res[:, 0], res[:, 4], res[:, 5]
    ax.plot(AL, EC, color=SERIES[1], lw=2)
    ax.plot(AL, Q, color=SERIES[0], lw=2)
    ax.axhline(168, color=SERIES[2], lw=2)
    i = int(np.argmax(Q))
    ax.plot(AL[i], Q[i], "o", ms=8, color=SERIES[0], mec=SURFACE, mew=2, zorder=5)
    ax.annotate(f"¹³²Sn + ¹⁰⁴Mo : {Q[i]:.0f} MeV", (AL[i], Q[i]), xytext=(-8, 10),
                textcoords="offset points", fontsize=10, color=ENCRE, ha="right")
    ax.axvline(95, color=AXE, lw=0.8)
    ax.text(95.5, 150, "pic léger\nobservé\nA ≈ 95", fontsize=9.5, color=ENCRE_2)
    ax.set_xlim(AL.min() - 1, AL.max() + 1)
    ax.set_ylim(140, 270)
    ax.set_xlabel("masse du fragment léger (U-236 → deux fragments)")
    ax.set_ylabel("énergie (MeV)")
    etiquettes(ax, None, [float(EC[-1]), float(Q[-1]), 168.0],
               ["Coulomb au contact", "Q (masses mesurées)", "énergie cinétique\ndes fragments"],
               [SERIES[1], SERIES[0], SERIES[2]], ecart=11)
    titre(ax, "Le partage : énergie disponible et répulsion",
          "Q culmine avec ¹³²Sn, doublement magique ; au contact, la répulsion dépasserait Q.")
    note(fig, "Q : AME2020, pour chaque masse la charge la plus favorable. Coulomb : sphères de rayon "
              "1,2 A^(1/3) fm au contact. 168 MeV : bilan de fission.py.")
    fig.savefig("fission_partage.png", dpi=120, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    t = table_liaison()
    lignes = noyau_compose(t)
    taille()
    resonances()
    lent_rapide()
    franchissement(t)
    res = partage(t)
    figures(lignes, res)
    print("\n  --> sections_taille.png, fission_seuil.png, fission_partage.png")
