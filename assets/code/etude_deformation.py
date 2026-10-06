"""
etude_deformation.py -- Pourquoi les noyaux se deforment : les pieces du
modele macroscopique-microscopique, regardees une a une et confrontees aux
mesures.

  1. schema des niveaux spheriques : oscillateur, + mu l^2, + spin-orbite
  2. pentes du diagramme de Nilsson : formule au premier ordre
  3. diagramme de Nilsson des actinides : Omega^pi, gaps deformes
  4. Strutinsky : l'escalier des niveaux et sa version lissee ;
     correction de couches en fonction du nombre de nucleons
  5. bilan d'Er-166 : goutte, couches, appariement (effet Jahn-Teller)
  6. BCS : occupations (Sn-116) et gaps face aux masses (Delta(3))
  7. transition de forme des samariums, face aux beta2 et E(2+) mesures
  8. moments quadrupolaires : le modele face aux B(E2) mesures
  9. moment d'inertie de l'U-238 : une signature de l'appariement

Donnees :
  * masses AME2020 (donnees_liaison.py) ;
  * B(E2; 0+ -> 2+) et E(2+) adoptes, NNDC : B. Pritychenko et al., At. Data
    Nucl. Data Tables 107 (2016) 1, et a defaut S. Raman et al., ADNDT 78
    (2001) 1. Extrait dans be2_adopte.csv (voir donnees_be2).

Prerequis : micro.npz (precalcul_micro.py).
Lancer :    python etude_deformation.py   ->  deformation_*.png  (~1 min)
"""

import csv
import json
import os
import numpy as np
from donnees_liaison import table_liaison
from modeles_masse import ModeleLineaire, colonnes_M2
from nilsson import (base, _blocs_constants, _facteur_volume, kappa_mu,
                     spectre, hbar_omega)
from strutinsky import Lissage, appariement_details
from frdm import ModeleMacMic, GRILLE_DELTA, rapport_axes
from trace_termes import (SURFACE, ENCRE, ENCRE_2, DISCRET, GRILLE, AXE, SERIES,
                          BLEUS, style, titre as titre_axe, ecarter, etiquettes)

LETTRES = "spdfghijklmnoqr"          # l = 0 a 14 ; l = 7 note j, usage nucleaire
SYMBOLES = {50: "Sn", 62: "Sm", 68: "Er", 82: "Pb", 92: "U"}
URL_BE2 = "https://www.nndc.bnl.gov/be2/data/adopted-entries.json"


def titre(t):
    print("\n" + "=" * 76 + "\n" + t + "\n" + "=" * 76)


def virgule(x, nd=2):
    return f"{x:.{nd}f}".replace(".", ",").replace("-", "−")


# ==========================================================================
# Donnees B(E2)
# ==========================================================================
def beta2_de_be2(BE2, Z, A):
    """beta2 = (4 pi / 3 Z R0^2) sqrt(B(E2)up / e^2), R0 = 1,2 A^(1/3) fm :
    la convention des deux compilations. B(E2) en e^2 b^2 (1 b = 100 fm^2)."""
    R0_carre = (1.2 * A ** (1 / 3)) ** 2 / 100.0
    return 4 * np.pi / (3 * Z * R0_carre) * np.sqrt(BE2)


def donnees_be2(chemin="be2_adopte.csv", source="adopted-entries.json"):
    """{(Z, N): dict(E2_keV, BE2, beta2, beta2_table, reference)}, noyaux pair-pair.

    Lit be2_adopte.csv s'il existe. Sinon le reconstruit a partir du fichier
    JSON des valeurs adoptees du NNDC (telecharge depuis URL_BE2 s'il manque).
    Quand un noyau figure dans les deux compilations, on garde Pritychenko
    et al. (2016), plus recente que Raman et al. (2001).

    beta2 est RECALCULE a partir de B(E2), la grandeur mesuree, avec la
    formule des tables (beta2_de_be2). La colonne beta2_table garde la valeur
    imprimee : elles different de plus de 2 % pour quelques noyaux, dont
    l'U-238 de la table de 2016 (0,2741 imprime, 0,2875 d'apres son B(E2))."""
    if not os.path.exists(chemin):
        if not os.path.exists(source):
            import urllib.request
            urllib.request.urlretrieve(URL_BE2, source)
        entrees = json.load(open(source))
        garde = {}
        for e in entrees:
            if "be2" not in e:
                continue
            z, n = e["nuclide"]["z"], e["nuclide"]["n"]
            recent = "ADNDT_107" in e["reference"]
            if (z, n) in garde and not recent:
                continue
            garde[(z, n)] = dict(
                Z=z, N=n, E2_keV=e.get("transitionEnergy", {}).get("value", ""),
                BE2=e["be2"]["value"], BE2_err=e["be2"].get("error", ""),
                beta2_table=e["deformationParameter"]["value"],
                reference="Pritychenko 2016" if recent else "Raman 2001")
        with open(chemin, "w", newline="") as f:
            f.write(f"# B(E2;0+->2+) adoptes, NNDC ({URL_BE2})\n")
            f.write("# Pritychenko et al. ADNDT 107 (2016) 1, sinon Raman et al. ADNDT 78 (2001) 1\n")
            f.write("# B(E2) en e^2 b^2, E(2+) en keV ; beta2_table : valeur imprimee "
                    "(R0 = 1,2 A^(1/3) fm)\n")
            w = csv.DictWriter(f, lineterminator="\n",
                               fieldnames=["Z", "N", "E2_keV", "BE2", "BE2_err",
                                           "beta2_table", "reference"])
            w.writeheader()
            for k in sorted(garde):
                w.writerow(garde[k])
    table = {}
    with open(chemin) as f:
        lignes = [l for l in f if not l.startswith("#")]
    for r in csv.DictReader(lignes):
        z, n = int(r["Z"]), int(r["N"])
        table[(z, n)] = dict(
            E2_keV=float(r["E2_keV"]) if r["E2_keV"] else np.nan,
            BE2=float(r["BE2"]), beta2=float(beta2_de_be2(float(r["BE2"]), z, z + n)),
            beta2_table=float(r["beta2_table"]), reference=r["reference"])
    return table


def beta2_equivalent(delta):
    """Le beta2 que les tables de B(E2) attribueraient a notre forme.

    Le modele donne un ellipsoide uniformement charge, de demi-axes
    a = R q^(2/3) (axe de symetrie) et b = R q^(-1/3), q = a/b. Son moment
    quadrupolaire intrinseque est exact : Q0 = (2/5) Z (a^2 - b^2). Les tables
    convertissent Q0 en beta2 par Q0 = (3 / sqrt(5 pi)) Z R0^2 beta2. Avec le
    meme R0, le rayon s'elimine :
        beta2_eq = (2 sqrt(5 pi) / 15) (q^(4/3) - q^(-2/3))
    Au premier ordre, beta2_eq = delta / 0,946 (delta_vers_beta2 de nilsson.py)."""
    q = rapport_axes(np.asarray(delta, float))
    return 2 * np.sqrt(5 * np.pi) / 15 * (q ** (4 / 3) - q ** (-2 / 3))


# ==========================================================================
# Outils Nilsson : spectre etiquete
# ==========================================================================
def nom_orbitale(N, l, jx2):
    """Notation spectroscopique : 1i13/2, 2g9/2... (n = nombre de noeuds + 1)."""
    return f"{(N - l) // 2 + 1}{LETTRES[l]}{jx2}/2"


def lignes_nilsson(deltas, type_nucleon, N_max=14):
    """Niveaux suivis en continu : dans un bloc (Omega, parite), deux niveaux ne
    se croisent jamais (regle de non-croisement), donc le k-ieme niveau du bloc
    est une ligne continue. Retourne une liste de dict(E, om_x2, parite, parent)
    ou E[i] est l'energie a deltas[i] et parent l'orbitale spherique d'origine."""
    lignes = []
    for om_x2 in range(1, 2 * N_max + 2, 2):
        etats = base(N_max, om_x2)
        diag_osc, so, Q = _blocs_constants(N_max, om_x2, type_nucleon)
        par = np.array([N % 2 for N, _, _ in etats])
        for p in (0, 1):
            idx = np.where(par == p)[0]
            if len(idx) == 0:
                continue
            E = []
            parent = None
            for d in deltas:
                c = _facteur_volume(d)
                H = (c * (np.diag(diag_osc[idx]) - (2 / 3) * d * Q[np.ix_(idx, idx)])
                     + np.diag(so[idx]))
                w, v = np.linalg.eigh(H)
                E.append(w)
                if abs(d) < 1e-12:
                    parent = [nom_orbitale(*etats[idx[np.argmax(v[:, k] ** 2)]])
                              for k in range(len(w))]
            E = np.array(E)
            for k in range(E.shape[1]):
                lignes.append(dict(E=E[:, k], om_x2=om_x2, parite="+" if p == 0 else "−",
                                   parent=parent[k] if parent else ""))
    return lignes


# ==========================================================================
def etape1_niveaux():
    titre("1. LES NIVEAUX SPHERIQUES : D'OU VIENNENT LES NOMBRES MAGIQUES")
    print("""
  Trois etages, neutrons, en unites de hbar*w0 :
    oscillateur pur       e = N + 3/2              (couches tres degenerees)
    + mu (l^2 - <l^2>_N)  separe les l d'une meme couche N
    + 2 l.s               separe j = l +- 1/2, d'un ecart kappa (2l+1)
""")
    niveaux = []
    for N in range(0, 7):
        kap, mu = kappa_mu("neutron", N)
        for l in range(N % 2, N + 1, 2):
            e_l = N + 1.5 - kap * mu * (l * (l + 1) - N * (N + 3) / 2)
            for jx2 in (2 * l - 1, 2 * l + 1):
                if jx2 < 0:
                    continue
                j = jx2 / 2
                ls = (j * (j + 1) - l * (l + 1) - 0.75) / 2
                niveaux.append(dict(N=N, l=l, jx2=jx2, e_osc=N + 1.5, e_l=e_l,
                                    e=e_l - 2 * kap * ls, deg=jx2 + 1))
    niveaux.sort(key=lambda x: x["e"])
    cumul = 0
    print(f"  {'orbitale':<10} {'e (hw0)':>8} {'places':>7} {'cumul':>6} {'ecart au suivant':>17}")
    for i, v in enumerate(niveaux):
        cumul += v["deg"]
        v["cumul"] = cumul
        v["ecart"] = niveaux[i + 1]["e"] - v["e"] if i + 1 < len(niveaux) else 0.0
        if v["N"] <= 6:
            marque = "  <== couche fermee" if v["ecart"] > 0.25 else ""
            print(f"  {nom_orbitale(v['N'], v['l'], v['jx2']):<10} {v['e']:>8.3f} "
                  f"{v['deg']:>7} {cumul:>6} {v['ecart']:>17.3f}{marque}")
    print("""
  Oscillateur pur : couches fermees a 2, 8, 20, 40, 70, 112, 168.
  Avec le spin-orbite, l'orbitale de j maximal de chaque couche N >= 3
  (1f7/2, 1g9/2, 1h11/2, 1i13/2) descend dans la couche du dessous :
  apparaissent 28, 50, 82, 126. Les ecarts a 6 et 16 sont reels dans ce
  jeu de parametres (kappa = 0,105 pour N <= 2), mais ne concernent que des
  noyaux tres legers, hors du domaine du modele.
""")
    return niveaux


# ==========================================================================
def etape2_pentes():
    titre("2. LES PENTES DU DIAGRAMME DE NILSSON")
    print("""
  Au premier ordre en delta (theoreme de Hellmann-Feynman), un niveau
  |N l j Omega> se deplace de -(2/3) delta <rho^2 P2> :
      de/ddelta = -(2/3) (N + 3/2) [j(j+1) - 3 Omega^2] / [4 j(j+1)]   (hbar*w0)
  Comparaison a la diagonalisation complete (difference finie, delta = 0,01) :
""")
    lignes = lignes_nilsson([0.0, 0.01], "neutron")
    res = []
    print(f"  {'orbitale':<10} {'Omega':>6} {'formule':>9} {'calcul':>9}")
    for nom, N, jx2 in (("1i13/2", 6, 13),):
        j = jx2 / 2
        for om_x2 in range(1, jx2 + 1, 2):
            om = om_x2 / 2
            f = -(2 / 3) * (N + 1.5) * (j * (j + 1) - 3 * om * om) / (4 * j * (j + 1))
            L = [x for x in lignes if x["parent"] == nom and x["om_x2"] == om_x2][0]
            c = (L["E"][1] - L["E"][0]) / 0.01
            res.append((om_x2, f, c))
            print(f"  {nom:<10} {om_x2:>4}/2 {f:>9.3f} {c:>9.3f}")
    print("""
  Accord a 2 % pres (l'ecart est l'effet du second ordre, melange des
  couches N et N +- 2). Omega petit : le nucleon orbite autour de l'axe
  long ; sa densite s'aligne avec la deformation allongee et son energie
  BAISSE. Omega = j : orbite equatoriale, energie qui MONTE. La somme des
  pentes d'une couche j pleine est nulle : une couche complete est
  spherique et ne gagne rien a se deformer au premier ordre.
""")
    return res


# ==========================================================================
def etape3_nilsson(Bd):
    titre("3. DIAGRAMME DE NILSSON DES ACTINIDES : LES GAPS DEFORMES")
    print("\n  Ecart (hbar*w0) entre le dernier niveau occupe et le suivant :")
    gaps = {}
    for typ, ns in (("neutron", (142, 144, 146, 148, 150, 152, 154)),
                    ("proton", (88, 90, 92, 94, 96, 98, 100, 102))):
        print(f"  {typ}s")
        for d in (0.0, 0.15, 0.2, 0.25, 0.3):
            e = spectre(d, typ)
            g = {n: e[n // 2] - e[n // 2 - 1] for n in ns}
            gaps[(typ, d)] = g
            print(f"    delta={d:4.2f} " + " ".join(f"{n}:{g[n]:.3f}" for n in ns))
    print("""
  A delta = 0, ces nombres tombent au milieu de couches tres degenerees :
  aucun ecart. Vers delta = 0,25 s'ouvrent des gaps DEFORMES, a N = 152
  (0,11 hbar*w0, soit ~0,7 MeV) et Z = 100 (0,15 hbar*w0). Ils restent 3 a 5
  fois plus petits que les gaps spheriques (0,4 a 0,6 hbar*w0 a 50, 82, 126).

  Le test dans les masses MESUREES (AME2020) : un gap se voit comme un saut
  de l'energie de separation de deux neutrons S2n. On calcule
      delta2n(N) = S2n(N) - S2n(N+2),   S2n(N) = B(N) - B(N-2)
  qui culmine a une fermeture de couche :""")
    for z, nom in ((96, "Cm"), (98, "Cf"), (100, "Fm")):
        ligne = f"    {nom} (Z = {z}) :"
        for n in range(144, 156, 2):
            k = [(z, n - 2), (z, n), (z, n + 2)]
            if all(x in Bd for x in k):
                d2n = (Bd[k[1]] - Bd[k[0]]) - (Bd[k[2]] - Bd[k[1]])
                ligne += f"  {n}:{d2n:.2f}"
        print(ligne + "  MeV")
    print("""
  Le pic a N = 152, net pour Fm et visible pour Cf, confirme le gap deforme
  predit. Z = 100 ne peut pas etre teste ici de la meme facon : il faudrait
  les masses de Z = 102, et la table de donnees_liaison.py s'arrete a Z = 100.
""")
    return gaps


# ==========================================================================
def etape4_strutinsky():
    titre("4. STRUTINSKY : CORRECTION DE COUCHES EN FONCTION DE n")
    res = {}
    ns = np.arange(8, 186, 2)
    for d in (0.0, 0.25):
        L = Lissage(spectre(d, "neutron"), 1.2)
        res[d] = (ns, np.array([L.correction_couches(int(n)) for n in ns]))
    print(f"\n  neutrons, correction de couches (hbar*w0), extrema locaux :")
    for d, (n, dE) in res.items():
        mins = [(int(n[k]), dE[k]) for k in range(1, len(n) - 1)
                if dE[k] < dE[k - 1] and dE[k] < dE[k + 1] and dE[k] < -0.15]
        maxs = [(int(n[k]), dE[k]) for k in range(1, len(n) - 1)
                if dE[k] > dE[k - 1] and dE[k] > dE[k + 1] and dE[k] > 0.4]
        print(f"    delta={d:4.2f}  minima " + ", ".join(f"{a}:{b:+.2f}" for a, b in mins))
        print(f"               maxima " + ", ".join(f"{a}:{b:+.2f}" for a, b in maxs))
    hw = hbar_omega(82, 126, "neutron")
    m = res[0.0][1][list(res[0.0][0]).index(126)]
    print(f"""
  Spherique : minima aux nombres magiques, maxima a mi-couche (n = 64, 100,
  150). Deforme (delta = 0,25) : les minima se deplacent a mi-couche
  spherique (n ~ 100-106, 148-152). La ou la sphere est penalisee, une
  forme allongee est favorisee : c'est tout le mecanisme de la deformation.
  Ordre de grandeur : a N = 126, {m:.2f} hbar*w0 = {m * hw:.1f} MeV (hbar*w0 = {hw:.2f} MeV
  pour les neutrons de Pb-208).
""")
    L = Lissage(spectre(0.0, "neutron"), 1.2)
    return res, L


# ==========================================================================
def ajuster():
    t = table_liaison(mesurees_seulement=True)
    t = t[(t["Z"] >= 8) & (t["N"] >= 8)]
    Z, N, B = t["Z"].astype(float), t["N"].astype(float), t["B_MeV"]
    c0 = ModeleLineaire("M2", colonnes_M2).ajuster(Z, N, B).coefs
    mm = ModeleMacMic(appariement=True).ajuster(Z, N, B, c0)
    mm0 = ModeleMacMic(appariement=False)
    mm0.coefs = mm.coefs
    return Z, N, B, mm, mm0


def etape5_bilan(mm, mm0):
    titre("5. BILAN D'ER-166 : QUI VEUT LA DEFORMATION, QUI S'Y OPPOSE")
    z, n = 68, 98
    tot, mac = mm.energie_forme(z, n)
    tot0, _ = mm0.energie_forme(z, n)
    couches = tot0 - mac
    appar = tot - tot0
    G = list(GRILLE_DELTA)
    print(f"\n  {'delta':>7} {'goutte':>8} {'couches':>8} {'appar.':>8} {'total':>8}   (MeV)")
    for d in (-0.1, -0.05, -0.025, 0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3):
        k = G.index(round(d, 4))
        print(f"  {d:>+7.3f} {mac[k]:>8.2f} {couches[k]:>8.2f} {appar[k]:>8.2f} {tot[k]:>8.2f}")
    k0, k1, k25 = G.index(0.0), G.index(0.025), G.index(0.25)
    pente = (couches[k1] - couches[k0]) / 0.025
    pente_m = (couches[k0 - 1] - couches[k0]) / 0.025
    c = mm.coefs
    A, I = z + n, n - z
    Es = c[1] * A ** (2 / 3) - c[5] * I ** 2 / A ** (4 / 3)
    Ec = c[2] * z * (z - 1) / A ** (1 / 3) - c[6] * z ** (4 / 3) / A ** (1 / 3)
    x = Ec / (2 * Es)
    C_th = 8 / 45 * Es * (1 - x)
    C_num = mac[G.index(0.1)] / 0.1 ** 2
    print(f"""
  * La goutte : E ~ C delta^2, C = (8/45) E_s (1 - x) = {C_th:.1f} MeV
    (calcul exact a delta = 0,1 : {C_num:.1f} MeV). Elle prefere la sphere.
  * Les couches seules : la sphere est un SOMMET. En partant de delta = 0,
    l'energie de couches baisse LINEAIREMENT des deux cotes ({pente:.0f} et
    {pente_m:.0f} MeV par unite de delta) : une couche partiellement remplie
    leve sa degenerescence en se deformant. C'est l'effet Jahn-Teller (1937),
    transpose au noyau. Contre une goutte quadratique, un gain lineaire
    l'emporte toujours a petite deformation.
  * L'appariement : maximal ({appar[k0]:.2f} MeV) a la sphere, ou la densite de
    niveaux au niveau de Fermi est la plus forte ; il s'efface en se
    deformant ({appar[k25]:+.2f} MeV a delta = 0,25). Il S'OPPOSE a la deformation
    et arrondit le sommet en un plateau.
  * Bilan a delta = 0,25 : couches {couches[k25] - couches[k0]:+.1f}, appariement
    {appar[k25] - appar[k0]:+.1f}, goutte {mac[k25]:+.1f} -> total {tot[k25] - tot[k0]:+.1f} MeV.
""")
    return dict(mac=mac, couches=couches, appar=appar, tot=tot)


# ==========================================================================
def delta_3(Bd, z, n):
    """Difference de masse pair-impair a trois points, en neutrons :
        Delta3(N) = -((-1)^N / 2) [B(N+1) - 2 B(N) + B(N-1)]   (B > 0)
    A N impair, c'est la mesure standard du gap d'appariement des neutrons
    (Satula, Dobaczewski, Nazarewicz 1998) ; a N pair, elle melange gap et
    ecart entre niveaux individuels."""
    k = [(z, n - 1), (z, n), (z, n + 1)]
    if not all(x in Bd for x in k):
        return None
    return -((-1) ** n) / 2 * (Bd[k[2]] - 2 * Bd[k[1]] + Bd[k[0]])


def details_bcs(mm, z, n, typ="neutron"):
    d = mm.forme(z, n)[0]
    dg = GRILLE_DELTA[np.argmin(np.abs(GRILLE_DELTA - d))]
    L = Lissage(spectre(dg, typ), 1.2)
    hw = hbar_omega(z, n, typ)
    return appariement_details(L, n if typ == "neutron" else z, hw, z + n), dg


def etape6_bcs(Z, N, B, mm):
    titre("6. APPARIEMENT : OCCUPATIONS ET GAPS FACE AUX MASSES")
    r, dg = details_bcs(mm, 50, 66)
    print(f"""
  Sn-116, neutrons (delta = {dg:g}) : Delta = {r['Delta']:.3f} MeV, G = {r['G']:.3f} MeV,
  {len(r['eps'])} niveaux de paires dans la fenetre (+- hbar*w0 autour du niveau de Fermi).
  Les occupations v^2 passent de 1 a 0 sur une largeur ~ 2 Delta autour de
  lambda : la surface de Fermi est diffuse.
""")
    occ = r
    Bd = {(int(z), int(n)): b for z, n, b in zip(Z, N, B)}
    chaines = {}
    for z, plage in ((50, range(50, 87)), (82, range(100, 133))):
        mod, exp = [], []
        for n in plage:
            if n % 2 == 0 and (z, n) in Bd:
                mod.append((n, details_bcs(mm, z, n)[0]["Delta"]))
            elif n % 2 == 1:
                e = delta_3(Bd, z, n)
                if e is not None:
                    exp.append((n, e))
        chaines[z] = (np.array(mod), np.array(exp))
        m, e = chaines[z]
        print(f"  {SYMBOLES[z]} (Z = {z}) :  N pair -> Delta_n BCS ; N impair -> Delta3 mesure")
        ligne = "   "
        for n_, v in sorted(list(map(tuple, m)) + list(map(tuple, e))):
            ligne += f" {int(n_)}:{v:.2f}"
            if len(ligne) > 70:
                print(ligne)
                ligne = "   "
        if ligne.strip():
            print(ligne)
    # accord quantitatif, hors couches fermees
    print()
    for z, (m, e) in chaines.items():
        dm = {int(a): b for a, b in m}
        paires = [(dm[int(n) - 1], dm[int(n) + 1], v) for n, v in e
                  if int(n) - 1 in dm and int(n) + 1 in dm
                  and dm[int(n) - 1] > 0 and dm[int(n) + 1] > 0]
        moy = np.array([(a + b) / 2 for a, b, _ in paires])
        ex = np.array([v for _, _, v in paires])
        print(f"  {SYMBOLES[z]} : {len(ex)} isotopes impairs ; moyenne Delta BCS (voisins) "
              f"{moy.mean():.2f}, Delta3 {ex.mean():.2f} MeV ; ecart RMS "
              f"{np.sqrt(np.mean((moy - ex) ** 2)):.2f} MeV")
    print("""
  Le gap BCS suit les gaps mesures sans aucun ajustement, y compris sa
  chute vers les couches fermees. Mais a N = 82 et N = 126, il tombe
  exactement a ZERO (equation du gap sans solution) : c'est le defaut de
  BCS deja identifie dans la lecon, que corrige Lipkin-Nogami.
""")
    return occ, chaines


# ==========================================================================
def etape7_samarium(mm, be2):
    titre("7. LA TRANSITION DE FORME DES SAMARIUMS (N = 82 -> 92)")
    courbes = {}
    i0 = list(GRILLE_DELTA).index(0.0)
    print(f"\n  {'noyau':<8} {'N':>3} {'delta*':>7} {'beta2_eq':>9} {'gain (MeV)':>11} "
          f"{'beta2 mes.':>11} {'E(2+) keV':>10}")
    for n in range(82, 93, 2):
        tot, _ = mm.energie_forme(62, n)
        d, E = mm.forme(62, n)
        courbes[n] = tot - tot[i0]
        m = be2[(62, n)]
        print(f"  Sm-{62 + n:<4} {n:>3} {d:>7.3f} {float(beta2_equivalent(d)):>9.3f} "
              f"{E - tot[i0]:>11.2f} {m['beta2']:>11.3f} {m['E2_keV']:>10.1f}")
    c = courbes[88]
    i = i0
    while i > 0 and c[i - 1] < 0.4:
        i -= 1
    j = i0
    while j < len(c) - 1 and c[j + 1] < 0.4:
        j += 1
    print(f"""
  Le modele reproduit la TRANSITION : spherique jusqu'a N = 88, deforme des
  N = 90, comme les mesures (E(2+) divise par 2,7 entre Sm-150 et Sm-152).
  Sm-150 est un cas limite : sa courbe E(delta) reste a moins de 0,4 MeV de
  la sphere de delta = {GRILLE_DELTA[i]:+.3f} a {GRILLE_DELTA[j]:+.3f} (noyau "mou"). Les beta2
  mesures des isotopes spheriques (0,08 a 0,19) ne sont pas des
  deformations statiques : ce sont des amplitudes de VIBRATION autour de la
  sphere, que ce modele statique ne decrit pas.
""")
    return courbes


def etape8_quadrupole(Z, N, mm, be2):
    titre("8. MOMENTS QUADRUPOLAIRES : LE MODELE FACE AUX B(E2) MESURES")
    mesures = set(zip(Z.astype(int), N.astype(int)))
    pts = []
    for (z, n), m in be2.items():
        if (z, n) in mesures and z + n >= 40 and z % 2 == 0 and n % 2 == 0:
            d = mm.forme(z, n)[0]
            pts.append((z, n, m["beta2"], abs(float(beta2_equivalent(d))), d))
    P = np.array(pts)
    A = P[:, 0] + P[:, 1]
    print(f"""
  On compare des quantites MESURABLES : le moment quadrupolaire Q0 de la
  forme du modele (ellipsoide uniformement charge, exact) et celui deduit
  du B(E2) mesure, B(E2)up = (5/16 pi) e^2 Q0^2. Exprime en beta2 avec la
  convention des tables (beta2_equivalent), le rayon s'elimine.

  Conversion au premier ordre (lecon, version precedente) contre exacte :""")
    for d in (0.1, 0.2, 0.25, -0.1, -0.2):
        print(f"    delta = {d:+.2f} : premier ordre {d / 0.946:+.3f}, exact "
              f"{float(beta2_equivalent(d)):+.3f}")
    stats = {}
    for lab, s in (("A >= 40", A >= 40), ("A >= 150", A >= 150)):
        x, y, d = P[s, 2], P[s, 3], P[s, 4]
        g = x >= 0.2
        p = x < 0.12
        stats[lab] = dict(n=int(s.sum()), r=np.corrcoef(x, y)[0, 1], n_def=int(g.sum()),
                          ok_def=int((np.abs(d[g]) > 0.05).sum()),
                          rapport=float(np.median(y[g] / x[g])), n_sph=int(p.sum()),
                          ok_sph=int((np.abs(d[p]) <= 0.05).sum()))
        st = stats[lab]
        print(f"""
  {lab} : {st['n']} noyaux pair-pair (masse et B(E2) mesures)
    correlation beta2 modele / mesure : {st['r']:.2f}
    beta2 mesure >= 0,2 : {st['n_def']} noyaux, dont {st['ok_def']} predits deformes ;
      rapport median modele / mesure : {st['rapport']:.2f}
    beta2 mesure < 0,12 : {st['n_sph']} noyaux, dont {st['ok_sph']} predits spheriques""")
    ecarts = [(z, n, m["beta2_table"], m["beta2"], m["reference"]) for (z, n), m in
              sorted(be2.items()) if abs(m["beta2_table"] / m["beta2"] - 1) > 0.02]
    print(f"""
  beta2 recalcule a partir de B(E2) pour tous les noyaux. Ecarts > 2 % avec
  la valeur imprimee dans les tables : {len(ecarts)} noyaux sur {len(be2)}, dont :""")
    for z, n, bt, b, ref in ecarts:
        if z + n >= 150:
            print(f"    Z={z} N={n} : imprime {bt:.4f}, recalcule {b:.4f}  ({ref})")
    i = np.where((P[:, 0] == 92) & (P[:, 1] == 146))[0][0]
    z, n = 92, 146
    R0 = 1.2 * (z + n) ** (1 / 3) / 10                     # en sqrt(barn)
    Q0_mes = np.sqrt(16 * np.pi / 5 * be2[(z, n)]["BE2"])
    q = rapport_axes(P[i, 4])
    Q0_mod = 0.4 * z * R0 ** 2 * (q ** (4 / 3) - q ** (-2 / 3))
    print(f"""
  U-238 : delta* = {P[i, 4]:.3f}. Premier ordre : beta2 = {P[i, 4] / 0.946:.3f} ; exact :
  beta2_eq = {P[i, 3]:.3f}. Mesure ({be2[(z, n)]['reference']}) : {P[i, 2]:.3f}.
  Q0 : modele {Q0_mod:.2f} b, mesure {Q0_mes:.2f} b (B(E2) = {be2[(z, n)]['BE2']:.2f} e^2 b^2).
  La conversion au premier ordre expliquait {(P[i, 3] - P[i, 4] / 0.946) / (P[i, 2] - P[i, 4] / 0.946):.0%} de l'ecart.
  Il reste une sous-estimation de ~{100 * (1 - P[i, 3] / P[i, 2]):.0f} % : forme limitee aux
  ellipsoides (pas d'eps4), potentiel de Nilsson etendu a grand delta.
""")
    rates = (A < 150) & (P[:, 2] >= 0.2) & (np.abs(P[:, 4]) <= 0.05)
    print(f"  40 <= A < 150 : {rates.sum()} noyaux mesures deformes (beta2 >= 0,2) predits "
          "spheriques, par region :")
    for lab, s_ in (("A < 70", A < 70), ("70 <= A < 90", (A >= 70) & (A < 90)),
                    ("90 <= A < 120", (A >= 90) & (A < 120)), ("120 <= A < 150", A >= 120)):
        sel = rates & s_
        noms = ", ".join(f"{int(z)}-{int(z + n)}" for z, n in P[sel, :2])
        print(f"    {lab:<15} {sel.sum():>3}  (Z-A : {noms})")
    return P, stats


def etape9_inertie(be2, mm):
    titre("9. LE MOMENT D'INERTIE DE L'U-238 : UNE SIGNATURE DE L'APPARIEMENT")
    z, n = 92, 146
    A = z + n
    hb2_m = 197.327 ** 2 / 938.92                       # hbar^2/m, MeV fm^2
    R = 1.2 * A ** (1 / 3)
    q = rapport_axes(mm.forme(z, n)[0])
    a2b2 = R ** 2 * (q ** (4 / 3) + q ** (-2 / 3))       # a^2 + b^2, fm^2
    J_rig = A * a2b2 / 5                                 # en unites de m (fm^2)
    E2_rig = 6 * hb2_m / (2 * J_rig) * 1000              # keV
    E2 = be2[(z, n)]["E2_keV"]
    print(f"""
  Bande de rotation : E(I) = (hbar^2 / 2J) I(I+1), donc E(2+) = 6 hbar^2 / 2J.
  Rotateur rigide (ellipsoide du modele, J = (m A / 5)(a^2 + b^2)) :
      E(2+) = {E2_rig:.1f} keV
  Mesure : E(2+) = {E2:.1f} keV, soit J / J_rigide = {E2_rig / E2:.2f}.
  Le noyau tourne avec un moment d'inertie moitie moindre que celui d'un
  solide : les paires de nucleons, comme un superfluide, ne suivent pas
  toutes la rotation (Bohr et Mottelson, vol. II). Un modele sans
  appariement predirait le moment d'inertie rigide.
""")
    return E2_rig, E2


# ==========================================================================
# Figures
# ==========================================================================
def figures(niveaux, strut, Lsph, bilan, occ, chaines, courbes_sm, formes_sm, be2, P):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def figure(h=5.6, gauche=0.10, droite=0.80, bas=0.14, haut=0.83, w=9):
        fig, ax = plt.subplots(figsize=(w, h), facecolor=SURFACE)
        fig.subplots_adjust(left=gauche, right=droite, top=haut, bottom=bas)
        return fig, ax

    def note(fig, texte):
        fig.text(0.02, 0.015, texte, fontsize=8.5, color=DISCRET)

    def sauver(fig, nom):
        fig.savefig(nom, dpi=130, facecolor=SURFACE)
        plt.close(fig)
        print(f"  --> {nom}")

    # ---------------------------------------------------------------- 1
    fig, ax = figure(h=8.4, gauche=0.08, droite=0.86, bas=0.07, haut=0.85)
    style(ax)
    ax.grid(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    X = {"osc": 0.0, "l": 1.15, "j": 2.3}
    L = 0.55
    ymin, ymax = 1.2, 8.15
    niv = [v for v in niveaux if v["N"] <= 6]
    intrus = lambda v: v["N"] >= 3 and v["l"] == v["N"] and v["jx2"] == 2 * v["l"] + 1
    for N in range(7):
        ax.plot([X["osc"], X["osc"] + L], [N + 1.5] * 2, color=ENCRE, lw=2.2,
                solid_capstyle="butt")
        ax.text(X["osc"] - 0.05, N + 1.5, f"N = {N}", ha="right", va="center",
                fontsize=10, color=ENCRE_2)
    vus = set()
    for v in niv:
        col = SERIES[1] if intrus(v) else ENCRE
        if (v["N"], v["l"]) not in vus:
            vus.add((v["N"], v["l"]))
            ax.plot([X["osc"] + L, X["l"]], [v["e_osc"], v["e_l"]], color=AXE, lw=0.7)
            ax.plot([X["l"], X["l"] + L], [v["e_l"]] * 2, color=ENCRE, lw=2.0,
                    solid_capstyle="butt")
        ax.plot([X["l"] + L, X["j"]], [v["e_l"], v["e"]], color=AXE, lw=0.7)
        ax.plot([X["j"], X["j"] + L], [v["e"]] * 2, color=col, lw=2.0,
                solid_capstyle="butt", zorder=3 if intrus(v) else 2)
    # etiquettes j a droite, ecartees
    ys = ecarter([v["e"] for v in niv], 0.135)
    for v, y in zip(niv, ys):
        ax.plot([X["j"] + L, X["j"] + L + 0.08], [v["e"], y], color=AXE, lw=0.6)
        ax.text(X["j"] + L + 0.1, y, nom_orbitale(v["N"], v["l"], v["jx2"]),
                va="center", fontsize=8.6, color=ENCRE,
                fontweight="semibold" if intrus(v) else None)
        ax.text(X["j"] + L + 0.62, y, f"{v['deg']}", va="center", ha="right",
                fontsize=8.0, color=DISCRET)
    # nombres magiques aux grands ecarts
    for v in niv:
        if v["ecart"] > 0.25 and v["cumul"] <= 126:
            y = v["e"] + v["ecart"] / 2
            magique = v["cumul"] in (2, 8, 20, 28, 50, 82, 126)
            ax.text(X["j"] + L + 0.95, y, str(v["cumul"]), va="center", ha="left",
                    fontsize=11.5 if magique else 9.5,
                    color=ENCRE if magique else DISCRET,
                    fontweight="bold" if magique else None)
    # couches de l'oscillateur
    cum = 0
    for N in range(7):
        cum += (N + 1) * (N + 2)
        ax.text(X["osc"] + L / 2, N + 1.5 + 0.5, str(cum), ha="center", va="center",
                fontsize=9.5, color=DISCRET)
    for k, (cle, lab) in enumerate((("osc", "oscillateur\nharmonique"),
                                    ("l", "+ μ (l² − ⟨l²⟩)"),
                                    ("j", "+ spin-orbite\n2 l·s"))):
        ax.text(X[cle] + L / 2, ymax + 0.05, lab, ha="center", va="bottom", fontsize=10.5,
                color=ENCRE_2)
    ax.text(X["j"] + L + 0.62, ymax + 0.05, "places", ha="right", va="bottom",
            fontsize=8.5, color=DISCRET)
    ax.set_xlim(-0.55, 4.0)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks([])
    ax.set_ylabel("énergie (ħω₀)")
    fig.text(0.08, 0.968, "D'où viennent les nombres magiques", fontsize=13.5,
             color=ENCRE, fontweight="semibold")
    fig.text(0.08, 0.944, "Neutrons, κ et μ de Bengtsson et Ragnarsson. En orange : "
             "l'orbitale de j maximal de chaque couche,", fontsize=10, color=ENCRE_2)
    fig.text(0.08, 0.926, "que le spin-orbite fait descendre dans la couche inférieure.",
             fontsize=10, color=ENCRE_2)
    note(fig, "Calcul : etude_deformation.py. Nombres gris dans la 1re colonne : "
         "couches fermées de l'oscillateur pur.")
    sauver(fig, "deformation_niveaux.png")

    # ---------------------------------------------------------------- 2
    fig, ax = figure(h=7.6, gauche=0.09, droite=0.85, bas=0.09, haut=0.86)
    style(ax)
    dmax = 0.32
    ds = np.round(np.linspace(0.0, dmax, 65), 5)
    lignes = lignes_nilsson(ds, "neutron")
    emin, emax = 7.2, 7.95
    signe = {"+": "⁺", "−": "⁻"}
    fin = []
    for Lg in lignes:
        E = Lg["E"]
        if E.max() < emin - 0.3 or E.min() > emax + 0.3:
            continue
        col = SERIES[0] if Lg["parite"] == "+" else SERIES[1]
        ax.plot(ds, E, color=col, lw=1.5, solid_joinstyle="round")
        if emin < E[-1] < emax:
            fin.append((E[-1], f"{Lg['om_x2']}/2{signe[Lg['parite']]}", col))
    # niveau de Fermi de N = 146 et gap deforme N = 152
    tous = np.array([spectre(d, "neutron") for d in ds])
    f = 0.5 * (tous[:, 72] + tous[:, 73])
    ax.plot(ds, f, color=ENCRE, lw=1.6, ls="--")
    k = int(np.argmin(np.abs(ds - 0.10)))
    ax.annotate("niveau de Fermi\npour N = 146 (U-238)", (ds[k], f[k]), xytext=(0.035, 7.27),
                fontsize=9.5, color=ENCRE,
                arrowprops=dict(arrowstyle="-", color=ENCRE_2, lw=0.8))
    k25 = int(np.argmin(np.abs(ds - 0.25)))
    bas_, haut_ = tous[k25, 75], tous[k25, 76]
    ax.annotate("", (0.25, haut_), xytext=(0.25, bas_),
                arrowprops=dict(arrowstyle="<->", color=ENCRE, lw=1.2,
                                shrinkA=0, shrinkB=0))
    ax.text(0.254, 0.5 * (bas_ + haut_), "gap\nN = 152", fontsize=10, color=ENCRE,
            fontweight="semibold", va="center",
            bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.9))
    ax.axvline(0.193, color=DISCRET, lw=0.9, ls=":")
    ax.text(0.191, emax - 0.012, "δ* de l'U-238", fontsize=9, color=ENCRE_2, va="top",
            ha="right")
    # orbitales spheriques a gauche, dans le cadre
    e0 = sorted({(round(Lg["E"][0], 4), Lg["parent"]) for Lg in lignes
                 if emin < Lg["E"][0] < emax})
    groupes = []                                   # orbitales quasi degenerees
    for y, nom in e0:
        if groupes and y - groupes[-1][0] < 0.03:
            if nom not in groupes[-1][1]:
                groupes[-1][1].append(nom)
        else:
            groupes.append((y, [nom]))
    for y, noms in groupes:
        ax.text(0.004, y + 0.016, " et ".join(noms), ha="left", va="bottom",
                fontsize=9.5, color=ENCRE, fontweight="semibold",
                bbox=dict(boxstyle="round,pad=0.12", fc=SURFACE, ec="none", alpha=0.85))
    # Omega^pi a droite
    fin.sort()
    ys = ecarter([f_[0] for f_ in fin], 0.026)
    tr = ax.get_yaxis_transform()
    for (y0, lab, col), y1 in zip(fin, ys):
        ax.plot([1.0, 1.025], [y0, y1], color=DISCRET, lw=0.5, clip_on=False, transform=tr)
        ax.plot([1.025, 1.06], [y1, y1], color=col, lw=2.2, clip_on=False, transform=tr)
        ax.text(1.07, y1, lab, va="center", fontsize=9, color=ENCRE, transform=tr)
    ax.set_xlim(0, dmax)
    ax.set_ylim(emin, emax)
    ax.set_xlabel("déformation δ (allongée)")
    ax.set_ylabel("énergie du niveau (ħω₀)")
    fig.text(0.09, 0.955, "Diagramme de Nilsson des actinides (neutrons)", fontsize=13.5,
             color=ENCRE, fontweight="semibold")
    fig.text(0.09, 0.93, "À droite, Ω et la parité de chaque niveau (bleu : +, orange : −). "
             "À gauche, l'orbitale sphérique d'origine.", fontsize=10, color=ENCRE_2)
    fig.text(0.09, 0.908, "Ω petit : l'énergie baisse quand le noyau s'allonge. "
             "Ω grand : elle monte.", fontsize=10, color=ENCRE_2)
    note(fig, "Calcul : nilsson.py, diagonalisation par bloc (Ω, parité). "
         "Deux niveaux d'un même bloc ne se croisent jamais.")
    sauver(fig, "deformation_nilsson.png")

    # ---------------------------------------------------------------- 3
    fig, ax = figure(droite=0.95)
    style(ax)
    e = np.linspace(1.0, 8.4, 3000)
    escalier = np.array([2 * np.sum(Lsph.niveaux <= x) for x in e])
    lisse = np.interp(e, Lsph.e, Lsph.Ncum)
    ax.plot(e, escalier, color=SERIES[0], lw=2, drawstyle="steps-post")
    ax.plot(e, lisse, color=SERIES[1], lw=2)
    for n in (8, 20, 28, 50, 82, 126):
        k = np.where(escalier == n)[0]
        if len(k):
            xm = 0.5 * (e[k[0]] + e[k[-1]])
            ax.text(xm, n + 6, str(n), ha="center", fontsize=10, color=ENCRE,
                    fontweight="semibold")
    ax.set_xlim(1.0, 8.4)
    ax.set_ylim(0, 200)
    ax.set_xlabel("énergie e (ħω₀)")
    ax.set_ylabel("nombre de neutrons sous e")
    ax.legend(["niveaux réels (escalier)", "lissage de Strutinsky"], loc="upper left",
              fontsize=10, frameon=False, labelcolor=ENCRE)
    titre_axe(ax, "L'escalier des niveaux et sa version lisse",
              "Neutrons, δ = 0. Les marches longues sont les couches fermées ; "
              "γ = 1,2 ħω₀.")
    note(fig, "δE_couches = (somme des niveaux remplis) − (même somme sur la courbe lisse). "
         "Calcul : strutinsky.py.")
    sauver(fig, "deformation_escalier.png")

    # ---------------------------------------------------------------- 4
    fig, ax = figure(droite=0.76)
    style(ax)
    for k, (d, lab) in enumerate(((0.0, "sphère, δ = 0"), (0.25, "allongé, δ = 0,25"))):
        n, dE = strut[d]
        ax.plot(n, dE, color=SERIES[k], lw=2, marker="o", ms=3)
    for m in (28, 50, 82, 126):
        ax.axvline(m, color=AXE, lw=0.9, ls=":")
        ax.text(m, 1.72, str(m), ha="center", fontsize=9.5, color=ENCRE_2)
    ax.axhline(0, color=AXE, lw=0.8)
    ax.set_xlim(8, 184)
    ax.set_ylim(-1.3, 1.85)
    ax.set_xlabel("nombre de neutrons n")
    ax.set_ylabel("correction de couches (ħω₀)")
    etiquettes(ax, 184, [strut[0.0][1][-1], strut[0.25][1][-1]],
               ["sphère, δ = 0", "allongé, δ = 0,25"], SERIES[:2], ecart=0.25)
    titre_axe(ax, "Où la sphère coûte, la déformation rapporte",
              "Neutrons. Pour A ≈ 200, 1 ħω₀ ≈ 7 MeV.")
    note(fig, "Calcul : etude_deformation.py (Nilsson + Strutinsky, γ = 1,2 ħω₀).")
    sauver(fig, "deformation_couches_n.png")

    # ---------------------------------------------------------------- 5
    fig, ax = figure(droite=0.80)
    style(ax)
    g = GRILLE_DELTA
    series = [(bilan["mac"], "goutte"), (bilan["couches"], "couches"),
              (bilan["appar"], "appariement"), (bilan["tot"], "total")]
    cols = [SERIES[0], SERIES[1], SERIES[2], ENCRE]
    for (y, lab), c in zip(series, cols):
        ax.plot(g, y, color=c, lw=2.6 if lab == "total" else 2)
    k = int(np.argmin(bilan["tot"]))
    ax.plot(g[k], bilan["tot"][k], "o", ms=8, color=ENCRE, mec=SURFACE, mew=2, zorder=5)
    ax.axhline(0, color=AXE, lw=0.8)
    ax.axvline(0, color=AXE, lw=0.8)
    ax.set_xlim(g[0], g[-1])
    ax.set_xlabel("déformation δ  (< 0 aplati, > 0 allongé)")
    ax.set_ylabel("contribution à l'énergie (MeV)")
    etiquettes(ax, g[-1], [s[0][-1] for s in series], [s[1] for s in series], cols, ecart=1.6)
    titre_axe(ax, "Erbium 166 : qui veut la déformation, qui s'y oppose",
              "La sphère est un sommet pour les couches, un creux pour l'appariement. "
              "Point : minimum du total.")
    note(fig, "Couches et appariement : protons + neutrons. Calcul : etude_deformation.py.")
    sauver(fig, "deformation_bilan.png")

    # ---------------------------------------------------------------- 6
    fig, ax = figure(droite=0.80)
    style(ax)
    eps, v2, D = occ["eps"] - occ["lam"], occ["v2"], occ["Delta"]
    lam = 0.0
    x = np.linspace(eps.min() - 0.5, eps.max() + 0.5, 600)
    ax.plot(x, (x < lam).astype(float), color=DISCRET, lw=1.4, drawstyle="steps-post")
    ax.plot(x, 0.5 * (1 - (x - lam) / np.sqrt((x - lam) ** 2 + D * D)), color=SERIES[0], lw=2)
    ax.plot(eps, v2, "o", ms=7, color=SERIES[0], mec=SURFACE, mew=1.5, zorder=5)
    ax.axvspan(lam - D, lam + D, color=SERIES[0], alpha=0.08, lw=0)
    ax.axvline(lam, color=ENCRE_2, lw=0.8, ls=":")
    ax.text(lam, 1.07, "λ", ha="center", fontsize=11, color=ENCRE)
    ax.text(lam + D, 1.07, "λ + Δ", ha="center", fontsize=9.5, color=ENCRE_2)
    ax.text(lam - D, 1.07, "λ − Δ", ha="center", fontsize=9.5, color=ENCRE_2)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlabel("énergie du niveau par rapport au potentiel chimique, ε − λ (MeV)")
    ax.set_ylabel("probabilité d'occupation v²")
    etiquettes(ax, x[-1], [0.18, 0.0], ["BCS", "sans\nappariement"], [SERIES[0], DISCRET],
               ecart=0.22)
    titre_axe(ax, f"Étain 116 : la surface de Fermi devient floue",
              f"Neutrons. Gap BCS Δ = {virgule(D)} MeV, sans paramètre ajusté.")
    note(fig, "Chaque point est un niveau de Nilsson (δ = 0, plusieurs points superposés "
         "par couche j). Calcul : strutinsky.py.")
    sauver(fig, "deformation_occupations.png")

    # ---------------------------------------------------------------- 7
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 5.4), facecolor=SURFACE, sharey=True)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.80, bottom=0.14, wspace=0.08)
    for ax, (z, magique) in zip(axs, ((50, 82), (82, 126))):
        style(ax)
        m, e = chaines[z]
        ax.plot(m[:, 0], m[:, 1], "-o", color=SERIES[0], lw=2, ms=5, mec=SURFACE, mew=1)
        ax.plot(e[:, 0], e[:, 1], "s", color=SERIES[1], ms=6, mec=SURFACE, mew=1)
        nn = np.linspace(m[:, 0].min(), m[:, 0].max(), 100)
        ax.plot(nn, 12 / np.sqrt(z + nn), color=DISCRET, lw=1.3, ls="--")
        ax.axvline(magique, color=AXE, lw=0.9, ls=":")
        ax.text(magique - 0.6, 2.3, f"N = {magique}", ha="right", fontsize=9.5,
                color=ENCRE_2)
        ax.set_xlabel("nombre de neutrons N")
        ax.text(0.0, 1.03, f"{SYMBOLES[z]} (Z = {z})", transform=ax.transAxes, fontsize=11.5,
                color=ENCRE, fontweight="semibold")
        ax.set_ylim(-0.08, 2.45)
    axs[0].set_ylabel("gap d'appariement des neutrons (MeV)")
    axs[0].legend(["Δ BCS du modèle (N pair)", "Δ⁽³⁾ mesuré (N impair, AME2020)",
                   "gap moyen 12/√A"], loc="upper left", fontsize=9.5, frameon=False,
                  labelcolor=ENCRE)
    fig.text(0.08, 0.93, "Le gap BCS face aux masses mesurées", fontsize=13.5, color=ENCRE,
             fontweight="semibold")
    fig.text(0.08, 0.895, "Bon accord en milieu de couche ; aux couches fermées, BCS "
             "tombe exactement à zéro.", fontsize=10.5, color=ENCRE_2)
    note(fig, "Δ⁽³⁾(N) = −(−1)^N [B(N+1) − 2B(N) + B(N−1)] / 2. "
         "Calcul : etude_deformation.py.")
    sauver(fig, "deformation_gaps.png")

    # ---------------------------------------------------------------- 8
    fig, ax = figure(droite=0.95)
    style(ax)
    rampe = BLEUS[1:]
    ns = sorted(courbes_sm)
    for k, n in enumerate(ns):
        y = courbes_sm[n]
        ax.plot(GRILLE_DELTA, y, color=rampe[k], lw=2)
        i = int(np.argmin(y))
        ax.plot(GRILLE_DELTA[i], y[i], "o", ms=6, color=rampe[k], mec=SURFACE, mew=1.5)
    ax.axhline(0, color=AXE, lw=0.8)
    ax.axvline(0, color=AXE, lw=0.8)
    ax.set_xlim(-0.3, 0.45)
    ax.set_ylim(-4, 8)
    ax.set_xlabel("déformation δ")
    ax.set_ylabel("E(δ) − E(0)  (MeV)")
    from matplotlib.lines import Line2D
    ax.legend([Line2D([], [], color=rampe[k], lw=2.4) for k in range(len(ns))],
              [f"Sm-{62 + n}  (N = {n})" for n in ns], loc="lower left", fontsize=9.5,
              frameon=True, facecolor=SURFACE, edgecolor="none", framealpha=0.95,
              labelcolor=ENCRE)
    titre_axe(ax, "Samarium : la sphère cède à N = 90",
              "Énergie de déformation du modèle, du Sm-144 (magique) au Sm-154. "
              "Point : minimum.")
    note(fig, "Calcul : etude_deformation.py (modèle mac-mic de la leçon, 8 coefficients "
         "ajustés sur AME2020).")
    sauver(fig, "deformation_samarium.png")

    # ---------------------------------------------------------------- 9
    fig, axs = plt.subplots(2, 1, figsize=(9, 7.4), facecolor=SURFACE, sharex=True)
    fig.subplots_adjust(left=0.11, right=0.74, top=0.87, bottom=0.14, hspace=0.18)
    n = np.array(ns)
    bm = np.array([be2[(62, k)]["beta2"] for k in ns])
    e2 = np.array([be2[(62, k)]["E2_keV"] for k in ns])
    bmod = np.array([abs(float(beta2_equivalent(formes_sm[k]))) for k in ns])
    ax = axs[0]
    style(ax)
    ax.plot(n, bm, "-s", color=SERIES[1], lw=2, ms=7, mec=SURFACE, mew=1.5)
    ax.plot(n, bmod, "-o", color=SERIES[0], lw=2, ms=7, mec=SURFACE, mew=1.5)
    ax.set_ylabel("β₂")
    ax.set_ylim(-0.02, 0.4)
    etiquettes(ax, n[-1], [bm[-1], bmod[-1]], ["mesuré (B(E2))", "modèle (δ*, via Q₀)"],
               [SERIES[1], SERIES[0]], ecart=0.06)
    ax.text(0, 1.04, "Déformation quadrupolaire", transform=ax.transAxes, fontsize=11.5,
            color=ENCRE, fontweight="semibold")
    ax = axs[1]
    style(ax)
    ax.plot(n, e2, "-s", color=SERIES[1], lw=2, ms=7, mec=SURFACE, mew=1.5)
    ax.set_yscale("log")
    ax.set_ylim(50, 3000)
    from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
    ax.yaxis.set_major_locator(FixedLocator([50, 100, 200, 500, 1000, 2000]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    for k, v in zip(n, e2):
        ax.annotate(f"{v:.0f}", (k, v), xytext=(0, 9), textcoords="offset points",
                    ha="center", fontsize=9, color=ENCRE_2)
    ax.set_ylabel("E(2⁺) (keV)")
    ax.set_xlabel("nombre de neutrons N")
    ax.set_xticks(n)
    ax.set_xticklabels([f"{k}\nSm-{62 + k}" for k in n])
    ax.text(0, 1.04, "Énergie du premier état 2⁺ (mesurée)", transform=ax.transAxes,
            fontsize=11.5, color=ENCRE, fontweight="semibold")
    fig.text(0.11, 0.955, "Samarium : la transition vue par les mesures", fontsize=13.5,
             color=ENCRE, fontweight="semibold")
    fig.text(0.11, 0.925, "β₂ saute et E(2⁺) s'effondre entre N = 88 et 90 : le noyau "
             "se met à tourner comme un rotateur.", fontsize=10.5, color=ENCRE_2)
    note(fig, "Mesures : Pritychenko et al., ADNDT 107 (2016), via NNDC. "
         "Calcul : etude_deformation.py.")
    sauver(fig, "deformation_samarium_mesures.png")

    # ---------------------------------------------------------------- 10
    fig, ax = figure(h=6.4, gauche=0.11, droite=0.95, bas=0.12, haut=0.85)
    style(ax)
    A = P[:, 0] + P[:, 1]
    for k, (s, lab) in enumerate(((A < 150, "40 ≤ A < 150"), (A >= 150, "A ≥ 150"))):
        ax.plot(P[s, 2], P[s, 3], "o", ms=5.5, color=SERIES[k], mec=SURFACE, mew=0.8,
                alpha=0.9, label=lab)
    ax.plot([0, 0.5], [0, 0.5], color=DISCRET, lw=1, ls="--")
    for z, n_, nom, dxy in ((92, 146, "U-238", (8, -12)), (68, 98, "Er-166", (-8, 10)),
                            (62, 90, "Sm-152", (8, -4)), (82, 126, "Pb-208", (8, 6)),
                            (62, 88, "Sm-150", (6, 8))):
        i = np.where((P[:, 0] == z) & (P[:, 1] == n_))[0]
        if len(i):
            i = i[0]
            ax.annotate(nom, (P[i, 2], P[i, 3]), xytext=dxy, textcoords="offset points",
                        fontsize=9.5, color=ENCRE, ha="left" if dxy[0] > 0 else "right")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(-0.01, 0.42)
    p1, p2 = ax.transData.transform([(0.3, 0.3), (0.4, 0.4)])
    angle = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
    ax.text(0.385, 0.392, "modèle = mesure", rotation=angle, rotation_mode="anchor",
            fontsize=9, color=DISCRET, ha="center")
    ax.set_xlabel("β₂ mesuré (B(E2), tables NNDC)")
    ax.set_ylabel("β₂ du modèle (forme δ*, via Q₀)")
    ax.legend(loc="upper left", fontsize=10, frameon=False, labelcolor=ENCRE)
    titre_axe(ax, "Le modèle face aux moments quadrupolaires mesurés",
              "Noyaux pair-pair, A ≥ 40, masse et B(E2) mesurées.")
    note(fig, "Les noyaux sphériques du modèle (axe horizontal) ont un β₂ mesuré non nul : "
         "vibrations, pas déformation statique.")
    sauver(fig, "deformation_quadrupole.png")


def main():
    niveaux = etape1_niveaux()
    etape2_pentes()
    Z, N, B, mm, mm0 = ajuster()
    etape3_nilsson({(int(z), int(n)): b for z, n, b in zip(Z, N, B)})
    strut, Lsph = etape4_strutinsky()
    bilan = etape5_bilan(mm, mm0)
    occ, chaines = etape6_bcs(Z, N, B, mm)
    be2 = donnees_be2()
    courbes_sm = etape7_samarium(mm, be2)
    formes_sm = {n: mm.forme(62, n)[0] for n in courbes_sm}
    P, _ = etape8_quadrupole(Z, N, mm, be2)
    etape9_inertie(be2, mm)
    titre("FIGURES")
    figures(niveaux, strut, Lsph, bilan, occ, chaines, courbes_sm, formes_sm, be2, P)


if __name__ == "__main__":
    main()
