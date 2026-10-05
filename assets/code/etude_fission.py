"""
etude_fission.py -- De l'uranium au reacteur.

Sept questions, dans l'ordre logique :

  1. Pourquoi un noyau lourd libere-t-il de l'energie en se cassant ?
  2. Pourquoi l'U235 et pas l'U238 ?
  3. Combien d'energie, et sous quelle forme ?
  4. Pourquoi faut-il ralentir les neutrons, et avec quoi ?
  5. Le bilan complet : la formule des quatre facteurs.
  6. Pourquoi un REP ne peut pas marcher a l'uranium naturel.
  7. Du microscopique aux constantes a deux groupes (bouclage avec coeur.py).

Lancer :  python etude_fission.py
"""

import numpy as np
from fission import (energie_liaison, liaison_par_nucleon, energie_fission,
                     energie_separation_neutron, BILAN_ENERGIE, NUCLIDES,
                     BARN, densite_atomique, xi_moyen, nombre_de_chocs,
                     pouvoir_ralentissement, rapport_moderation,
                     facteur_eta, facteur_f, facteur_p_homogene,
                     facteur_p_heterogene, k_infini_quatre_facteurs)

# Densites de reference
N_U_UO2 = densite_atomique(10.4, 270.03)     # noyaux d'U par cm3 d'UO2
N_H2O = densite_atomique(0.72, 18.015)       # molecules d'eau par cm3 (a 310 C)
N_D2O = densite_atomique(1.10, 20.03)
RAYON_PASTILLE = 0.41                        # cm


def titre(n, texte):
    print()
    print("=" * 74)
    print(f"{n}. {texte}")
    print("=" * 74)


# ==========================================================================
def etape1_pourquoi_la_fission_libere():
    titre(1, "POURQUOI LA FISSION LIBERE-T-ELLE DE L'ENERGIE ?")
    print("""
  Tout tient dans une seule courbe : l'energie de liaison par nucleon B/A.
  Plus B/A est grand, plus le noyau est solidement lie, donc stable. La
  courbe monte de l'hydrogene au fer, puis REDESCEND.

  Un noyau lourd qui se casse en deux morceaux plus proches du fer devient
  donc plus lie. La difference d'energie de liaison est liberee.

  On la calcule avec la formule semi-empirique de masse (Bethe-Weizsacker),
  qui traite le noyau comme une goutte de liquide chargee.
""")
    print(f"  {'noyau':>8} {'A':>5} {'Z':>4} {'B/A calcule':>13} {'B/A mesure':>12}")
    reference = {"He4": 7.074, "Fe56": 8.790, "Kr92": 8.513,
                 "Ba141": 8.326, "U235": 7.591, "U238": 7.570}
    for nom, A, Z in (("He4", 4, 2), ("Fe56", 56, 26), ("Kr92", 92, 36),
                      ("Ba141", 141, 56), ("U235", 235, 92), ("U238", 238, 92)):
        print(f"  {nom:>8} {A:>5} {Z:>4} {liaison_par_nucleon(A, Z):>13.3f} "
              f"{reference[nom]:>12.3f}")
    print("""
  Excellent pour les noyaux lourds (U235 : 7.600 calcule contre 7.591
  mesure). Mauvais pour l'helium 4 (5.49 contre 7.07) : la goutte liquide
  ignore les couches nucleaires, et He4 est doublement magique, donc
  anormalement lie. Retenir la limite du modele.

  LE TERME COUPABLE. Dans B = a_V*A - a_S*A^(2/3) - a_C*Z^2/A^(1/3) - ...,
  c'est le terme coulombien en Z^2/A^(1/3) qui fait redescendre la courbe :
  les protons se repoussent tous mutuellement (en Z^2), alors que la force
  nucleaire ne lie qu'aux voisins immediats (en A). Passe le fer, la
  repulsion l'emporte. L'uranium est deja presque instable ; il suffit de
  le deformer un peu.
""")


# ==========================================================================
def etape2_pourquoi_u235():
    titre(2, "POURQUOI L'U235 ET PAS L'U238 ?")
    print("""
  Les deux sont fissionnables en principe. Un seul est FISSILE, c'est-a-dire
  cassable par un neutron LENT. La difference est un des plus beaux
  resultats de la physique nucleaire, et elle tient a l'appariement des
  nucleons.

  Quand un noyau capture un neutron, il recoit l'energie de liaison de ce
  neutron, S_n, sous forme d'excitation. S'il en recoit plus que la
  barriere de fission, il casse. Sinon, il se contente d'emettre un gamma.
""")
    barriere = {"U236": 6.2, "U239": 6.6}
    print(f"  {'reaction':<26} {'noyau':>7} {'S_n calcule':>12} "
          f"{'barriere':>10} {'fissile ?':>12}")
    for reaction, A, Z, comp in (("U235 + n -> U236*", 236, 92, "U236"),
                                 ("U238 + n -> U239*", 239, 92, "U239")):
        sn = float(energie_separation_neutron(A, Z))
        b = barriere[comp]
        verdict = "OUI" if sn > b else "non"
        print(f"  {reaction:<26} {comp:>7} {sn:>10.2f} MeV {b:>8.1f} MeV "
              f"{verdict:>12}")
    print("""
  Le mecanisme, en une phrase : l'U235 a un nombre IMPAIR de neutrons (143).
  Le neutron incident vient completer une paire, et l'energie d'appariement
  ainsi liberee (~0.6 MeV) suffit a franchir la barriere. L'U238 a un
  nombre PAIR de neutrons (146) ; le neutron incident reste celibataire, il
  ne recoit pas ce bonus, et l'excitation est insuffisante.

  Un neutron thermique (0.025 eV) apporte une energie cinetique
  parfaitement negligeable. Tout vient de S_n. C'est pourquoi la fissilite
  ne depend pas de la vitesse du neutron mais de la PARITE du noyau cible.

  Meme regle pour les autres fissiles : Pu239 (N=145, impair), U233
  (N=141, impair). Aucune exception.

  L'U238 fissionne quand meme, mais il faut lui apporter la difference en
  energie cinetique : environ 1 MeV de seuil. C'est le facteur epsilon de
  la formule des quatre facteurs, qui vaut ~1.03 : 3 % des fissions d'un
  REP ont lieu sur l'U238, par des neutrons encore rapides.
""")


# ==========================================================================
def etape3_bilan_energie():
    titre(3, "COMBIEN D'ENERGIE, ET SOUS QUELLE FORME ?")
    print("""
  On calcule Q = B(fragment1) + B(fragment2) - B(U236). Les neutrons libres
  ont une energie de liaison nulle et ne comptent pas.
""")
    print(f"  {'partition de U236':<32} {'Q (goutte liquide)':>20}")
    partitions = [((118, 46), (118, 46), 0, "symetrique"),
                  ((133, 52), (100, 40), 3, ""),
                  ((137, 54), (97, 38), 2, ""),
                  ((141, 56), (92, 36), 3, "la plus citee")]
    for (A1, Z1), (A2, Z2), nn, note in partitions:
        Q = float(energie_fission(236, 92, A1, Z1, A2, Z2))
        label = f"{A1}/{Z1} + {A2}/{Z2} + {nn}n"
        print(f"  {label:<32} {Q:>14.1f} MeV   {note}")
    print("""
  POURQUOI LA FISSION REELLE EST-ELLE ASYMETRIQUE ? Le calcul donne un
  maximum pour la fission SYMETRIQUE (185 MeV), alors que la fission reelle
  de l'U235 par neutron lent fait deux bosses vers A=95 et A=139.

  Attention au piege : ce n'est PAS une erreur du modele. Les masses
  MESUREES donnent elles aussi plus d'energie au partage symetrique
  (193 MeV contre 167 pour Ba141+Kr92) -- voir analyse_liaison.py, qui
  refait le calcul sur donnees experimentales.

  L'asymetrie n'est donc pas un effet de bilan energetique. Elle vient de
  la DYNAMIQUE du noyau au point de scission : la surface d'energie
  potentielle au point selle favorise des fragments proches des couches
  fermees Z=50 et N=82. Le noyau ne choisit pas le partage le plus
  exothermique, il suit le chemin le plus facile.

  Pour Ba141 + Kr92 + 3n, le calcul donne 156 MeV, la mesure 167 MeV.
  Erreur de 6 %, honorable pour un modele a cinq parametres. La vraie
  limite du modele est ailleurs : voir analyse_liaison.py, qui confronte la
  formule aux 2825 nuclides mesures.
""")
    print("  Repartition des ~200 MeV souvent cites :")
    print(f"  {'forme':<44} {'MeV':>7} {'recuperable':>13}")
    tot_fission = rec_total = 0.0
    for nom, e, recuperable, de_fission in BILAN_ENERGIE:
        if de_fission:
            tot_fission += e
        if recuperable:
            rec_total += e
        print(f"  {nom:<44} {e:>7.1f} {'oui' if recuperable else 'NON':>13}")
    print(f"  {'-' * 66}")
    print(f"  {'energie de fission proprement dite':<44} {tot_fission:>7.1f}")
    print(f"  {'total recuperable dans le coeur':<44} {rec_total:>7.1f}")
    print(f"""
  Les antineutrinos ({12.0:.0f} MeV, soit {12.0 / tot_fission * 100:.0f} % du total) traversent la Terre sans
  interagir : cette energie est definitivement perdue. C'est pourquoi on
  retient ~200 MeV recuperables et non ~207 MeV liberes.

  ORDRE DE GRANDEUR UTILE : 1 g d'U235 fissionne libere environ 1 MW-jour,
  soit l'equivalent de 2 a 3 tonnes de charbon. Le rapport est de l'ordre
  du million, et c'est simplement le rapport entre energies nucleaires
  (MeV) et chimiques (eV).
""")


# ==========================================================================
def etape4_ralentissement():
    titre(4, "POURQUOI ET COMMENT RALENTIR LES NEUTRONS ?")
    print("""
  Un neutron nait a ~2 MeV. Or la section de fission de l'U235 vaut 583
  barns a 0.025 eV, contre environ 1 barn a 2 MeV : elle est ~500 fois plus
  grande une fois le neutron thermalise. Sur uranium peu enrichi, il faut
  donc ralentir, et vite, car entre les deux se trouvent les resonances de
  capture de l'U238 (6 a 200 eV), qui sont des pieges.

  Un choc elastique fait gagner en moyenne xi en lethargie u = ln(E0/E),
  independamment de l'energie. D'ou le nombre de chocs necessaires :
""")
    print(f"  {'noyau':>8} {'A':>8} {'xi':>9} {'chocs 2 MeV -> 0.025 eV':>25}")
    for nom, A in (("H", 1.008), ("D", 2.014), ("C", 12.011),
                   ("O", 15.999), ("U238", 238.05)):
        print(f"  {nom:>8} {A:>8.2f} {xi_moyen(A):>9.4f} "
              f"{nombre_de_chocs(A):>25.0f}")
    print("""
  L'hydrogene est imbattable : xi = 1 exactement, 18 chocs suffisent. La
  raison est mecanique -- un neutron et un proton ont la meme masse, donc
  un choc frontal transfere TOUTE l'energie. C'est le billard. Sur l'U238,
  le neutron rebondit comme une bille sur un mur : 2171 chocs.

  Mais ralentir vite ne suffit pas : il faut ralentir SANS ABSORBER. Le
  vrai critere est le rapport de moderation xi*Sigma_s / Sigma_a.
""")
    donnees = [("eau legere H2O", "H", 2 * N_H2O, N_H2O),
               ("eau lourde D2O", "D", 2 * N_D2O, N_D2O),
               ("graphite", "C", densite_atomique(1.60, 12.011), 0.0)]
    print(f"  {'moderateur':<18} {'xi*Sigma_s [1/cm]':>18} {'rapport xi*Ss/Sa':>19}")
    for nom, cle, N, _ in donnees:
        print(f"  {nom:<18} {pouvoir_ralentissement(cle, N):>18.4f} "
              f"{rapport_moderation(cle, N):>19.0f}")
    print("""
  Le classement s'INVERSE. L'eau legere ralentit 6 fois plus vite que l'eau
  lourde (0.987 contre 0.162 par cm, noyau moderateur seul), mais son
  rapport de moderation est 76 fois plus mauvais (62 contre 4720), parce
  que l'hydrogene capture des neutrons (H + n -> deuterium).

  Toute l'architecture des filieres decoule de ce seul tableau :

  * EAU LEGERE : ralentit vite -> reseau compact, coeur petit, forte
    densite de puissance. Mais elle absorbe -> il FAUT enrichir. C'est le
    REP.
  * EAU LOURDE : n'absorbe presque pas -> marche a l'uranium NATUREL,
    donc pas d'usine d'enrichissement. Mais elle ralentit mal -> il en faut
    beaucoup, le reseau est tres ouvert et le coeur enorme. C'est CANDU.
  * GRAPHITE : intermediaire, coeur encore plus grand. Filieres UNGG,
    RBMK, et les reacteurs a haute temperature.

  Ce n'est pas un choix technique de detail : c'est ce qui decide si un
  pays a besoin ou non d'une usine d'enrichissement.
""")


# ==========================================================================
def _quatre_facteurs(enrichissement, Vm_Vf, moderateur="H", Nm=None,
                     heterogene=True):
    """Calcule (k_inf, eta, epsilon, p, f) pour un reseau donne."""
    Nm = N_H2O if Nm is None else Nm
    e = enrichissement
    fv_comb = 1.0 / (1.0 + Vm_Vf)
    fv_mod = Vm_Vf / (1.0 + Vm_Vf)

    N_U = N_U_UO2 * fv_comb                     # homogeneise
    N_mod = 2 * Nm * fv_mod
    N_O = Nm * fv_mod + 2 * N_U_UO2 * fv_comb   # oxygene de l'eau + de l'UO2

    sigma_a_comb = (e * NUCLIDES["U235"]["sa"] + (1 - e) * NUCLIDES["U238"]["sa"])
    f = facteur_f(N_U, sigma_a_comb, [(moderateur, N_mod), ("O", N_O)])

    if heterogene:
        xi_Ss_mod = (pouvoir_ralentissement(moderateur, 2 * Nm)
                     + pouvoir_ralentissement("O", Nm))
        p, I = facteur_p_heterogene(N_U_UO2 * (1 - e), fv_comb, fv_mod,
                                    xi_Ss_mod, RAYON_PASTILLE)
    else:
        xi_Ss = (pouvoir_ralentissement(moderateur, N_mod)
                 + pouvoir_ralentissement("O", N_O))
        Ss = (N_mod * NUCLIDES[moderateur]["ss"] * BARN
              + N_O * NUCLIDES["O"]["ss"] * BARN)
        p, I = facteur_p_homogene(N_U * (1 - e), Ss, xi_Ss)

    eta, eps = facteur_eta(e), 1.03
    return k_infini_quatre_facteurs(eta, eps, p, f), eta, eps, p, f, I


def etape5_quatre_facteurs():
    titre(5, "LE BILAN COMPLET : LA FORMULE DES QUATRE FACTEURS")
    print("""
      k_inf = eta * epsilon * p * f

  Le cycle du neutron thermique, lu a l'envers :
    f       utilisation thermique : il est absorbe dans le COMBUSTIBLE
            plutot que dans l'eau ou les structures
    p       echappement aux resonances : il a survecu a la traversee des
            resonances de l'U238 en ralentissant
    epsilon fission rapide : quelques fissions sur U238 l'ont precede
    eta     reproduction : son absorption produit eta nouveaux neutrons

  eta ne depend QUE de l'enrichissement, pas de la geometrie. C'est le
  plafond absolu du reacteur :
""")
    print(f"  {'enrichissement':>16} {'eta':>9}")
    for e in (0.0072, 0.02, 0.03, 0.045, 0.20, 1.0):
        etq = " (naturel)" if abs(e - 0.0072) < 1e-9 else (
            " (REP)" if abs(e - 0.03) < 1e-9 else "")
        print(f"  {e * 100:>14.2f} % {facteur_eta(e):>9.4f}{etq}")

    print("""
  HOMOGENE CONTRE HETEROGENE. Si l'on dissolvait l'uranium dans l'eau au
  lieu de le mettre en crayons, que se passerait-il ? Comparaison a 3 %
  d'enrichissement :
""")
    print(f"  {'Vm/Vf':>6} {'p homogene':>12} {'p heterogene':>14} "
          f"{'k_inf hom.':>12} {'k_inf het.':>12}")
    for v in (1.0, 1.5, 2.0, 3.0, 4.0):
        kh, _, _, ph, _, _ = _quatre_facteurs(0.03, v, heterogene=False)
        ke, _, _, pe, _, _ = _quatre_facteurs(0.03, v, heterogene=True)
        print(f"  {v:>6.1f} {ph:>12.4f} {pe:>14.4f} {kh:>12.4f} {ke:>12.4f}")
    print("""
  Le reseau heterogene gagne systematiquement, et c'est POUR CELA que le
  combustible est en crayons. Deux auto-protections se cumulent :

  * ENERGETIQUE : a l'energie exacte d'une resonance, les premiers atomes
    d'U238 absorbent tout, creusent le flux, et les suivants ne voient plus
    rien. L'integrale de resonance effective tombe de ~280 barns a
    dilution infinie a ~23 barns.
  * SPATIALE (heterogene seulement) : les neutrons ralentissent dans l'EAU,
    loin de l'uranium. Ils ne traversent la zone dangereuse 6-200 eV qu'a
    l'exterieur du crayon, ou il n'y a pas d'U238 pour les capturer.

  Le premier reacteur de Fermi, en 1942, etait un empilement heterogene de
  blocs d'uranium dans du graphite. Homogeneise, il n'aurait jamais
  diverge.
""")


# ==========================================================================
def etape6_moderation_et_uranium_naturel():
    titre(6, "SOUS-MODERATION ET URANIUM NATUREL")
    print("""
  Combien d'eau faut-il ? Deux effets opposes quand on en ajoute :
    p augmente  (on ralentit mieux, moins de captures resonnantes)
    f diminue   (l'eau absorbe une part croissante des neutrons)
  Il existe donc un optimum.
""")
    grille = np.linspace(0.5, 12.0, 90)
    ks = [_quatre_facteurs(0.03, v)[0] for v in grille]
    i_max = int(np.argmax(ks))
    print(f"  {'Vm/Vf':>7} {'p':>8} {'f':>8} {'k_inf':>9}")
    for v in (1.0, 2.0, 3.0, 4.4, 6.0, 8.0, 12.0):
        k, _, _, p, f, _ = _quatre_facteurs(0.03, v)
        note = "  <-- optimum" if abs(v - 4.4) < 0.05 else (
            "  <-- REP reel" if abs(v - 2.0) < 0.05 else "")
        print(f"  {v:>7.1f} {p:>8.4f} {f:>8.4f} {k:>9.4f}{note}")
    print(f"""
  Optimum calcule : Vm/Vf = {grille[i_max]:.2f}, k_inf = {ks[i_max]:.4f}.
  Or un REP fonctionne vers Vm/Vf = 2, nettement A GAUCHE de l'optimum.
  Ce n'est pas une contrainte subie, c'est un CHOIX DE SURETE.

  Un REP est volontairement SOUS-MODERE. Si la temperature monte, l'eau se
  dilate, Vm/Vf diminue -- et comme on est a gauche du maximum, k_inf
  DIMINUE. La puissance baisse d'elle-meme. Le coefficient de temperature
  moderateur est negatif, le reacteur est intrinsequement stable.

  A droite de l'optimum, ce serait l'inverse : dilatation -> plus de
  reactivite -> plus de puissance -> plus de dilatation. Emballement. Un
  coeur sur-modere serait interdit.

  C'est aussi la vraie raison de la limite en bore vue dans
  NEUTRONIQUE_COEUR.md : trop de bore rend ce coefficient positif.
""")
    print("  URANIUM NATUREL : quelle filiere fonctionne ?\n")
    print(f"  {'cas':<32} {'Vm/Vf':>7} {'p':>8} {'f':>8} {'k_inf':>9}")
    cas = [("U 3.0 % + H2O (REP)", 0.03, 2.0, "H", N_H2O),
           ("U naturel + H2O", 0.0072, 2.0, "H", N_H2O),
           ("U naturel + H2O, tres dilue", 0.0072, 6.0, "H", N_H2O),
           ("U naturel + D2O, reseau serre", 0.0072, 2.0, "D", N_D2O),
           ("U naturel + D2O (CANDU)", 0.0072, 16.0, "D", N_D2O)]
    for nom, e, v, mod, Nm in cas:
        k, _, _, p, f, _ = _quatre_facteurs(e, v, moderateur=mod, Nm=Nm)
        verdict = "  DIVERGE" if k > 1.0 else "  impossible"
        print(f"  {nom:<32} {v:>7.1f} {p:>8.4f} {f:>8.4f} {k:>9.4f}{verdict}")
    print("""
  Resultat central : l'uranium naturel avec de l'eau ordinaire donne
  k_inf < 1, quelle que soit la quantite d'eau. Le reacteur est
  IMPOSSIBLE, et aucun raffinement de conception n'y changera rien. Il faut
  soit enrichir, soit changer de moderateur.

  Avec de l'eau lourde et un reseau tres ouvert (Vm/Vf ~ 16), k_inf > 1 :
  c'est CANDU, qui se passe d'usine d'enrichissement. Noter la ligne
  precedente : la MEME eau lourde dans un reseau serre donne k_inf = 0.32.
  On ne change pas de moderateur sans redessiner tout le reseau.
""")


# ==========================================================================
def etape7_bouclage():
    titre(7, "BOUCLAGE : DU MICROSCOPIQUE AUX CONSTANTES A DEUX GROUPES")
    print("""
  coeur.py utilisait Sigma_a2 = 0.082 et nu_Sigma_f2 = 0.1375 sans jamais
  dire d'ou elles venaient. On peut maintenant les estimer.

  Premiere tentative naive : Sigma = N * sigma, avec sigma pris a 2200 m/s
  et N homogeneise sur l'assemblage.
""")
    e, Vm_Vf = 0.03, 2.0
    fv_comb = 1.0 / (1.0 + Vm_Vf)
    N_U = N_U_UO2 * fv_comb
    sa_2200 = e * NUCLIDES["U235"]["sa"] + (1 - e) * NUCLIDES["U238"]["sa"]
    sf_2200 = e * NUCLIDES["U235"]["sf"]
    Sa_naif = N_U * sa_2200 * BARN
    Sf_naif = N_U * NUCLIDES["U235"]["nu"] * sf_2200 * BARN

    print(f"    N_U homogeneise      = {N_U:.3e} /cm3")
    print(f"    Sigma_a2  (naif)     = {Sa_naif:.4f} 1/cm   "
          f"[modele coeur.py : 0.082]")
    print(f"    nuSigma_f2 (naif)    = {Sf_naif:.4f} 1/cm   "
          f"[modele coeur.py : 0.1375]")

    T = 580.0
    g = np.sqrt(np.pi) / 2 * np.sqrt(293.6 / T)
    print(f"""
  Facteur 2 d'ecart. Ce n'est pas une erreur, il manque deux etapes.

  (a) MOYENNE SUR LE SPECTRE. La valeur 2200 m/s est une convention, pas
      une moyenne. Les neutrons thermiques suivent une distribution de
      Maxwell a la temperature du moderateur. Pour une section en 1/v, la
      moyenne vaut sigma(2200) * (sqrt(pi)/2) * sqrt(293.6/T).

      A T = {T:.0f} K (REP en fonctionnement) : facteur {g:.4f}
""")
    print(f"    Sigma_a2  corrige    = {Sa_naif * g:.4f} 1/cm   "
          f"[modele coeur.py : 0.082]")
    print(f"    nuSigma_f2 corrige   = {Sf_naif * g:.4f} 1/cm   "
          f"[modele coeur.py : 0.1375]")
    print("""
  (b) FACTEUR DE DESAVANTAGE. L'homogeneisation ci-dessus suppose le flux
      uniforme. Il ne l'est pas : le crayon absorbe, donc le flux thermique
      y est DEPRIME par rapport a l'eau. Les sections homogeneisees
      correctes sont ponderees par le flux, pas par le volume. Le rapport
      phi_moderateur / phi_combustible (facteur de desavantage) vaut
      typiquement 1.1 a 1.3 dans un REP, ce qui reduit encore les sections
      combustible homogeneisees.

  Apres correction spectrale, on est dans un facteur 1.3 des valeurs de
  coeur.py, avec un calcul de coin de table. C'est l'ordre de grandeur
  attendu -- et c'est precisement pour cela que la generation des
  constantes de groupe est un metier a part entiere (calcul de reseau :
  APOLLO, CASMO, WIMS), et non une multiplication N*sigma.

  RETENIR : le chemin complet est
      donnees nucleaires evaluees (JEFF, ENDF/B)
        -> traitement des resonances et autoprotection
        -> calcul de reseau 2D fin (assemblage, multigroupe)
        -> condensation a 2 groupes + homogeneisation ponderee par le flux
        -> constantes utilisees par coeur.py
  Chaque fleche est un domaine de recherche.
""")


# ==========================================================================
def graphiques():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"  (graphiques ignores : {exc})")
        return

    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))

    # (a) courbe B/A
    ax = axes[0][0]
    A = np.arange(2, 260)
    Z = np.round(A / (1.98 + 0.0155 * A ** (2.0 / 3.0)))   # vallee de stabilite
    ax.plot(A, liaison_par_nucleon(A, Z), lw=1.8, color="#1f4e79")
    for nom, a, z in (("He4", 4, 2), ("Ni62", 62, 28), ("U235", 235, 92)):
        ax.plot(a, liaison_par_nucleon(a, z), "o", color="crimson")
        ax.annotate(nom, (a, liaison_par_nucleon(a, z)),
                    textcoords="offset points", xytext=(6, -10), fontsize=9)
    ax.axvline(62, color="grey", ls=":", lw=1)   # Ni-62 : vrai maximum mesure
    ax.set_xlabel("nombre de masse A")
    ax.set_ylabel("B/A [MeV]")
    ax.set_title("Energie de liaison par nucleon\n(fusion a gauche, fission a droite)")
    ax.grid(alpha=0.3)

    # (b) nombre de chocs
    ax = axes[0][1]
    noms = ["H", "D", "C", "O", "U238"]
    masses = [1.008, 2.014, 12.011, 15.999, 238.05]
    ax.bar(noms, [nombre_de_chocs(m) for m in masses], color="#1f4e79")
    ax.set_yscale("log")
    ax.set_ylabel("chocs pour 2 MeV -> 0.025 eV")
    ax.set_title("Efficacite de ralentissement")
    ax.grid(alpha=0.3, axis="y")
    for i, m in enumerate(masses):
        ax.text(i, nombre_de_chocs(m) * 1.15, f"{nombre_de_chocs(m):.0f}",
                ha="center", fontsize=9)

    # (c) courbe de moderation
    ax = axes[1][0]
    v = np.linspace(0.5, 12, 120)
    for e, coul in ((0.0072, "#c00000"), (0.02, "#e08214"), (0.03, "#1f4e79")):
        k = [_quatre_facteurs(e, x)[0] for x in v]
        ax.plot(v, k, color=coul, lw=1.8, label=f"{e * 100:.2f} % U235")
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.axvline(2.0, color="green", ls=":", lw=1.6)
    ax.annotate("REP\n(sous-modere)", (2.0, 0.72), color="green", fontsize=9,
                ha="center")
    kk = [_quatre_facteurs(0.03, x)[0] for x in v]
    ax.plot(v[int(np.argmax(kk))], max(kk), "k*", ms=13)
    ax.annotate("optimum", (v[int(np.argmax(kk))], max(kk)),
                textcoords="offset points", xytext=(8, 6), fontsize=9)
    ax.set_xlabel("rapport moderateur / combustible  Vm/Vf")
    ax.set_ylabel("k_inf")
    ax.set_title("Courbe de moderation (eau legere)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # (d) les quatre facteurs
    ax = axes[1][1]
    ens = np.linspace(0.005, 0.05, 60)
    res = [_quatre_facteurs(e, 2.0) for e in ens]
    ax.plot(ens * 100, [r[1] for r in res], label="eta (reproduction)", lw=1.8)
    ax.plot(ens * 100, [r[3] for r in res], label="p (resonances)", lw=1.8)
    ax.plot(ens * 100, [r[4] for r in res], label="f (utilisation th.)", lw=1.8)
    ax.plot(ens * 100, [r[0] for r in res], "k", label="k_inf", lw=2.4)
    ax.axhline(1.0, color="grey", ls="--", lw=1)
    ax.axvline(0.72, color="#c00000", ls=":", lw=1.6)
    ax.annotate("naturel", (0.72, 2.0), color="#c00000", fontsize=9)
    ax.set_xlabel("enrichissement [% U235]")
    ax.set_title("Les quatre facteurs (Vm/Vf = 2)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("De l'uranium au reacteur", fontsize=14)
    fig.tight_layout()
    fig.savefig("fission_bilan.png", dpi=130)
    print("  --> fission_bilan.png ecrit")


# ==========================================================================
def main():
    etape1_pourquoi_la_fission_libere()
    etape2_pourquoi_u235()
    etape3_bilan_energie()
    etape4_ralentissement()
    etape5_quatre_facteurs()
    etape6_moderation_et_uranium_naturel()
    etape7_bouclage()
    print("=" * 74)
    graphiques()


if __name__ == "__main__":
    main()
