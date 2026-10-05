"""
fission.py -- De la structure du noyau aux sections efficaces du reacteur.

On remonte a la source de tout ce qui precede : pourquoi l'uranium fissionne,
pourquoi seul l'U235 le fait avec des neutrons lents, combien d'energie sort,
et comment on passe des donnees microscopiques (barns, noyaux/cm3) aux
sections macroscopiques (1/cm) utilisees dans diffusion.py et coeur.py.

Unites : energies en MeV, sections microscopiques en barns (1 b = 1e-24 cm2),
densites atomiques en noyaux/cm3, sections macroscopiques en 1/cm.
"""

import numpy as np

BARN = 1.0e-24          # cm2
N_AVOGADRO = 6.02214e23


# ==========================================================================
# 1. Formule semi-empirique de masse (Bethe-Weizsacker)
# ==========================================================================
# Coefficients ajustes sur les masses mesurees, en MeV.
A_VOLUME = 15.75
A_SURFACE = 17.80
A_COULOMB = 0.711
A_ASYMETRIE = 23.70
A_APPARIEMENT = 11.18


def energie_liaison(A, Z):
    """Energie de liaison totale d'un noyau (A, Z), en MeV.

    B = a_V*A - a_S*A^(2/3) - a_C*Z^2/A^(1/3) - a_A*(A-2Z)^2/A + delta

    Chaque terme raconte quelque chose :
      volume     : chaque nucleon est lie a ses voisins -> proportionnel a A
      surface    : les nucleons du bord ont moins de voisins -> -A^(2/3)
      Coulomb    : les protons se repoussent -> -Z^2/A^(1/3). C'est CE terme
                   qui rend les noyaux lourds fragiles et rend la fission
                   exothermique.
      asymetrie  : il faut autant de n que de p (principe de Pauli)
      appariement: les nucleons s'apparient. Terme decisif pour comprendre
                   pourquoi l'U235 est fissile et pas l'U238.
    """
    A = np.asarray(A, dtype=float)
    Z = np.asarray(Z, dtype=float)
    N = A - Z

    B = (A_VOLUME * A
         - A_SURFACE * A ** (2.0 / 3.0)
         - A_COULOMB * Z ** 2 / A ** (1.0 / 3.0)
         - A_ASYMETRIE * (A - 2 * Z) ** 2 / A)

    # terme d'appariement : + si N et Z pairs, - si tous deux impairs, 0 sinon
    pair_Z = (Z % 2 == 0)
    pair_N = (N % 2 == 0)
    delta = np.where(pair_Z & pair_N, A_APPARIEMENT / np.sqrt(A),
                     np.where(~pair_Z & ~pair_N, -A_APPARIEMENT / np.sqrt(A), 0.0))
    return B + delta


def liaison_par_nucleon(A, Z):
    """B/A en MeV : la courbe qui explique fission ET fusion."""
    return energie_liaison(A, Z) / np.asarray(A, dtype=float)


def energie_separation_neutron(A, Z):
    """Energie a fournir pour arracher un neutron : S_n = B(A,Z) - B(A-1,Z).

    Lue a l'envers, c'est l'energie LIBEREE quand le noyau (A-1,Z) capture
    un neutron. C'est la grandeur qui decide de la fissilite.
    """
    return energie_liaison(A, Z) - energie_liaison(A - 1, Z)


# ==========================================================================
# 2. Bilan energetique de la fission
# ==========================================================================
def energie_fission(A_lourd, Z_lourd, A1, Z1, A2, Z2):
    """Energie liberee par une fission, en MeV.

    Q = B(fragment 1) + B(fragment 2) - B(noyau compose)

    Les neutrons libres ont une energie de liaison nulle, ils ne comptent
    donc pas dans le bilan. L'energie sort parce que les fragments, plus
    proches du fer, sont PLUS lies que l'uranium : la difference est
    liberee, essentiellement en energie cinetique de repulsion coulombienne
    des deux fragments.
    """
    return (energie_liaison(A1, Z1) + energie_liaison(A2, Z2)
            - energie_liaison(A_lourd, Z_lourd))


# Repartition mesuree des ~200 MeV d'une fission d'U235.
# Seule la partie "recuperable" chauffe l'eau : les antineutrinos
# traversent la Terre sans interagir.
# (libelle, MeV, recuperable, issu de la fission elle-meme)
BILAN_ENERGIE = [
    ("energie cinetique des fragments", 168.0, True,  True),
    ("energie cinetique des neutrons", 5.0, True,  True),
    ("rayons gamma prompts", 7.0, True,  True),
    ("beta des produits de fission", 8.0, True,  True),
    ("gamma des produits de fission", 7.0, True,  True),
    ("antineutrinos", 12.0, False, True),
    ("captures (n,gamma) hors fission", 10.0, True,  False),
]


# ==========================================================================
# 3. Donnees nucleaires : sections efficaces microscopiques
# ==========================================================================
# Valeurs a 2200 m/s (0.0253 eV), le point de reference conventionnel.
# sigma_a = absorption totale (capture + fission), sigma_f = fission,
# sigma_s = diffusion elastique, nu = neutrons par fission.
NUCLIDES = {
    #  nom        M       sigma_a  sigma_f  sigma_s   nu
    "U235":  dict(M=235.04, sa=681.0,  sf=583.0, ss=10.0, nu=2.42),
    "U238":  dict(M=238.05, sa=2.68,   sf=0.0,   ss=9.0,  nu=0.0),
    "Pu239": dict(M=239.05, sa=1017.0, sf=748.0, ss=8.0,  nu=2.87),
    "H":     dict(M=1.008,  sa=0.332,  sf=0.0,   ss=20.5, nu=0.0),
    "D":     dict(M=2.014,  sa=0.00052, sf=0.0,  ss=3.4,  nu=0.0),
    "O":     dict(M=15.999, sa=0.00019, sf=0.0,  ss=3.8,  nu=0.0),
    "C":     dict(M=12.011, sa=0.0035, sf=0.0,   ss=4.75, nu=0.0),
    "B10":   dict(M=10.013, sa=3840.0, sf=0.0,   ss=4.0,  nu=0.0),
    "Xe135": dict(M=134.91, sa=2.65e6, sf=0.0,   ss=4.0,  nu=0.0),
}


def densite_atomique(masse_volumique, masse_molaire, atomes_par_molecule=1):
    """N = rho * N_A / M, en noyaux/cm3.

    C'est le pont entre la chimie (g/cm3) et la neutronique (1/cm) :
        Sigma [1/cm] = N [1/cm3] * sigma [cm2]
    """
    return masse_volumique * N_AVOGADRO / masse_molaire * atomes_par_molecule


# ==========================================================================
# 4. Ralentissement des neutrons
# ==========================================================================
def xi_moyen(A):
    """Gain moyen de lethargie par choc elastique.

    La lethargie u = ln(E0/E) mesure le chemin parcouru vers les basses
    energies. Un choc fait gagner en moyenne xi, INDEPENDAMMENT de
    l'energie : c'est ce qui rend le comptage si simple.

        xi = 1 + (A-1)^2/(2A) * ln((A-1)/(A+1))

    Pour A = 1 (hydrogene) la formule est singuliere mais la limite vaut
    exactement 1 : un neutron peut perdre TOUTE son energie en un seul choc
    sur un proton, parce qu'ils ont la meme masse. Aucun autre noyau ne
    fait cela.
    """
    A = float(A)
    if abs(A - 1.0) < 1e-9:
        return 1.0
    return 1.0 + (A - 1.0) ** 2 / (2.0 * A) * np.log((A - 1.0) / (A + 1.0))


def nombre_de_chocs(A, E_debut=2.0e6, E_fin=0.0253):
    """Nombre moyen de chocs pour ralentir de E_debut a E_fin (en eV)."""
    return np.log(E_debut / E_fin) / xi_moyen(A)


def pouvoir_ralentissement(nom, N):
    """xi * Sigma_s : vitesse a laquelle le milieu ralentit les neutrons."""
    d = NUCLIDES[nom]
    return xi_moyen(d["M"]) * N * d["ss"] * BARN


def rapport_moderation(nom, N):
    """xi*Sigma_s / Sigma_a : le vrai critere de choix d'un moderateur.

    Ralentir vite ne suffit pas : il faut ralentir sans absorber. C'est ce
    rapport, et non le pouvoir de ralentissement seul, qui explique
    pourquoi l'eau lourde permet un reacteur a l'uranium naturel et pas
    l'eau ordinaire.
    """
    d = NUCLIDES[nom]
    return pouvoir_ralentissement(nom, N) / (N * d["sa"] * BARN)


# ==========================================================================
# 5. Formule des quatre facteurs
# ==========================================================================
def facteur_eta(enrichissement):
    """eta = neutrons produits par neutron thermique absorbe DANS le combustible.

        eta = nu * Sigma_f(comb) / Sigma_a(comb)

    Plafond absolu du reacteur : meme sans aucune fuite ni absorption
    parasite, on ne peut pas faire mieux. Depend uniquement de
    l'enrichissement, pas de la geometrie.
    """
    e = enrichissement
    u5, u8 = NUCLIDES["U235"], NUCLIDES["U238"]
    # proportions atomiques (approximation : enrichissement massique ~ atomique)
    f5, f8 = e, 1.0 - e
    sigma_f = f5 * u5["sf"] + f8 * u8["sf"]
    sigma_a = f5 * u5["sa"] + f8 * u8["sa"]
    return u5["nu"] * sigma_f / sigma_a


def facteur_f(N_comb, sigma_a_comb, autres):
    """f = utilisation thermique = Sigma_a(combustible) / Sigma_a(total).

    Fraction des neutrons thermiques qui finissent dans le combustible
    plutot que dans le moderateur, les structures ou le bore.
    """
    Sa_comb = N_comb * sigma_a_comb * BARN
    Sa_autres = sum(N * NUCLIDES[nom]["sa"] * BARN for nom, N in autres)
    return Sa_comb / (Sa_comb + Sa_autres)


def facteur_p_homogene(N_U238, Sigma_s_total, xi_Sigma_s):
    """Probabilite d'echapper aux resonances, melange HOMOGENE.

        p = exp( - N_28 * I_eff / (xi * Sigma_s) )

    I_eff est l'integrale de resonance effective de l'U238. A dilution
    infinie elle vaut ~280 barns, mais elle chute fortement par
    AUTO-PROTECTION : les premieres couches d'U238 absorbent tellement a
    l'energie d'une resonance qu'elles creusent le flux, et les atomes
    suivants ne voient plus rien a absorber. On utilise la correlation
    empirique classique :

        I_eff = 3.0 * sigma_p^0.55     (barns)

    ou sigma_p est la section de diffusion du milieu ramenee a un atome
    d'U238. C'est une CORRELATION, pas un calcul : la vraie integrale de
    resonance demande un traitement fin des resonances (Nordheim, methode
    des sous-groupes). Ordre de grandeur seulement.

    Arguments en unites macroscopiques (1/cm) ; la conversion en barns est
    faite ici, justement pour qu'on ne puisse pas l'oublier.
    """
    sigma_p = Sigma_s_total / (N_U238 * BARN)          # barns par atome d'U238
    I_eff = 3.0 * sigma_p ** 0.55                      # barns
    return float(np.exp(-N_U238 * I_eff * BARN / xi_Sigma_s)), I_eff


def facteur_p_heterogene(N_U238_comb, V_comb, V_mod, xi_Sigma_s_mod,
                         rayon_pastille, densite_UO2=10.4):
    """Probabilite d'echapper aux resonances, reseau HETEROGENE (crayons).

        p = exp( - N_28 * V_comb * I_eff / (xi*Sigma_s_mod * V_mod) )

    Ici N_28 est la densite dans le COMBUSTIBLE, pas dans le melange
    homogeneise. L'integrale de resonance suit la correlation classique
    en surface sur masse du crayon :

        I_eff = 4.45 + 26.6 * sqrt(S/M)     (barns, S/M en cm2/g)

    Pour un cylindre, S/M = 2 / (rayon * densite).

    Le point physique : dans un reseau, les neutrons ralentissent dans
    l'eau, LOIN de l'uranium. Ils ne traversent la zone dangereuse
    (6-200 eV) qu'a l'exterieur du crayon, ou il n'y a pas d'U238 pour les
    capturer. C'est l'auto-protection SPATIALE, qui s'ajoute a
    l'auto-protection energetique.
    """
    S_sur_M = 2.0 / (rayon_pastille * densite_UO2)     # cm2/g
    I_eff = 4.45 + 26.6 * np.sqrt(S_sur_M)             # barns
    expo = N_U238_comb * V_comb * I_eff * BARN / (xi_Sigma_s_mod * V_mod)
    return float(np.exp(-expo)), float(I_eff)


def k_infini_quatre_facteurs(eta, epsilon, p, f):
    """k_inf = eta * epsilon * p * f.

    Le cycle du neutron en quatre etapes, lu de droite a gauche :
      f       : il est absorbe dans le combustible plutot qu'ailleurs
      p       : il a survecu aux resonances de l'U238 en ralentissant
      epsilon : quelques fissions rapides sur l'U238 l'ont precede
      eta     : son absorption produit eta nouveaux neutrons
    """
    return eta * epsilon * p * f
