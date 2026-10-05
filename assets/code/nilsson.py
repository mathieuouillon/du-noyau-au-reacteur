"""
nilsson.py -- Niveaux d'energie individuels des nucleons dans un noyau deforme.

C'est la piece MICROSCOPIQUE du modele macroscopique-microscopique : au lieu
de compter des nucleons de valence (M3), on calcule les vrais niveaux.

LE MODELE (S.G. Nilsson, Mat. Fys. Medd. Dan. Vid. Selsk. 29, no. 16, 1955)
--------------------------------------------------------------------------
Chaque nucleon se deplace dans un potentiel moyen cree par tous les autres :
un oscillateur harmonique a symetrie axiale, plus un couplage spin-orbite
et une correction qui aplatit le fond du puits.

    h = p^2/2m + (m/2) [ w_perp^2 (x^2 + y^2) + w_z^2 z^2 ]
        - kappa * hbar*w0_sph * [ 2 l.s + mu (l^2 - <l^2>_N) ]

  * deformation delta : w_z^2 = w0^2 (1 - 4 delta/3), w_perp^2 = w0^2 (1 + 2 delta/3)
    delta > 0 : noyau allonge (prolate, ballon de rugby)
    delta < 0 : noyau aplati (oblate, galette)
  * conservation du volume : w_perp^2 w_z = cte, d'ou
        w0(delta) = w0_sph * (1 - 4 delta^2/3 - 16 delta^3/27)^(-1/6)
  * l.s : couplage spin-orbite (Mayer et Jensen, 1949, prix Nobel 1963). C'est
    LUI qui produit les nombres magiques 28, 50, 82, 126 ; un oscillateur
    pur donnerait 2, 8, 20, 40, 70, 112.
  * mu l^2 : rend le potentiel plus plat au centre, plus proche d'un puits
    de Woods-Saxon realiste.

Avec la decomposition
    (m/2)[w_perp^2 (x^2+y^2) + w_z^2 z^2] = (m/2) w0^2 r^2 - (2/3) delta m w0^2 r^2 P2(cos theta)
on diagonalise dans la base spherique |N l j Omega> de l'oscillateur
isotrope de frequence w0(delta). Omega (projection de j sur l'axe) et la
parite sont conserves ; chaque niveau est doublement degenere (+-Omega).

ECART A LA FORMULATION D'ORIGINE. Nilsson ecrit l.s et l^2 en coordonnees
"etirees" (l_t) qui suivent la deformation. Nous les ecrivons en
coordonnees spheriques ordinaires et gardons TOUT le couplage
quadrupolaire, y compris entre couches N et N+-2. Les deux formulations
coincident a petite deformation et divergent progressivement au-dela de
|delta| ~ 0.3 : c'est une approximation assumee de ce code.

PARAMETRES kappa, mu : voir PARAMETRES ci-dessous, avec leur provenance.

Unites : toutes les energies sont en unites de hbar*w0_sph (oscillateur
spherique), ce qui rend le spectre INDEPENDANT du noyau : il ne depend que de
delta et du type de nucleon. Le passage en MeV se fait par
    hbar*w0 = 41 A^(-1/3) * (1 +- (N-Z)/3A) MeV     (+ neutrons, - protons)
(Bohr et Mottelson, Nuclear Structure vol. I ; correction d'isospin de
Nilsson et al. 1969).
"""

from functools import lru_cache
import numpy as np
from sympy import S
from sympy.physics.wigner import clebsch_gordan

# ==========================================================================
# Parametres (kappa, mu) par couche majeure N de l'oscillateur
# ==========================================================================
# Jeu "standard" de T. Bengtsson et I. Ragnarsson, Nucl. Phys. A436 (1985) 14,
# lui-meme issu de S.G. Nilsson et al., Nucl. Phys. A131 (1969) 1.
#
# Provenance de CHAQUE valeur (voir MODELE_FRDM.md) :
#   N = 2      : verifie (table standard reproduite dans arXiv nucl-th/0002063)
#   N = 3,4,5  : verifie (table standard reproduite dans arXiv nucl-th/9810061)
#   N = 6 n    : verifie (kappa = 0.062, mu = 0.34, cite par plusieurs travaux)
#   N = 6 p,
#   N = 7 n    : jeu ajuste pour la region Z ~ 100 (arXiv 1110.6134, table I),
#                faute d'avoir pu verifier la valeur standard
#   N = 0, 1   : prolonges depuis N = 2 (non verifie ; couches profondes,
#                influence faible sauf pour les noyaux legers)
#   N >= 7 p,
#   N >= 8 n   : prolonges depuis la derniere couche connue (non verifie ;
#                couches tres au-dessus du niveau de Fermi, qui ne servent
#                qu'au lissage de Strutinsky)
PARAMETRES = {
    "proton": {0: (0.105, 0.00), 1: (0.105, 0.00), 2: (0.105, 0.00),
               3: (0.090, 0.30), 4: (0.065, 0.57), 5: (0.060, 0.65),
               6: (0.057, 0.654)},
    "neutron": {0: (0.105, 0.00), 1: (0.105, 0.00), 2: (0.105, 0.00),
                3: (0.090, 0.25), 4: (0.070, 0.39), 5: (0.062, 0.43),
                6: (0.062, 0.34), 7: (0.0634, 0.318)},
}


def kappa_mu(type_nucleon, N):
    table = PARAMETRES[type_nucleon]
    return table[min(N, max(table))]


# ==========================================================================
# Elements de matrice
# ==========================================================================
_RHO = np.linspace(0.0, 14.0, 6001)


@lru_cache(maxsize=None)
def _radiale(n, l):
    """Fonction radiale normalisee de l'oscillateur (b = 1), echantillonnee.
    R_nl(rho) ~ rho^l exp(-rho^2/2) L_n^(l+1/2)(rho^2), n = nombre de noeuds."""
    from scipy.special import eval_genlaguerre
    f = _RHO ** l * np.exp(-_RHO ** 2 / 2) * eval_genlaguerre(n, l + 0.5, _RHO ** 2)
    f = f / np.sqrt(np.trapezoid(f ** 2 * _RHO ** 2, _RHO))
    return f


@lru_cache(maxsize=None)
def rho2(N1, l1, N2, l2):
    """<N1 l1 | rho^2 | N2 l2>, calcule par quadrature numerique."""
    f1 = _radiale((N1 - l1) // 2, l1)
    f2 = _radiale((N2 - l2) // 2, l2)
    return float(np.trapezoid(f1 * f2 * _RHO ** 4, _RHO))


@lru_cache(maxsize=None)
def _cg(j1, m1, j2, m2, J, M):
    return float(clebsch_gordan(S(j1), S(j2), S(J), S(m1), S(m2), S(M)))


@lru_cache(maxsize=None)
def p2(l1, j1x2, l2, j2x2, om_x2):
    """<l1 j1 Omega | P2(cos theta) | l2 j2 Omega> pour des etats couples l+s.
    j et Omega passes en double (entiers) pour l'exactitude.
    On developpe sur les etats |l m_l>|sigma> et on utilise
      <l1 m|P2|l2 m> = sqrt((2 l2+1)/(2 l1+1)) <l2 0 2 0|l1 0><l2 m 2 0|l1 m>
    """
    j1, j2, om = S(j1x2) / 2, S(j2x2) / 2, S(om_x2) / 2
    total = 0.0
    for sig in (S(1) / 2, -S(1) / 2):
        m = om - sig
        if abs(m) > l1 or abs(m) > l2:
            continue
        c1 = _cg(l1, m, S(1) / 2, sig, j1, om)
        c2 = _cg(l2, m, S(1) / 2, sig, j2, om)
        if c1 == 0 or c2 == 0:
            continue
        ang = (np.sqrt((2 * l2 + 1) / (2 * l1 + 1))
               * _cg(l2, 0, 2, 0, l1, 0) * _cg(l2, m, 2, 0, l1, m))
        total += c1 * c2 * ang
    return total


def base(N_max, om_x2):
    """Etats |N l j Omega> de la base, pour un Omega donne."""
    etats = []
    for N in range(N_max + 1):
        for l in range(N % 2, N + 1, 2):
            for jx2 in (2 * l - 1, 2 * l + 1):
                if jx2 > 0 and jx2 >= om_x2:
                    etats.append((N, l, jx2))
    return etats


# ==========================================================================
# Diagonalisation
# ==========================================================================
def _facteur_volume(delta):
    return (1 - 4 * delta ** 2 / 3 - 16 * delta ** 3 / 27) ** (-1 / 6)


@lru_cache(maxsize=None)
def _blocs_constants(N_max, om_x2, type_nucleon):
    """Parties de la matrice independantes de delta, pour un bloc Omega."""
    etats = base(N_max, om_x2)
    n = len(etats)
    diag_osc = np.array([N + 1.5 for N, _, _ in etats])
    so = np.zeros(n)
    Q = np.zeros((n, n))
    for i, (N, l, jx2) in enumerate(etats):
        kap, mu = kappa_mu(type_nucleon, N)
        j = jx2 / 2
        ls = (j * (j + 1) - l * (l + 1) - 0.75) / 2
        so[i] = -kap * (2 * ls + mu * (l * (l + 1) - N * (N + 3) / 2))
    for i, (N1, l1, j1) in enumerate(etats):
        for k, (N2, l2, j2) in enumerate(etats):
            if k < i or abs(N1 - N2) > 2 or abs(l1 - l2) > 2:
                continue
            v = rho2(N1, l1, N2, l2) * p2(l1, j1, l2, j2, om_x2)
            Q[i, k] = Q[k, i] = v
    return diag_osc, so, Q


def spectre(delta, type_nucleon, N_max=14):
    """Niveaux individuels (en unites de hbar*w0_sph), tries.
    Chaque valeur represente un niveau doublement degenere."""
    c = _facteur_volume(delta)
    niveaux = []
    for om_x2 in range(1, 2 * N_max + 2, 2):
        diag_osc, so, Q = _blocs_constants(N_max, om_x2, type_nucleon)
        H = c * (np.diag(diag_osc) - (2.0 / 3.0) * delta * Q) + np.diag(so)
        niveaux.append(np.linalg.eigvalsh(H))
    return np.sort(np.concatenate(niveaux))


def hbar_omega(Z, N, type_nucleon):
    """hbar*w0 en MeV pour ce type de nucleon dans le noyau (Z, N)."""
    A = Z + N
    signe = 1.0 if type_nucleon == "neutron" else -1.0
    return 41.0 * A ** (-1.0 / 3.0) * (1.0 + signe * (N - Z) / (3.0 * A))


def delta_vers_beta2(delta):
    """Conversion approchee vers le beta2 des tables experimentales :
    delta ~ (3/2) sqrt(5/4pi) beta2 ~ 0.946 beta2 (au premier ordre)."""
    return np.asarray(delta) / (1.5 * np.sqrt(5.0 / (4.0 * np.pi)))
