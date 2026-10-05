"""
modeles_masse.py -- Modeles d'energie de liaison au-dela de la goutte liquide.

Une hierarchie, chaque etage ajoutant une piece de physique :

  M0  goutte liquide, coefficients de manuel (ceux de fission.py)
  M1  goutte liquide, coefficients REAJUSTES sur AME2020
  M2  goutte liquide etendue : symetrie de surface, echange coulombien,
      terme de Wigner
  M3  M2 + correction de couches (comptage des nucleons de valence)
  M4  M3 + correction statistique du residu (regression a noyau)

M0 a M3 sont LINEAIRES dans leurs coefficients : l'ajustement est un simple
moindres carres, sans aucun hasard. M4 est un apprentissage automatique.

Convention : B > 0 pour un noyau lie, en MeV. La metrique standard du domaine
est l'ecart RMS sur l'energie de liaison TOTALE B (pas B/A).
"""

import numpy as np
from fission import energie_liaison as B_manuel

# Nombres magiques, plus 184 (fermeture de couche predite, jamais observee)
MAGIQUES = np.array([0, 2, 8, 20, 28, 50, 82, 126, 184])


# ==========================================================================
# Variables de base
# ==========================================================================
def _base(Z, N):
    Z = np.asarray(Z, float)
    N = np.asarray(N, float)
    A = Z + N
    I = N - Z
    pair = np.where((Z % 2 == 0) & (N % 2 == 0), 1.0,
                    np.where((Z % 2 == 1) & (N % 2 == 1), -1.0, 0.0))
    return Z, N, A, I, pair


# ==========================================================================
# Le terme de couches
# ==========================================================================
def valence(n):
    """Mesure de l'eloignement d'une couche fermee, pour un type de nucleon.

    Si n est dans la couche [M_bas, M_haut[ de degenerescence D = M_haut-M_bas,
    avec nu = n - M_bas particules de valence et D - nu trous :

        x = nu * (D - nu) / D

    x vaut 0 sur une couche fermee et est maximal a mi-couche. C'est la
    forme de Dieperink et Van Isacker (Eur. Phys. J. A 42, 2009), qui
    reprend l'idee de Casten : ce qui compte n'est pas n mais le nombre de
    particules OU de trous de valence, le plus petit des deux.
    """
    n = np.asarray(n, float)
    i = np.searchsorted(MAGIQUES, n, side="right") - 1
    bas = MAGIQUES[i]
    haut = MAGIQUES[np.minimum(i + 1, len(MAGIQUES) - 1)]
    D = haut - bas
    nu = n - bas
    return nu * (D - nu) / D


# ==========================================================================
# Matrices de conception (une colonne = un terme = un coefficient)
# ==========================================================================
def colonnes_M1(Z, N):
    """Goutte liquide a 5 termes. Z(Z-1) plutot que Z^2 : un proton ne se
    repousse pas lui-meme."""
    Z, N, A, I, pair = _base(Z, N)
    return {
        "volume":    A,
        "surface":   -A ** (2 / 3),
        "Coulomb":   -Z * (Z - 1) / A ** (1 / 3),
        "asymetrie": -I ** 2 / A,
        "appariement": pair / A ** 0.5,
    }


def colonnes_M2(Z, N):
    """Goutte liquide etendue.

    symetrie de surface : l'exces de neutrons coute MOINS cher en surface
                          qu'au coeur (il y forme une peau de neutrons). Le
                          coefficient de symetrie effectif d'un noyau fini
                          est a_sym(A) = J - Q*A^(-1/3), plus petit que
                          celui de la matiere nucleaire infinie J.
                          -> + (N-Z)^2 / A^(4/3)
    echange coulombien  : correction quantique a la repulsion (principe de
                          Pauli entre protons) -> + Z^(4/3) / A^(1/3)
    Wigner              : surliaison des noyaux N ~ Z, ou neutrons et
                          protons occupent les memes orbitales -> -|N-Z|/A
    """
    Z_, N_, A, I, pair = _base(Z, N)
    c = colonnes_M1(Z, N)
    c["sym. surface"] = I ** 2 / A ** (4 / 3)
    c["echange coul."] = Z_ ** (4 / 3) / A ** (1 / 3)
    c["Wigner"] = -np.abs(I) / A
    return c


def colonnes_M3(Z, N):
    """M2 + couches. S = x_p + x_n ; trois termes :
        a1 * S          effet lineaire de l'eloignement des couches
        a2 * S^2        courbure
        a3 * x_p * x_n  interaction proton-neutron de valence, moteur de la
                        deformation des noyaux en milieu de couche
    """
    c = colonnes_M2(Z, N)
    xp, xn = valence(Z), valence(N)
    S = xp + xn
    c["couches S"] = S
    c["couches S^2"] = S ** 2
    c["couches p-n"] = xp * xn
    return c


# ==========================================================================
# Ajustement lineaire
# ==========================================================================
class ModeleLineaire:
    """Modele B = somme_k a_k * terme_k, ajuste par moindres carres."""

    def __init__(self, nom, fabrique):
        self.nom = nom
        self.fabrique = fabrique
        self.coefs = None

    def _X(self, Z, N):
        cols = self.fabrique(Z, N)
        self.termes = list(cols)
        return np.column_stack([cols[k] for k in self.termes])

    def ajuster(self, Z, N, B):
        X = self._X(Z, N)
        self.coefs, *_ = np.linalg.lstsq(X, B, rcond=None)
        return self

    def predire(self, Z, N):
        return self._X(Z, N) @ self.coefs

    def nombre_parametres(self):
        return len(self.coefs)


class ModeleManuel:
    """M0 : la formule de fission.py, coefficients de manuel, rien a ajuster."""
    nom = "M0  goutte liquide (manuel)"

    def ajuster(self, Z, N, B):
        return self

    def predire(self, Z, N):
        Z = np.asarray(Z, float)
        N = np.asarray(N, float)
        return B_manuel(Z + N, Z)

    def nombre_parametres(self):
        return 5


# ==========================================================================
# M4 : correction statistique du residu
# ==========================================================================
class ModeleCorrige:
    """M3 + regression a noyau gaussien (kernel ridge) sur le residu.

    Idee : M3 capture la physique connue ; ce qu'il rate n'est pas du bruit
    mais une structure lisse sur la carte (Z, N) -- deformations, couches
    secondaires. Une regression a noyau apprend cette structure a partir des
    voisins. Chaque prediction est une moyenne ponderee des erreurs de M3 sur
    les noyaux proches, la ponderation decroissant comme exp(-d^2 / 2 l^2).

    l (longueur de correlation, en nombre de nucleons) et lambda
    (regularisation) sont choisis par validation croisee SUR LES DONNEES
    D'ENTRAINEMENT seulement.

    C'est la methode de correction par fonctions de base radiales (RBF)
    utilisee dans la litterature pour corriger FRDM, HFB ou les modeles WS
    (Wang et Liu, 2011).
    """
    nom = "M4  M3 + regression a noyau"

    def __init__(self, longueurs=(1.5, 2.0, 3.0, 4.0, 6.0), lambdas=(0.03, 0.3)):
        self.base = ModeleLineaire("M3", colonnes_M3)
        self.longueurs = longueurs
        self.lambdas = lambdas

    @staticmethod
    def _noyau(P, Q, l):
        d2 = (P[:, :1] - Q[:, 0]) ** 2 + (P[:, 1:] - Q[:, 1]) ** 2
        return np.exp(-d2 / (2 * l * l))

    def _apprendre(self, P, r, l, lam):
        """Resout (K + lam I) alpha = r par Cholesky, K = L L^T.
        On utilise la factorisation de numpy : scipy.linalg.cho_factor
        (lower=True) echoue a tort sur cette matrice avec SciPy 1.18."""
        from scipy.linalg import solve_triangular
        K = self._noyau(P, P, l)
        K[np.diag_indices_from(K)] += lam
        L = np.linalg.cholesky(K)
        y = solve_triangular(L, r, lower=True, check_finite=False)
        return solve_triangular(L.T, y, lower=False, check_finite=False)

    def ajuster(self, Z, N, B, graine=0):
        self.base.ajuster(Z, N, B)
        r = B - self.base.predire(Z, N)
        P = np.column_stack([Z, N]).astype(float)

        # choix des hyperparametres par validation croisee a 3 plis
        rng = np.random.default_rng(graine)
        plis = rng.integers(0, 3, len(P))
        meilleur = (np.inf, None)
        for l in self.longueurs:
            for lam in self.lambdas:
                err = []
                for k in range(3):
                    a, t = plis != k, plis == k
                    alpha = self._apprendre(P[a], r[a], l, lam)
                    err.append(r[t] - self._noyau(P[t], P[a], l) @ alpha)
                rms = np.sqrt(np.mean(np.concatenate(err) ** 2))
                if rms < meilleur[0]:
                    meilleur = (rms, (l, lam))
        self.l, self.lam = meilleur[1]
        self.P = P
        self.alpha = self._apprendre(P, r, self.l, self.lam)
        return self

    def predire(self, Z, N):
        Q = np.column_stack([Z, N]).astype(float)
        return self.base.predire(Z, N) + self._noyau(Q, self.P, self.l) @ self.alpha

    def nombre_parametres(self):
        # un poids par noyau d'entrainement : c'est la nature d'une methode
        # a noyau, et la raison pour laquelle elle peut sur-apprendre
        return self.base.nombre_parametres() + len(self.alpha) + 2


def tous_les_modeles():
    return [
        ModeleManuel(),
        ModeleLineaire("M1  goutte liquide reajustee", colonnes_M1),
        ModeleLineaire("M2  goutte liquide etendue", colonnes_M2),
        ModeleLineaire("M3  M2 + couches", colonnes_M3),
        ModeleCorrige(),
    ]
