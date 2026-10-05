"""
frdm.py -- Modele macroscopique-microscopique, dans l'esprit de FRDM.

    B(Z, N) = B_goutte_spherique(Z, N)
              - min sur la forme delta de [ E_def_macro(delta)
                                           + dE_couches(delta)
                                           + dE_appariement(delta) ]

C'est l'architecture de FRDM (P. Moller, J.R. Nix, W.D. Myers,
W.J. Swiatecki, At. Data Nucl. Data Tables 59 (1995) 185 ; version 2012 :
P. Moller, A.J. Sierk, T. Ichikawa, H. Sagawa, ADNDT 109-110 (2016) 1),
en version tres simplifiee. Ce qui est garde, ce qui est simplifie :

  piece                  | FRDM                          | ce code
  -----------------------|-------------------------------|-------------------------
  macroscopique          | gouttelette a portee finie,   | goutte liquide etendue
                         | ~30 termes                    | (M2, 8 termes)
  formes                 | eps2, eps3, eps4, eps6, gamma | spheroides seulement (delta)
  potentiel individuel   | Yukawa replie (folded-Yukawa) | Nilsson (oscillateur
                         |                               | + l.s + l^2)
  correction de couches  | Strutinsky                    | Strutinsky (identique)
  appariement            | Lipkin-Nogami                 | BCS + Strutinsky
  parametres ajustes     | ~ 40, sur toute la table      | 8 coefficients lineaires

LA DEFORMATION MACROSCOPIQUE (N. Bohr, J.A. Wheeler, Phys. Rev. 56 (1939) 426)
---------------------------------------------------------------------------
Deformer une goutte a volume constant change deux energies :
  * la SURFACE augmente  -> cout
  * la repulsion COULOMB diminue (les charges s'eloignent) -> gain
Pour un ellipsoide de revolution d'excentricite e, ces rapports a la sphere
ont des expressions exactes (formules classiques, voir R.W. Hasse et
W.D. Myers, Geometrical Relationships of Macroscopic Nuclear Physics,
Springer 1988) ; a petite deformation
    B_s ~ 1 + (2/45) e^4,   B_c ~ 1 - (1/45) e^4
d'ou une energie de deformation proportionnelle a (2 E_s - E_c) = 2 E_s (1 - x)
ou x = E_c / (2 E_s) est le PARAMETRE DE FISSILITE. Pour x >= 1, la sphere
n'est plus stable : la goutte fissionne spontanement. L'uranium a x ~ 0.7 :
stable, mais pas de beaucoup. C'est l'explication de Bohr et Wheeler de la
fission, publiee neuf mois apres sa decouverte.

On prend comme forme de la goutte l'ellipsoide equipotentiel du potentiel de
Nilsson, de rapport d'axes a/b = w_perp/w_z : forme macroscopique et
potentiel microscopique sont ainsi decrits par le MEME delta.
"""

import numpy as np
from nilsson import spectre, hbar_omega
from strutinsky import Lissage, correction_appariement

GRILLE_DELTA = np.round(np.arange(-0.30, 0.4501, 0.025), 4)


# ==========================================================================
# Formes : ellipsoides de revolution
# ==========================================================================
def rapport_axes(delta):
    """a/b = (axe de symetrie)/(axe equatorial) = w_perp / w_z."""
    return np.sqrt((1 + 2 * delta / 3) / (1 - 4 * delta / 3))


def surface_coulomb(delta):
    """(B_s, B_c) : surface et energie coulombienne rapportees a la sphere,
    pour l'ellipsoide de revolution de meme volume (formules exactes)."""
    q = rapport_axes(delta)
    if abs(q - 1) < 1e-9:
        return 1.0, 1.0
    if q > 1:                                          # allonge (prolate)
        e = np.sqrt(1 - 1 / q ** 2)
        Bs = 0.5 * q ** (-2 / 3) * (1 + q * np.arcsin(e) / e)
        Bc = (1 - e * e) ** (1 / 3) / (2 * e) * np.log((1 + e) / (1 - e))
    else:                                              # aplati (oblate)
        e = np.sqrt(1 - q ** 2)
        Bs = 0.5 * q ** (-2 / 3) * (1 + q * q * np.arctanh(e) / e)
        Bc = (1 - e * e) ** (1 / 6) * np.arcsin(e) / e
    return float(Bs), float(Bc)


_BS_BC = {d: surface_coulomb(d) for d in GRILLE_DELTA}


# ==========================================================================
# Tables microscopiques (calculees une fois pour toutes)
# ==========================================================================
class TablesMicro:
    """Pour chaque delta de la grille et chaque type de nucleon, le spectre de
    Nilsson et son lissage de Strutinsky. En unites de hbar*w0, donc valables
    pour tous les noyaux."""

    def __init__(self, gamma=1.2, N_max=14):
        self.lissages = {}
        for d in GRILLE_DELTA:
            for t in ("proton", "neutron"):
                self.lissages[(d, t)] = Lissage(spectre(d, t, N_max), gamma)

    def energie_micro(self, Z, N, delta, appariement=True):
        """dE_couches + dE_appariement en MeV, protons + neutrons."""
        A = Z + N
        total = 0.0
        for t, n in (("proton", Z), ("neutron", N)):
            L = self.lissages[(delta, t)]
            hw = hbar_omega(Z, N, t)
            total += L.correction_couches(n) * hw
            if appariement:
                total += correction_appariement(L, n, hw, A)
        return total


# ==========================================================================
# Le modele
# ==========================================================================
def colonnes_macro(Z, N):
    """Termes de la goutte liquide spherique (ceux de M2, sans couches)."""
    Z = np.asarray(Z, float)
    N = np.asarray(N, float)
    A = Z + N
    I = N - Z
    pair = np.where((Z % 2 == 0) & (N % 2 == 0), 1.0,
                    np.where((Z % 2 == 1) & (N % 2 == 1), -1.0, 0.0))
    return np.column_stack([
        A,                                   # volume
        -A ** (2 / 3),                       # surface
        -Z * (Z - 1) / A ** (1 / 3),         # Coulomb
        -I ** 2 / A,                         # asymetrie
        pair / A ** 0.5,                     # appariement moyen
        I ** 2 / A ** (4 / 3),               # symetrie de surface
        Z ** (4 / 3) / A ** (1 / 3),         # echange coulombien
        -np.abs(I) / A,                      # Wigner
    ])


TERMES = ["volume", "surface", "Coulomb", "asymetrie", "appariement",
          "sym. surface", "echange coul.", "Wigner"]


class ModeleMacMic:
    """Ajustement auto-coherent :
      1. avec les coefficients courants, on calcule pour chaque noyau la
         forme delta* qui minimise E_def_macro + E_micro ;
      2. a forme fixee, B_exp + E(delta*) est lineaire dans les coefficients
         macroscopiques : moindres carres ;
      3. on recommence jusqu'a stabilite (3 ou 4 tours suffisent).
    La partie microscopique n'a AUCUN parametre ajuste : kappa, mu et
    l'appariement viennent de la litterature.
    """

    def __init__(self, tables=None, appariement=True, cache="micro.npz"):
        """tables : TablesMicro (calcul a la volee), ou None pour relire le
        cache produit par precalcul_micro.py."""
        self.tables = tables
        self.appariement = appariement
        self.coefs = None
        self._cache = {}
        if tables is None:
            d = np.load(cache)
            assert np.allclose(d["delta"], GRILLE_DELTA), "cache d'une autre grille"
            col = d["micro"] if appariement else d["couches_seules"]
            self._cache = {(int(z), int(n)): col[i]
                           for i, (z, n) in enumerate(zip(d["Z"], d["N"]))}

    def _micro(self, Z, N):
        cle = (int(Z), int(N))
        if cle not in self._cache:
            if self.tables is None:
                self.tables = TablesMicro()
            self._cache[cle] = np.array([
                self.tables.energie_micro(int(Z), int(N), d, self.appariement)
                for d in GRILLE_DELTA])
        return self._cache[cle]

    def energie_forme(self, Z, N, coefs=None):
        """E(delta) = E_def_macro(delta) + E_micro(delta), sur la grille."""
        c = self.coefs if coefs is None else coefs
        A = Z + N
        I = N - Z
        Es0 = c[1] * A ** (2 / 3) - c[5] * I ** 2 / A ** (4 / 3)
        Ec0 = c[2] * Z * (Z - 1) / A ** (1 / 3) - c[6] * Z ** (4 / 3) / A ** (1 / 3)
        macro = np.array([Es0 * (_BS_BC[d][0] - 1) + Ec0 * (_BS_BC[d][1] - 1)
                          for d in GRILLE_DELTA])
        return macro + self._micro(Z, N), macro

    def forme(self, Z, N):
        """(delta*, E(delta*)) avec interpolation parabolique autour du minimum."""
        E, _ = self.energie_forme(Z, N)
        i = int(np.argmin(E))
        if 0 < i < len(E) - 1:
            a, b, c = E[i - 1], E[i], E[i + 1]
            den = a - 2 * b + c
            if den > 0:
                x = 0.5 * (a - c) / den
                h = GRILLE_DELTA[1] - GRILLE_DELTA[0]
                return GRILLE_DELTA[i] + x * h, b - 0.25 * (a - c) * x
        return GRILLE_DELTA[i], E[i]

    def ajuster(self, Z, N, B, coefs_initiaux, tours=4, verbose=False):
        self.coefs = np.array(coefs_initiaux, float)
        X = colonnes_macro(Z, N)
        for k in range(tours):
            Ef = np.array([self.forme(z, n)[1] for z, n in zip(Z, N)])
            self.coefs, *_ = np.linalg.lstsq(X, B + Ef, rcond=None)
            if verbose:
                r = B - self.predire(Z, N)
                print(f"    tour {k + 1} : RMS = {np.sqrt(np.mean(r ** 2)):.3f} MeV")
        return self

    def predire(self, Z, N):
        Z = np.atleast_1d(Z)
        N = np.atleast_1d(N)
        Ef = np.array([self.forme(z, n)[1] for z, n in zip(Z, N)])
        return colonnes_macro(Z, N) @ self.coefs - Ef

    def deformations(self, Z, N):
        return np.array([self.forme(z, n)[0] for z, n in zip(np.atleast_1d(Z),
                                                             np.atleast_1d(N))])
