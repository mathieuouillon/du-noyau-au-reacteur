"""
barrieres.py -- Energie de deformation d'une goutte chargee jusqu'a la
scission, et barrieres de fission.

LE PROBLEME
-----------
Le modele de frdm.py ne connait que les ellipsoides : il ne peut pas former
de col, donc pas decrire le chemin qui mene a la scission. Ici, on prend une
famille de formes a symetrie axiale assez riche pour aller de la sphere
jusqu'a deux fragments relies par un col qui se pince.

LA FAMILLE DE FORMES (c, h, alpha)
----------------------------------
Inspiree de la parametrisation "Funny Hills" de M. Brack et al.,
Rev. Mod. Phys. 44 (1972) 320. En coordonnees cylindriques (rho, z), unite
de longueur R0 (rayon de la sphere de meme volume) :

    rho^2(z) = (c^2 - z^2) * (A + B z^2/c^2 + alpha z/c),   -c <= z <= c

  * c     : elongation (demi-longueur ; c = 1 pour la sphere) ;
  * h     : formation du col, via B = 2h + (c - 1)/2 ;
  * alpha : asymetrie gauche-droite (masse des fragments) ;
  * A     : fixe par la conservation du volume, A = 1/c^3 - B/5
            (le terme en alpha, impair en z, ne change pas le volume).

c = 1, h = alpha = 0 donne la sphere ; h = -(c-1)/4 donne B = 0, donc un
spheroide de demi-axes c et c^(-1/2). La forme se pince (scission) quand
rho^2 s'annule a l'interieur de ]-c, c[.

LES TROIS ENERGIES DE FORME, RAPPORTEES A LA SPHERE
---------------------------------------------------
On ecrit f(z) = rho^2(z). Tout s'exprime avec f, f' et f'', ce qui evite les
singularites aux pointes (rho -> 0 mais rho*rho' = f'/2 reste fini).

  Surface   : S = 2 pi int sqrt(f + f'^2/4) dz              B_s = S / 4pi
  Courbure  : int (k1 + k2) dS = 2 pi int [1 - (2 f f'' - f'^2)/(4 f + f'^2)] dz
              rapportee a la sphere (8 pi)                  B_k
  Coulomb   : pour une charge uniforme, l'integrale de volume a 6 dimensions
              se ramene a une integrale DE SURFACE :
                  E_C = (rho_q^2/2) int_V int_V d3r d3r' / |r - r'|
                      = -(rho_q^2/4) oint oint (dS . dS') |r - r'|
              (utiliser lap|r-r'| = 2/|r-r'| puis deux fois Gauss).
              Avec la symetrie axiale il reste 3 integrales (z, z', angle
              relatif psi), sans singularite.                 B_c = E_C / E_C(sphere)

L'ENERGIE MACROSCOPIQUE
-----------------------
    E_def(forme) = E_s (B_s - 1) + E_k (B_k - 1) + E_c (B_c - 1)

avec deux jeux de parametres publies (K. Pomorski, J. Dudek, Phys. Rev. C 67
(2003) 044316, table I et eq. 3.5) :

  * MS-LD : goutte liquide de W.D. Myers, W.J. Swiatecki, Nucl. Phys. 81
    (1966) 1 ; Ark. Fys. 36 (1967) 343. Pas de terme de courbure.
  * LSD   : "Lublin-Strasbourg Drop" de Pomorski et Dudek, avec un terme de
    courbure en A^(1/3). Ajustee sur les masses seulement, elle reproduit
    pourtant les barrieres des noyaux Z > 70 a 0,9 MeV pres.

Unites : MeV, fm. e^2 = 1,44 MeV fm.
"""

import numpy as np

E2 = 1.44  # e^2 en MeV fm

# Parametres (Pomorski & Dudek 2003). Coulomb : (3/5) e^2 Z^2 / (r0 A^(1/3)).
PARAMETRES = {
    "MS-LD": dict(b_surf=18.56, k_surf=1.79, b_cour=0.0, k_cour=0.0, r0=1.2049),
    "LSD":   dict(b_surf=16.9707, k_surf=2.2938, b_cour=3.8602, k_cour=-2.3764,
                  r0=1.21725),
}


# ==========================================================================
# Formes
# ==========================================================================
class Forme:
    """Forme (c, h, alpha) ; f(z) = rho^2(z) et ses derivees."""

    def __init__(self, c, h=0.0, alpha=0.0):
        self.c, self.h, self.alpha = float(c), float(h), float(alpha)
        self.B = 2 * h + (c - 1) / 2
        self.A = 1 / c ** 3 - self.B / 5

    def _g(self, z):
        u = z / self.c
        return self.A + self.B * u * u + self.alpha * u

    def f(self, z):
        return (self.c ** 2 - z * z) * self._g(z)

    def df(self, z):
        c, u = self.c, z / self.c
        g = self._g(z)
        dg = (2 * self.B * u + self.alpha) / c
        return -2 * z * g + (c * c - z * z) * dg

    def d2f(self, z):
        c, u = self.c, z / self.c
        g = self._g(z)
        dg = (2 * self.B * u + self.alpha) / c
        d2g = 2 * self.B / c ** 2
        return -2 * g - 4 * z * dg + (c * c - z * z) * d2g

    def valide(self, n=401):
        """Forme d'un seul tenant : rho^2 > 0 strictement a l'interieur."""
        z = np.linspace(-self.c, self.c, n)[1:-1]
        return bool(np.all(self.f(z) > 0))

    def rayon_col(self, n=801):
        """Rayon du col (minimum local interieur de rho), en R0 ; NaN si la
        forme n'a pas de col."""
        z = np.linspace(-self.c, self.c, n)[1:-1]
        f = self.f(z)
        i = np.where((f[1:-1] < f[:-2]) & (f[1:-1] < f[2:]))[0] + 1
        if len(i) == 0:
            return float("nan")
        return float(np.sqrt(max(f[i].min(), 0.0)))

    def volume(self, n=2001):
        z = np.linspace(-self.c, self.c, n)
        return float(np.pi * np.trapezoid(self.f(z), z))

    def quadrupole(self, n=2001):
        """Moment quadrupolaire de la distribution uniforme, en unites de
        (3/4pi) R0^2 * volume : Q2 = int (2 z'^2 - rho^2) dV, z' mesure depuis
        le centre de masse."""
        z = np.linspace(-self.c, self.c, n)
        f = self.f(z)
        V = np.trapezoid(f, z)
        zc = np.trapezoid(z * f, z) / V
        zz = z - zc
        # int dV (2 z^2 - rho^2) = pi int dz [2 z^2 f - f^2/2]
        return float(np.trapezoid(2 * zz ** 2 * f - f * f / 2, z) / V)


# ==========================================================================
# Fonctions de forme B_s, B_k, B_c
# ==========================================================================
def _noeuds(c, n):
    """Changement de variable z = c sin(t) : il absorbe le comportement en
    racine carree de rho aux pointes, et la quadrature de Gauss-Legendre en t
    converge alors tres vite."""
    t, w = np.polynomial.legendre.leggauss(n)
    t = t * np.pi / 2
    w = w * np.pi / 2
    z = c * np.sin(t)
    dz = w * c * np.cos(t)
    return z, dz


def fonctions_forme(forme, n_z=96, n_psi=96):
    """(B_s, B_k, B_c) de la forme, rapportees a la sphere de meme volume."""
    z, dz = _noeuds(forme.c, n_z)
    f = np.maximum(forme.f(z), 0.0)
    fp = forme.df(z)
    fpp = forme.d2f(z)

    S = 2 * np.pi * np.sum(np.sqrt(f + fp * fp / 4) * dz)
    Bs = S / (4 * np.pi)

    courb = 1 - (2 * f * fpp - fp * fp) / (4 * f + fp * fp)
    Bk = 2 * np.pi * np.sum(courb * dz) / (8 * np.pi)

    # Coulomb : -(1/4) * 2pi * int dz int dz' int dpsi [rho rho' cos psi + f'f''/4] D
    rho = np.sqrt(f)
    psi = (np.arange(n_psi) + 0.5) * 2 * np.pi / n_psi     # regle du point milieu
    cpsi = np.cos(psi)
    total = 0.0
    for i in range(len(z)):
        rr = rho[i] * rho                                   # (n_z,)
        D2 = (f[i] + f)[:, None] - 2 * rr[:, None] * cpsi[None, :] \
            + ((z[i] - z) ** 2)[:, None]
        D = np.sqrt(np.maximum(D2, 0.0))
        integ = (rr[:, None] * cpsi[None, :] + (fp[i] * fp / 4)[:, None]) * D
        total += dz[i] * np.sum(integ.sum(axis=1) * dz) * (2 * np.pi / n_psi)
    Ec = -0.25 * 2 * np.pi * total                          # rho_q = 1, R0 = 1
    Bc = Ec / (16 * np.pi ** 2 / 15)                        # sphere : 16 pi^2 R^5 / 15
    return float(Bs), float(Bk), float(Bc)


# ==========================================================================
# Energie macroscopique
# ==========================================================================
def energies_sphere(Z, N, jeu="LSD"):
    """(E_s, E_k, E_c) de la sphere, en MeV, pour le jeu de parametres."""
    p = PARAMETRES[jeu]
    A = Z + N
    I = (N - Z) / A
    Es = p["b_surf"] * (1 - p["k_surf"] * I * I) * A ** (2 / 3)
    Ek = p["b_cour"] * (1 - p["k_cour"] * I * I) * A ** (1 / 3)
    Ec = 0.6 * E2 * Z * Z / (p["r0"] * A ** (1 / 3))
    return Es, Ek, Ec


def fissilite(Z, N, jeu="LSD"):
    """x = E_c / (2 E_s) (Bohr et Wheeler 1939)."""
    Es, _, Ec = energies_sphere(Z, N, jeu)
    return Ec / (2 * Es)


def energie_deformation(Z, N, B, jeu="LSD"):
    """E_def = E_s (B_s - 1) + E_k (B_k - 1) + E_c (B_c - 1), en MeV,
    a partir d'un triplet (B_s, B_k, B_c)."""
    Es, Ek, Ec = energies_sphere(Z, N, jeu)
    Bs, Bk, Bc = B
    return Es * (Bs - 1) + Ek * (Bk - 1) + Ec * (Bc - 1)


# ==========================================================================
# Carte (c, h) et col de la barriere
# ==========================================================================
class Carte:
    """Fonctions de forme tabulees sur une grille (c, h) a alpha fixe. Elles
    ne dependent pas du noyau : une seule carte sert pour tous les noyaux."""

    def __init__(self, cs, hs, alpha=0.0, **kw):
        self.cs, self.hs, self.alpha = np.asarray(cs), np.asarray(hs), alpha
        self.B = np.full((len(cs), len(hs), 3), np.nan)
        self.col = np.full((len(cs), len(hs)), np.nan)
        for i, c in enumerate(cs):
            for j, h in enumerate(hs):
                fo = Forme(c, h, alpha)
                if fo.valide():
                    self.B[i, j] = fonctions_forme(fo, **kw)
                    self.col[i, j] = fo.rayon_col()

    def energie(self, Z, N, jeu="LSD"):
        Es, Ek, Ec = energies_sphere(Z, N, jeu)
        return (Es * (self.B[..., 0] - 1) + Ek * (self.B[..., 1] - 1)
                + Ec * (self.B[..., 2] - 1))

    def vallee(self, Z, N, jeu="LSD"):
        """Pour chaque elongation c, le minimum en h : la vallee de fission.
        Retourne (energie, h du minimum)."""
        E = self.energie(Z, N, jeu)
        Ev = np.nanmin(E, axis=1)
        hv = self.hs[np.nanargmin(np.where(np.isnan(E), np.inf, E), axis=1)]
        return Ev, hv

    def barriere(self, Z, N, jeu="LSD"):
        """Hauteur de la barriere macroscopique : le maximum, le long de c, du
        minimum en h (c'est le point selle de la surface). Interpolation
        parabolique autour du maximum. Retourne (B_f, c_selle, h_selle)."""
        Ev, hv = self.vallee(Z, N, jeu)
        ok = np.isfinite(Ev)
        c, E, h = self.cs[ok], Ev[ok], hv[ok]
        i = int(np.argmax(E))
        if 0 < i < len(E) - 1:
            a, b, d = E[i - 1], E[i], E[i + 1]
            den = a - 2 * b + d
            if den < 0:
                x = 0.5 * (a - d) / den
                return float(b - 0.25 * (a - d) * x), float(c[i] + x * (c[1] - c[0])), float(h[i])
        return float(E[i]), float(c[i]), float(h[i])


def barriere_reduite(carte, x):
    """Barriere de la goutte liquide pure (sans courbure) en unites de E_s,
    pour un parametre de fissilite x quelconque :
        E_def / E_s = (B_s - 1) + 2 x (B_c - 1).
    Le point selle est le maximum le long de c du minimum en h."""
    e = (carte.B[..., 0] - 1) + 2 * x * (carte.B[..., 2] - 1)
    vallee = np.nanmin(e, axis=1)
    i = int(np.nanargmax(vallee))
    if 0 < i < len(vallee) - 1 and np.all(np.isfinite(vallee[i - 1:i + 2])):
        a, b, d = vallee[i - 1], vallee[i], vallee[i + 1]
        den = a - 2 * b + d
        if den < 0:
            u = 0.5 * (a - d) / den
            return float(b - 0.25 * (a - d) * u), float(carte.cs[i] + u * (carte.cs[1] - carte.cs[0]))
    return float(vallee[i]), float(carte.cs[i])


def bohr_wheeler(x):
    """Developpement de Bohr et Wheeler (1939), valable pour x proche de 1 :
    E_f / E_s = (98/135) (1 - x)^3."""
    return 98 / 135 * (1 - np.asarray(x)) ** 3


# ==========================================================================
# Penetrabilite d'une barriere parabolique (Hill et Wheeler 1953)
# ==========================================================================
def penetrabilite(E, Bf, hbar_omega):
    """T(E) = 1 / (1 + exp(2 pi (Bf - E) / hbar_omega)) : probabilite de
    traverser une barriere parabolique de hauteur Bf et de courbure
    hbar_omega (D.L. Hill, J.A. Wheeler, Phys. Rev. 89 (1953) 1102)."""
    return 1.0 / (1.0 + np.exp(2 * np.pi * (Bf - np.asarray(E)) / hbar_omega))
