"""
orbitales.py -- Les "orbitales" des nucleons, et la forme du noyau.

L'ANALOGIE AVEC L'ATOME
-----------------------
Dans l'atome d'hydrogene, l'electron occupe des orbitales 1s, 2p, 3d... dont
on trace la densite de probabilite |psi|^2. Dans le noyau, chaque nucleon
occupe de meme une orbitale du potentiel moyen. Mais trois differences :

  1. le potentiel n'est pas en 1/r (Coulomb) mais ressemble a un puits :
     les niveaux sont ceux d'un oscillateur, 1s, 1p, 1d, 2s, 1f... ;
  2. le spin-orbite est enorme (quelques MeV, contre des micro-eV dans
     l'atome) : on doit coupler l et s en j des le depart. D'ou les
     notations 1d5/2, 1g9/2, 1i13/2 ;
  3. le potentiel lui-meme peut etre DEFORME. l n'est alors plus un bon
     nombre quantique : seule reste la projection Omega de j sur l'axe du
     noyau. Les orbitales sont des melanges de plusieurs (l, j).

DE LA BASE AU NUAGE
-------------------
nilsson.py diagonalise h dans la base |N l j Omega>. Une orbitale est un
vecteur de coefficients c_i. En position, avec m = Omega - sigma :

    psi_sigma(r, theta, phi) = somme_i c_i R_nl(r/b) <l m 1/2 sigma | j Omega>
                               Y_lm(theta, phi) / b^(3/2)

    |psi|^2 = somme_sigma |psi_sigma|^2      (independant de phi : symetrie axiale)

b = sqrt(hbar / m w0) est la longueur de l'oscillateur, environ
1,006 A^(1/6) fm (de hbar*w0 = 41 A^(-1/3) MeV). La densite etant de
revolution autour de l'axe z, il suffit de la calculer dans le demi-plan
(rho, z) : une coupe suffit, comme pour les orbitales atomiques.

FORME DU NOYAU : on somme |psi|^2 sur toutes les orbitales occupees (deux
nucleons par niveau, +Omega et -Omega). On obtient la densite de nucleons
en fm^-3.
"""

import numpy as np
from scipy.special import eval_genlaguerre, gammaln, sph_harm_y
from sympy import S
from nilsson import _blocs_constants, base, _facteur_volume, hbar_omega, _cg

HBARC = 197.327          # MeV fm
MN = 938.92              # MeV, masse moyenne d'un nucleon
LETTRES = "spdfghijklmno"   # notation nucleaire usuelle (1j15/2 pour l = 7)


def longueur_oscillateur(Z, N, type_nucleon, delta=0.0):
    """b en fm pour ce type de nucleon, a la deformation delta.
    La base est celle de l'oscillateur de frequence w0(delta) = c(delta) w0_sph."""
    hw = hbar_omega(Z, N, type_nucleon) * _facteur_volume(delta)
    return HBARC / np.sqrt(MN * hw)


def nom_spherique(N, l, jx2):
    nr = (N - l) // 2
    return f"{nr + 1}{LETTRES[l]}{jx2}/2"


def etats(delta, type_nucleon, N_max=12):
    """Toutes les orbitales : liste triee de dict(energie, om_x2, base, vecteur)."""
    c = _facteur_volume(delta)
    res = []
    for om_x2 in range(1, 2 * N_max + 2, 2):
        diag_osc, so, Q = _blocs_constants(N_max, om_x2, type_nucleon)
        H = c * (np.diag(diag_osc) - (2.0 / 3.0) * delta * Q) + np.diag(so)
        e, V = np.linalg.eigh(H)
        b = base(N_max, om_x2)
        for k in range(len(e)):
            res.append(dict(energie=e[k], om_x2=om_x2, base=b, vecteur=V[:, k]))
    res.sort(key=lambda s: s["energie"])
    for s in res:
        i = int(np.argmax(s["vecteur"] ** 2))
        s["parent"] = nom_spherique(*s["base"][i])
        s["poids_parent"] = float(s["vecteur"][i] ** 2)
        s["parite"] = "+" if s["base"][i][0] % 2 == 0 else "-"
    return res


def _R(n, l, x):
    """Fonction radiale normalisee de l'oscillateur (b = 1)."""
    norme = np.sqrt(2.0 * np.exp(gammaln(n + 1) - gammaln(n + l + 1.5)))
    return norme * x ** l * np.exp(-x * x / 2) * eval_genlaguerre(n, l + 0.5, x * x)


def densite_orbitale(etat, rho, z, b):
    """|psi|^2 (fm^-3) sur une grille (rho, z) en fm, pour UN nucleon."""
    r = np.sqrt(rho ** 2 + z ** 2)
    theta = np.arctan2(rho, z)
    x = r / b
    om = etat["om_x2"] / 2
    psi = {0.5: np.zeros_like(r), -0.5: np.zeros_like(r)}
    for c_i, (N, l, jx2) in zip(etat["vecteur"], etat["base"]):
        if abs(c_i) < 1e-6:
            continue
        rad = c_i * _R((N - l) // 2, l, x)
        for sig in (0.5, -0.5):
            m = om - sig
            if abs(m) > l:
                continue
            mi = int(round(m))
            cg = _cg(l, mi, S(1) / 2, S(int(2 * sig)) / 2, S(jx2) / 2, S(etat["om_x2"]) / 2)
            if cg == 0:
                continue
            Y = np.real(sph_harm_y(l, int(round(m)), theta, 0.0))
            psi[sig] = psi[sig] + rad * cg * Y
    return (psi[0.5] ** 2 + psi[-0.5] ** 2) / b ** 3


def densite_noyau(Z, N, delta, rho, z, N_max=12):
    """Densites de protons et de neutrons (fm^-3) : somme des orbitales occupees."""
    dens = {}
    for t, n in (("proton", Z), ("neutron", N)):
        b = longueur_oscillateur(Z, N, t, delta)
        E = etats(delta, t, N_max)
        d = np.zeros_like(rho)
        for k in range(n // 2):
            d += 2 * densite_orbitale(E[k], rho, z, b)
        if n % 2:
            d += densite_orbitale(E[n // 2], rho, z, b)
        dens[t] = d
    return dens


def rayon_rms(densite, rho, z):
    """Rayon quadratique moyen d'une densite de revolution, en fm."""
    dr = rho[0, 1] - rho[0, 0]
    dz = z[1, 0] - z[0, 0]
    poids = 2 * np.pi * rho * densite * dr * dz
    return float(np.sqrt(np.sum(poids * (rho ** 2 + z ** 2)) / np.sum(poids))), float(poids.sum())
