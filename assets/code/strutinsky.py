"""
strutinsky.py -- Des niveaux individuels a la correction de couches.

LE PROBLEME
-----------
Sommer les energies des niveaux occupes donne une energie totale tres
mauvaise : le potentiel de Nilsson est fait pour avoir les bons NIVEAUX, pas
la bonne energie totale. La goutte liquide, elle, donne la bonne energie
MOYENNE mais ignore les couches.

L'IDEE DE STRUTINSKY (V.M. Strutinsky, Nucl. Phys. A95 (1967) 420 ;
Nucl. Phys. A122 (1968) 1)
-----------------------------------------------------------------------------
On garde de la somme des niveaux uniquement sa partie FLUCTUANTE :

    delta_E_couches = somme des niveaux occupes  -  meme somme "lissee"

La version lissee est obtenue en etalant chaque niveau sur une largeur
gamma ~ hbar*w0, plus grande que l'ecart entre couches : les couches
disparaissent, il ne reste que la tendance moyenne. La difference est la
part "couches" de l'energie, a ajouter a la goutte liquide :

    E = E_goutte_liquide(forme) + delta_E_couches(forme) + delta_E_appariement(forme)

C'est le theoreme de Strutinsky : sous des hypotheses justifiees par
Brack et al. (Rev. Mod. Phys. 44 (1972) 320, la revue dite "Funny Hills"),
cette decomposition reproduit l'energie d'un calcul Hartree-Fock au premier
ordre dans les fluctuations.

LISSAGE
-------
Densite lissee : g~(e) = somme_i d_i (1/gamma) f((e - e_i)/gamma), avec
    f(x) = exp(-x^2)/sqrt(pi) * L_M^(1/2)(x^2)
gaussienne corrigee par un polynome de Laguerre generalise d'ordre M
(correction de courbure, M = 3). Le polynome garantit que le lissage
n'altere pas une densite deja lisse (polynome de degre <= 2M).
d_i = 2 (degenerescence de Kramers).

CONDITION DE PLATEAU : le resultat ne doit pratiquement pas dependre de
gamma dans une plage autour de 1.0-1.4 hbar*w0. C'est le test de validite
de la methode ; il est verifie dans etude_frdm.py.

APPARIEMENT (BCS)
-----------------
Les nucleons s'apparient par paires de moments opposes, comme les electrons
d'un supraconducteur (Bardeen, Cooper, Schrieffer 1957 ; application aux
noyaux : Bohr, Mottelson, Pines 1958). L'appariement adoucit les effets de
couches. On le traite en BCS, avec la meme logique que Strutinsky :

    delta_P = P(niveaux reels) - P~(densite lissee)

La force G est fixee, noyau par noyau, pour que le modele "uniforme" (densite
lisse g~ au niveau de Fermi) donne le gap moyen empirique
    Delta~ = 12 / sqrt(A) MeV          (Bohr et Mottelson, vol. I)
C'est la prescription de Brack et al. (1972). FRDM utilise une variante plus
elaboree (Lipkin-Nogami, pour corriger la non-conservation du nombre de
particules de BCS) : ecart assume par ce code.
"""

import numpy as np
from scipy.special import eval_genlaguerre

M_COURBURE = 3


def _f(x):
    return np.exp(-x * x) / np.sqrt(np.pi) * eval_genlaguerre(M_COURBURE, 0.5, x * x)


class Lissage:
    """Quantites lissees d'un spectre, tabulees sur une grille d'energie.

    Tout est en unites de hbar*w0 : la table vaut pour TOUS les noyaux.
    """

    def __init__(self, niveaux, gamma=1.2, pas=0.004):
        self.niveaux = np.sort(niveaux)
        self.gamma = gamma
        e = np.arange(self.niveaux[0] - 6 * gamma, self.niveaux[-1], pas)
        g = np.zeros_like(e)
        for bloc in np.array_split(self.niveaux, max(1, len(self.niveaux) // 64)):
            g += (2.0 / gamma) * _f((e[:, None] - bloc[None, :]) / gamma).sum(1)
        self.e, self.g = e, g
        # integrales cumulees : nombre de particules et energie lissee
        de = np.diff(e)
        self.Ncum = np.concatenate([[0], np.cumsum(0.5 * (g[1:] + g[:-1]) * de)])
        eg = e * g
        self.Ecum = np.concatenate([[0], np.cumsum(0.5 * (eg[1:] + eg[:-1]) * de)])
        # somme exacte des niveaux occupes (2 nucleons par niveau)
        self.somme = np.concatenate([[0], np.cumsum(2 * self.niveaux)])

    def fermi_lisse(self, n):
        """Niveau de Fermi lisse : la ou la densite lissee contient n nucleons."""
        return float(np.interp(n, self.Ncum, self.e))

    def densite(self, e):
        return float(np.interp(e, self.e, self.g))

    def correction_couches(self, n):
        """delta_E = somme des n premiers etats - energie lissee (en hbar*w0)."""
        exacte = self.somme[n // 2] + (self.niveaux[n // 2] if n % 2 else 0.0)
        lisse = float(np.interp(n, self.Ncum, self.Ecum))
        return exacte - lisse


# ==========================================================================
# Appariement
# ==========================================================================
def _bcs(eps, n_paires, G):
    """Resout les equations BCS sur des niveaux de paires eps (en MeV).

    Equations du gap et du nombre :
        2/G = somme_k 1/E_k,   n_paires = somme_k v_k^2,
        E_k = sqrt((eps_k - lambda)^2 + Delta^2),  v_k^2 = (1 - (eps_k-lambda)/E_k)/2
    Resolution emboitee et robuste : pour Delta donne, lambda par dichotomie
    sur l'equation du nombre (monotone) ; puis Delta par dichotomie sur
    l'equation du gap. Si le gap ne s'ouvre pas (G trop faible devant
    l'ecart des niveaux), il n'y a pas d'appariement : Delta = 0.
    Retourne (Delta, energie BCS sans le terme -G sum v^4).
    """
    from scipy.optimize import brentq

    def lam_de(D):
        f = lambda l: np.sum(0.5 * (1 - (eps - l) / np.sqrt((eps - l) ** 2 + D * D))) - n_paires
        return brentq(f, eps[0] - 50.0, eps[-1] + 50.0, xtol=1e-10)

    def gap(D):
        l = lam_de(D)
        return np.sum(1 / np.sqrt((eps - l) ** 2 + D * D)) - 2 / G

    D_min, D_max = 1e-4, 20.0
    if gap(D_min) <= 0:
        return 0.0, None
    D = D_max if gap(D_max) > 0 else brentq(gap, D_min, D_max, xtol=1e-8)
    lam = lam_de(D)
    E = np.sqrt((eps - lam) ** 2 + D * D)
    v2 = 0.5 * (1 - (eps - lam) / E)
    return D, float(np.sum(2 * eps * v2) - D * D / G)


def correction_appariement(lissage, n, hw, A, demi_fenetre=1.0):
    """delta_P en MeV, pour n nucleons d'un type, avec hbar*w0 = hw MeV.

    Fenetre d'appariement : niveaux a moins de demi_fenetre*hbar*w0 du niveau
    de Fermi lisse. Nombre impair : le dernier nucleon "bloque" son niveau,
    qui sort de l'appariement (c'est l'origine microscopique de l'effet
    pair-impair des masses).
    """
    Dt = 12.0 / np.sqrt(A)                         # gap moyen empirique (MeV)
    lam_t = lissage.fermi_lisse(n)
    rho = 0.5 * lissage.densite(lam_t) / hw         # densite de PAIRES, par MeV
    W = demi_fenetre * hw
    if rho <= 0:
        return 0.0
    G = 1.0 / (rho * np.arcsinh(W / Dt))            # modele uniforme -> Delta = Delta~
    P_moyen = rho * W * (W - np.sqrt(W * W + Dt * Dt))   # ~ -rho Delta~^2 / 2

    niv = lissage.niveaux * hw                      # MeV, niveaux de paires
    bloque = None
    n_paires = n // 2
    if n % 2:
        bloque = n_paires                           # niveau a moitie occupe
    dedans = np.where(np.abs(niv - lam_t * hw) <= W)[0]
    if bloque is not None:
        dedans = dedans[dedans != bloque]
    if len(dedans) < 2:
        return -P_moyen
    eps = niv[dedans]
    paires_dessous = int(np.sum(dedans < n_paires))  # paires occupees dans la fenetre
    if paires_dessous == 0 or paires_dessous >= len(eps):
        return -P_moyen
    D, E_bcs = _bcs(eps, paires_dessous, G)
    E_sans = float(np.sum(2 * eps[:paires_dessous]))
    P = 0.0 if E_bcs is None else E_bcs - E_sans
    return P - P_moyen
