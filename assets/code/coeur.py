"""
coeur.py -- Modelisation d'un coeur de reacteur a eau pressurisee (REP).

S'appuie sur diffusion.py. On passe ici du "solveur de neutronique" au
"modele de coeur" : assemblages, plan de chargement, grappes de commande,
bore soluble, et les grandeurs que l'ingenieur de coeur regarde reellement
(reactivite, efficacite des grappes, facteur de point chaud).

ATTENTION SUR LES DONNEES
-------------------------
Les sections efficaces ci-dessous sont *representatives* d'un REP a deux
groupes, pas des donnees certifiees. Les ordres de grandeur sont bons
(k_inf, poids du bore, efficacite des grappes), mais ce n'est pas un
benchmark. Pour de la validation, utiliser une specification publiee
(benchmark IAEA 2D, BIBLIS) qui fournit ses propres sections.

Conventions : longueurs en cm, sections en 1/cm, groupe 0 = rapide.
"""

import numpy as np
from diffusion import Material, Mesh, solve

# --------------------------------------------------------------------------
# Geometrie type REP 900 MWe
# --------------------------------------------------------------------------
PAS_ASSEMBLAGE = 21.6      # cm, pas d'un assemblage 17x17
N_ASSEMBLAGES = 8          # par cote, en quart de coeur
EP_REFLECTEUR = 21.6 * 2   # cm d'eau + cloisonnement autour du coeur
HAUTEUR_ACTIVE = 366.0     # cm, hauteur active du coeur

# Le modele est 2D : il ne voit aucune fuite axiale. On la reintroduit par un
# buckling transverse B_z^2, qui ajoute une pseudo-absorption D*B_z^2 (voir
# Material.sigma_r). Si le flux axial est un cosinus, la fuite axiale vaut
# exactement D*B_z^2*phi : l'approximation est excellente pour un coeur nu
# axialement. Sans ce terme, k_eff est surestime d'environ 400 a 500 pcm
# (valeur mesuree par etude_coeur.py, pas une estimation).
BUCKLING_AXIAL = (np.pi / (HAUTEUR_ACTIVE + 4 * 0.37)) ** 2   # ~7.3e-5 1/cm2


# --------------------------------------------------------------------------
# Bibliotheque de sections efficaces a deux groupes
# --------------------------------------------------------------------------
def _combustible(nom, sigma_a2, nu_sigma_f2, kappa_sigma_f2):
    """Cree un combustible en ne faisant varier que le groupe thermique.

    C'est la ou l'enrichissement se voit le plus : plus d'U235 augmente
    nu_Sigma_f thermique (production) et, un peu, Sigma_a thermique
    (absorption). Le groupe rapide bouge tres peu d'une zone a l'autre.
    """
    return Material(
        nom,
        D=[1.43, 0.37],
        sigma_a=[0.0095, sigma_a2],
        sigma_s=[[0.0, 0.0175],
                 [0.0, 0.0]],
        nu_sigma_f=[0.0055, nu_sigma_f2],
        chi=[1.0, 0.0],
    )


def bibliotheque():
    """Les materiaux du coeur.

    Trois zones d'enrichissement croissant + un reflecteur. Le nombre de
    zones est un choix de conception : c'est le principal levier pour
    aplatir la puissance (voir etude_coeur.py).
    """
    mats = {
        # zone           Sigma_a2  nuSigma_f2  kappaSigma_f2
        "zone1": _combustible("zone1 (faible enr.)", 0.0800, 0.1250, 0.0400),
        "zone2": _combustible("zone2 (moyen enr.)", 0.0810, 0.1310, 0.0420),
        "zone3": _combustible("zone3 (fort enr.)", 0.0820, 0.1375, 0.0441),
        "reflecteur": Material(
            "reflecteur (eau + cloisonnement)",
            D=[1.50, 0.30],
            sigma_a=[0.0004, 0.0200],
            sigma_s=[[0.0, 0.0380],
                     [0.0, 0.0]],
            nu_sigma_f=[0.0, 0.0],
            chi=[0.0, 0.0],
        ),
    }
    # kappa*Sigma_f sert au calcul de puissance, pas a la neutronique.
    # On le range a part car Material ne le connait pas.
    mats["zone1"].kappa_sigma_f = np.array([0.0128, 0.0400])
    mats["zone2"].kappa_sigma_f = np.array([0.0131, 0.0420])
    mats["zone3"].kappa_sigma_f = np.array([0.0134, 0.0441])
    mats["reflecteur"].kappa_sigma_f = np.array([0.0, 0.0])
    return mats


def k_infini(m):
    """k_inf a deux groupes : milieu infini, aucune fuite.

    C'est le potentiel intrinseque du combustible. Le coeur reel aura
    toujours k_eff < k_inf a cause des fuites.
    """
    s12 = m.sigma_s[0, 1]
    return (m.nu_sigma_f[0] + m.nu_sigma_f[1] * s12 / m.sigma_a[1]) / (m.sigma_a[0] + s12)


# --------------------------------------------------------------------------
# Effets de controle de la reactivite
# --------------------------------------------------------------------------
# Section d'absorption thermique ajoutee par ppm de bore soluble.
# Le bore est un poison quasi purement thermique (1/v) : il n'agit
# pratiquement que sur le groupe 2.
# Valeur recalee sur la physique : N_B = 1e-6 * rho_eau * N_A / M_B noyaux/cm3
# par ppm, fois sigma_a(bore naturel, 2200 m/s) = 760 barns, corrige du
# spectre thermique REP (~0.55) et de la fraction volumique de moderateur
# dans l'assemblage homogeneise (~0.55). Donne ~9.2e-6.
SIGMA_A_PAR_PPM_BORE = 9.20e-6      # 1/cm par ppm

# Absorption thermique ajoutee par une grappe inseree (AIC ou B4C
# homogeneisee sur l'assemblage).
SIGMA_A_GRAPPE = 0.0300             # 1/cm


def avec_bore(mat, ppm):
    """Retourne une copie du materiau empoisonnee au bore soluble."""
    if ppm == 0.0:
        return mat
    m = Material(f"{mat.name} + {ppm:.0f} ppm B",
                 D=mat.D.copy(), sigma_a=mat.sigma_a.copy(),
                 sigma_s=mat.sigma_s.copy(),
                 nu_sigma_f=mat.nu_sigma_f.copy(), chi=mat.chi.copy())
    m.sigma_a[1] += ppm * SIGMA_A_PAR_PPM_BORE
    m.kappa_sigma_f = getattr(mat, "kappa_sigma_f", np.zeros(2))
    return m


def avec_grappe(mat):
    """Retourne une copie du materiau avec une grappe de commande inseree."""
    m = Material(f"{mat.name} + grappe",
                 D=mat.D.copy(), sigma_a=mat.sigma_a.copy(),
                 sigma_s=mat.sigma_s.copy(),
                 nu_sigma_f=mat.nu_sigma_f.copy(), chi=mat.chi.copy())
    m.sigma_a[1] += SIGMA_A_GRAPPE
    m.kappa_sigma_f = getattr(mat, "kappa_sigma_f", np.zeros(2))
    return m


# --------------------------------------------------------------------------
# Plans de chargement
# --------------------------------------------------------------------------
# Quart de coeur, 8x8 assemblages. L'indice [0,0] est au centre du coeur.
# '.' = hors coeur (reflecteur), 1/2/3 = zone d'enrichissement.

CHARGEMENT_UNIFORME = [
    "33333333",
    "33333333",
    "33333333",
    "33333333",
    "33333333",
    "333333..",
    "33333...",
    "333.....",
]

# Chargement "out-in" : le combustible le plus reactif en peripherie.
# C'est le levier classique d'aplatissement de la puissance.
CHARGEMENT_ZONE = [
    "11122233",
    "11122233",
    "11222333",
    "22223333",
    "22233333",
    "223333..",
    "33333...",
    "333.....",
]

# Positions de grappes (quart de coeur), en (ligne, colonne) assemblage.
GRAPPES = [(0, 0), (0, 4), (2, 2), (4, 0), (2, 6), (6, 2), (4, 4)]


# --------------------------------------------------------------------------
# Construction du maillage
# --------------------------------------------------------------------------
class Coeur:
    """Un coeur en quart de symetrie, pret a etre resolu.

    plan     : liste de chaines, voir CHARGEMENT_ZONE
    grappes  : liste de positions d'assemblages ou une grappe est inseree
    ppm_bore : concentration en bore soluble du moderateur
    n_par_ass: cellules de maillage par assemblage et par direction
    """

    def __init__(self, plan=CHARGEMENT_ZONE, grappes=(), ppm_bore=0.0,
                 n_par_ass=6, buckling_axial=BUCKLING_AXIAL):
        self.plan = [ligne for ligne in plan]
        self.grappes = set(map(tuple, grappes))
        self.ppm_bore = ppm_bore
        self.n_par_ass = n_par_ass
        self.buckling_axial = buckling_axial
        self.mats = bibliotheque()
        self._construire()

    # ---------------------------------------------------------------
    def _construire(self):
        n_ass = len(self.plan)
        n_refl_ass = int(round(EP_REFLECTEUR / PAS_ASSEMBLAGE))
        n_tot_ass = n_ass + n_refl_ass
        npa = self.n_par_ass
        h = PAS_ASSEMBLAGE / npa
        n = n_tot_ass * npa

        # Table des materiaux effectivement utilises
        materiaux, index = [], {}

        def idx(cle, base, grappe):
            signature = (cle, grappe)
            if signature not in index:
                m = avec_bore(self.mats[base], self.ppm_bore)
                if grappe:
                    m = avec_grappe(m)
                m.buckling_b2 = self.buckling_axial
                index[signature] = len(materiaux)
                materiaux.append(m)
            return index[signature]

        carte = np.zeros((n, n), int)
        self.carte_zone = np.full((n_tot_ass, n_tot_ass), ".", dtype="<U1")

        for ja in range(n_tot_ass):
            for ia in range(n_tot_ass):
                dans_plan = (ja < n_ass and ia < n_ass
                             and self.plan[ja][ia] != ".")
                if dans_plan:
                    z = self.plan[ja][ia]
                    base = f"zone{z}"
                    grappe = (ja, ia) in self.grappes
                    cle = f"z{z}"
                    self.carte_zone[ja, ia] = "G" if grappe else z
                else:
                    base, grappe, cle = "reflecteur", False, "refl"
                m = idx(cle, base, grappe)
                carte[ja * npa:(ja + 1) * npa, ia * npa:(ia + 1) * npa] = m

        self.materiaux = materiaux
        self.n_ass, self.n_tot_ass = n_ass, n_tot_ass
        self.mesh = Mesh(
            hx=np.full(n, h), hy=np.full(n, h),
            mat_map=carte, materials=materiaux,
            # Symetrie quart de coeur : reflexion au centre, vide a l'exterieur
            bc=dict(xmin="reflective", ymin="reflective",
                    xmax="vacuum", ymax="vacuum"),
        )

    # ---------------------------------------------------------------
    def resoudre(self, **kw):
        self.sol = solve(self.mesh, **kw)
        return self.sol

    # ---------------------------------------------------------------
    def puissance_assemblages(self):
        """Puissance moyenne par assemblage, normalisee a 1 sur le coeur.

        Puissance = somme_g kappa*Sigma_f,g * phi_g. Attention : ce n'est
        PAS le flux. Le flux thermique culmine dans le reflecteur, ou la
        puissance est rigoureusement nulle.
        """
        npa, npa_tot = self.n_par_ass, self.n_tot_ass
        flux = self.sol.flux                       # (G, ny, nx)
        V = self.mesh.volume

        dens = np.zeros(V.shape)
        for im, m in enumerate(self.materiaux):
            ksf = getattr(m, "kappa_sigma_f", np.zeros(2))
            if not np.any(ksf):
                continue
            masque = (self.mesh.mat_map == im)
            for g in range(flux.shape[0]):
                dens += masque * ksf[g] * flux[g]

        p_cellule = dens * V
        p_ass = p_cellule.reshape(npa_tot, npa, npa_tot, npa).sum(axis=(1, 3))

        actifs = p_ass > 0
        if actifs.any():
            p_ass = np.where(actifs, p_ass / p_ass[actifs].mean(), 0.0)
        return p_ass

    def facteur_point_chaud(self):
        """F_xy : puissance de l'assemblage le plus chaud / puissance moyenne.

        Grandeur de surete de premier ordre : elle fixe la marge a
        l'ebullition. Un coeur avec F_xy trop eleve doit etre deratine,
        meme si sa reactivite est parfaite.
        """
        p = self.puissance_assemblages()
        return float(p.max())


# --------------------------------------------------------------------------
# Reactivite
# --------------------------------------------------------------------------
def reactivite_pcm(k):
    """rho = (k - 1) / k, exprimee en pcm (1 pcm = 1e-5).

    On raisonne en reactivite plutot qu'en k parce que les effets
    s'additionnent (presque) lineairement en rho : le poids d'une grappe
    plus celui du bore donnent a peu pres le poids des deux ensemble.
    """
    return (k - 1.0) / k * 1e5


def poids_pcm(k_sans, k_avec):
    """Poids en reactivite d'un dispositif (grappe, bore, ...), en pcm.

    Par convention un absorbant a un poids positif : c'est la reactivite
    qu'il retire au coeur.
    """
    return reactivite_pcm(k_sans) - reactivite_pcm(k_avec)


def bore_critique(plan=CHARGEMENT_ZONE, grappes=(), n_par_ass=6,
                  ppm_min=0.0, ppm_max=4000.0, tol=0.5, verbose=False):
    """Cherche la concentration en bore qui rend le coeur critique (k = 1).

    C'est le calcul de base de l'exploitation : en debut de cycle le
    combustible est trop reactif, on compense avec du bore, puis on dilue
    progressivement a mesure que le combustible s'use.

    Methode de la secante sur rho(ppm), qui est quasi lineaire.
    """
    def k_de(ppm):
        c = Coeur(plan=plan, grappes=grappes, ppm_bore=ppm, n_par_ass=n_par_ass)
        return c.resoudre().k_eff

    a, b = ppm_min, ppm_max
    fa, fb = k_de(a) - 1.0, k_de(b) - 1.0
    if fa < 0:
        raise ValueError(f"coeur deja sous-critique sans bore (k={fa + 1:.5f})")
    if fb > 0:
        raise ValueError(f"toujours critique a {b:.0f} ppm (k={fb + 1:.5f})")

    for _ in range(40):
        c = b - fb * (b - a) / (fb - fa)
        fc = k_de(c) - 1.0
        if verbose:
            print(f"    {c:8.1f} ppm -> k = {fc + 1:.6f}")
        if fc > 0:
            a, fa = c, fc
        else:
            b, fb = c, fc
        if abs(b - a) < tol:
            break
    return 0.5 * (a + b)


# --------------------------------------------------------------------------
def afficher_carte(valeurs, plan, titre, fmt="{:6.3f}"):
    """Affiche une carte par assemblage, en quart de coeur."""
    n = len(plan)
    print(f"\n  {titre}")
    print("       " + "".join(f"{i:>7d}" for i in range(n)))
    for j in range(n):
        ligne = f"   {j:2d} |"
        for i in range(n):
            if plan[j][i] == ".":
                ligne += "      ."
            else:
                ligne += fmt.format(valeurs[j, i])
        print(ligne)
