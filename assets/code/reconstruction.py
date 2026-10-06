"""
reconstruction.py -- Des parametres de resonances aux sections efficaces.

Un fichier evalue ENDF-6 ne contient pas, dans la region des resonances
resolues, la section efficace elle-meme : il en donne les PARAMETRES
(energie E_r, spin J, largeurs Gamma_n, Gamma_gamma, Gamma_f) et un fond
tabule (MF=3). Il faut reconstruire la courbe, puis l'elargir par l'agitation
thermique des noyaux cibles. C'est le travail des modules RECONR et BROADR
du code NJOY ; on le refait ici, en Python, pour le formalisme utilise par
ENDF/B-VIII.0 pour l'U-235 et l'U-238 : Reich-Moore (LRF = 3).

LE FORMALISME DE REICH-MOORE (C.W. Reich, M.S. Moore, Phys. Rev. 111 (1958)
929 ; ENDF-102, annexe D)
--------------------------------------------------------------------------
Theorie de la matrice R ou les tres nombreuses voies de capture gamma sont
"eliminees" : elles n'apparaissent que par une largeur Gamma_gamma ajoutee
au denominateur. Pour chaque onde l et chaque spin J, on construit la
matrice (voie neutron + deux voies de fission)

    K_cc'(E) = sum_r  b_rc b_rc' / (E_r - E - i Gamma_gamma,r / 2)
    b_rn = sqrt(Gamma_n,r(E) / 2),  b_rf = signe(G) sqrt(|G_f,r| / 2)

et la matrice de collision U = e^{-2 i phi} [ 2 (I - i K)^{-1} - I ] (voie n).
D'ou, avec X = (I - i K)^{-1} et g_J = (2J + 1) / (2 (2I + 1)) :

    sigma_t = (2 pi / k^2) g_J (1 - Re U_nn)
    sigma_e = (pi / k^2)   g_J |1 - U_nn|^2
    sigma_f = (4 pi / k^2) g_J (|X_nfA|^2 + |X_nfB|^2)
    sigma_gamma = sigma_t - sigma_e - sigma_f

La largeur neutron depend de l'energie par le facteur de penetration
P_l(k a) : Gamma_n,r(E) = Gamma_n,r * P_l(E) / P_l(|E_r|). Le dephasage de
sphere dure phi_l(k a') donne la diffusion potentielle.

L'EFFET DOPPLER (noyau cible en gaz libre a la temperature T)
-------------------------------------------------------------
Avec x = sqrt(A E' / kT) et y = sqrt(A E / kT) :

    sigma(y, T) = 1 / (sqrt(pi) y^2) int_0^inf x^2 sigma(x, 0) [e^{-(x-y)^2} - e^{-(x+y)^2}] dx

(c'est la formule exacte du noyau "SIGMA1" de Cullen et Weisbin, 1976.)

Unites : energies en eV, longueurs en 10^-12 cm (convention ENDF), sections
efficaces en barns (pi/k^2 est alors directement en barns).
"""

import numpy as np

KB = 8.617333262e-5          # constante de Boltzmann, eV/K
# k = C_K * AWRI/(AWRI+1) * sqrt(E) en (10^-12 cm)^-1 : sqrt(2 m_n c^2) / (hbar c)
C_K = np.sqrt(2 * 939.56542052e6) / 1.973269804e-5 * 1e-12


# ==========================================================================
# Facteurs de penetration et dephasages (sphere dure)
# ==========================================================================
def penetration(l, rho):
    if l == 0:
        return rho
    if l == 1:
        return rho ** 3 / (1 + rho ** 2)
    if l == 2:
        return rho ** 5 / (9 + 3 * rho ** 2 + rho ** 4)
    raise ValueError("l > 2 non prevu")


def dephasage(l, rho):
    if l == 0:
        return rho
    if l == 1:
        return rho - np.arctan(rho)
    if l == 2:
        return rho - np.arctan(3 * rho / (3 - rho ** 2))
    raise ValueError("l > 2 non prevu")


# ==========================================================================
# Reich-Moore
# ==========================================================================
class ReichMoore:
    """Region des resonances resolues d'un fichier ENDF-6 (LRF = 3)."""

    def __init__(self, rrange, awr):
        self.EL, self.EH = rrange["EL"], rrange["EH"]
        self.I = rrange["SPI"]
        self.AP = rrange["AP"]
        self.NAPS = rrange["NAPS"]
        self.awr = awr
        self.sections = rrange["sections"]

    def _rayons(self, sec):
        awri = sec["AWRI"]
        if self.NAPS == 0:
            a_canal = 0.123 * awri ** (1 / 3) + 0.08
        else:
            a_canal = self.AP
        a_diff = sec["APL"] if sec["APL"] != 0 else self.AP
        return a_canal, a_diff

    def _spins(self, l):
        """(J, g_J, nombre de voies (s, J)) pour l'onde l."""
        I = self.I
        voies = {}
        for s in sorted({abs(I - 0.5), I + 0.5}):   # pour I = 0, une seule valeur
            J = abs(l - s)
            while J <= l + s + 1e-9:
                voies[round(J, 1)] = voies.get(round(J, 1), 0) + 1
                J += 1
        return [(J, (2 * J + 1) / (2 * (2 * I + 1)), n) for J, n in sorted(voies.items())]

    def sections_efficaces(self, E, morceau=1500):
        """(total, elastique, fission, capture) en barns, a 0 K."""
        E = np.asarray(E, float)
        out = np.zeros((4, len(E)))
        for i0 in range(0, len(E), morceau):
            e = E[i0:i0 + morceau]
            out[:, i0:i0 + morceau] = self._morceau(e)
        return out

    def _morceau(self, E):
        st = np.zeros_like(E)
        se = np.zeros_like(E)
        sf = np.zeros_like(E)
        for sec in self.sections:
            l = int(sec["L"])
            awri = sec["AWRI"]
            a_canal, a_diff = self._rayons(sec)
            k = C_K * awri / (awri + 1) * np.sqrt(E)
            P = penetration(l, k * a_canal)
            phi = dephasage(l, k * a_diff)
            pik2 = np.pi / k ** 2
            ER = np.asarray(sec["ER"])
            AJ = np.abs(np.asarray(sec["AJ"]))
            for J, g, nvoies in self._spins(l):
                sel = np.isclose(AJ, J)
                for _ in range(nvoies - (1 if sel.any() else 0)):
                    # voie sans resonance : diffusion potentielle seule
                    st += 4 * pik2 * g * np.sin(phi) ** 2
                    se += 4 * pik2 * g * np.sin(phi) ** 2
                if not sel.any():
                    continue
                Er = ER[sel]
                Gn = np.asarray(sec["GN"])[sel]
                Gg = np.asarray(sec["GG"])[sel]
                Gfa = np.asarray(sec["GFA"])[sel]
                Gfb = np.asarray(sec["GFB"])[sel]
                kr = C_K * awri / (awri + 1) * np.sqrt(np.abs(Er))
                Pr = penetration(l, kr * a_canal)
                Gn_E = Gn[None, :] * P[:, None] / Pr[None, :]
                inv = 1.0 / (Er[None, :] - E[:, None] - 0.5j * Gg[None, :])
                bn = np.sqrt(Gn_E / 2)
                ba = np.sign(Gfa) * np.sqrt(np.abs(Gfa) / 2)
                bb = np.sign(Gfb) * np.sqrt(np.abs(Gfb) / 2)
                K = np.empty((len(E), 3, 3), complex)
                K[:, 0, 0] = (bn * bn * inv).sum(1)
                K[:, 0, 1] = K[:, 1, 0] = (bn * ba * inv).sum(1)
                K[:, 0, 2] = K[:, 2, 0] = (bn * bb * inv).sum(1)
                K[:, 1, 1] = (ba * ba * inv).sum(1)
                K[:, 1, 2] = K[:, 2, 1] = (ba * bb * inv).sum(1)
                K[:, 2, 2] = (bb * bb * inv).sum(1)
                M = np.eye(3)[None, :, :] - 1j * K
                if np.all(Gfa == 0) and np.all(Gfb == 0):
                    Xnn = 1.0 / M[:, 0, 0]
                    Xna = Xnb = np.zeros_like(Xnn)
                else:
                    X = np.linalg.inv(M)
                    Xnn, Xna, Xnb = X[:, 0, 0], X[:, 0, 1], X[:, 0, 2]
                U = np.exp(-2j * phi) * (2 * Xnn - 1)
                st += 2 * pik2 * g * (1 - U.real)
                se += pik2 * g * np.abs(1 - U) ** 2
                sf += 4 * pik2 * g * (np.abs(Xna) ** 2 + np.abs(Xnb) ** 2)
        return np.array([st, se, sf, st - se - sf])

    def grille(self, fine=True, par_decade=300):
        """Grille d'energie : logarithmique, plus des points serres autour de
        chaque resonance, en fractions de sa largeur totale Gamma. La grille
        fine (~120 points par resonance) sert aux integrales et a l'effet
        Doppler ; la grille legere (~40 points) a l'affichage."""
        if fine:
            ecarts = np.concatenate([np.linspace(0, 3, 31), np.geomspace(3.3, 300, 30)])
        else:
            ecarts = np.array([0, .1, .2, .35, .5, .7, 1, 1.4, 2, 3, 5, 8, 14, 25, 50, 100])
        E = [np.logspace(np.log10(self.EL), np.log10(self.EH),
                         int(par_decade * np.log10(self.EH / self.EL)) + 1)]
        for sec in self.sections:
            Er = np.asarray(sec["ER"])
            G = (np.asarray(sec["GN"]) + np.asarray(sec["GG"]) + np.abs(sec["GFA"])
                 + np.abs(sec["GFB"]))
            ok = (Er > self.EL) & (Er < self.EH)
            for sg in (-1, 1):
                E.append((Er[ok][:, None] + sg * G[ok][:, None] * ecarts[None, :]).ravel())
        E = np.unique(np.concatenate(E))
        return E[(E >= self.EL) & (E <= self.EH)]


# ==========================================================================
# Elargissement Doppler
# ==========================================================================
def doppler(E0, sigma0, awr, T, E_cible=None, fenetre=5.0):
    """Section efficace a la temperature T a partir de la section a 0 K,
    tabulee sur la grille E0 (lineaire entre les points). Formule exacte du
    gaz libre (noyau SIGMA1), integree par trapezes en x = sqrt(A E / kT).
    sigma0 peut etre un tableau (n_reactions, len(E0)) : le noyau est alors
    calcule une seule fois pour toutes les reactions."""
    if E_cible is None:
        E_cible = E0
    sigma0 = np.atleast_2d(sigma0)
    a = awr / (KB * T)
    x = np.sqrt(a * E0)
    f = x[None, :] ** 2 * sigma0
    y = np.sqrt(a * np.asarray(E_cible))
    out = np.empty((sigma0.shape[0], len(y)))
    lo_all = np.searchsorted(x, np.maximum(y - fenetre, 0.0))
    hi_all = np.searchsorted(x, y + fenetre)
    for i, yi in enumerate(y):
        lo, hi = lo_all[i], hi_all[i]
        if hi - lo < 2:
            out[:, i] = [np.interp(yi, x, s) for s in sigma0]
            continue
        xs = x[lo:hi]
        noyau = np.exp(-(xs - yi) ** 2)
        if yi < fenetre:
            noyau = noyau - np.exp(-(xs + yi) ** 2)
        out[:, i] = np.trapezoid(f[:, lo:hi] * noyau[None, :], xs, axis=1) / (np.sqrt(np.pi) * yi * yi)
    return out if out.shape[0] > 1 else out[0]


# ==========================================================================
# Un noyau : parametres + fond tabule (MF=3)
# ==========================================================================
class Noyau:
    """Sections efficaces d'un fichier ENDF-6 : resonances resolues
    reconstruites (Reich-Moore) + fond MF=3 ; au-dessus, MF=3 tel quel (les
    fichiers traites ici ont LSSF = 1 : MF=3 contient deja les sections
    moyennes de la region non resolue)."""

    MT = {"total": 1, "elastique": 2, "fission": 18, "capture": 102}

    def __init__(self, chemin):
        import endf
        mat = endf.Material(chemin)
        mf2 = mat.section_data[2, 151]
        self.awr = mf2["AWR"]
        rr = mf2["isotopes"][0]["ranges"]
        resolu = [r for r in rr if r["LRU"] == 1][0]
        assert resolu["LRF"] == 3, "seul Reich-Moore est implemente"
        non_resolu = [r for r in rr if r["LRU"] == 2]
        assert all(r["LSSF"] == 1 for r in non_resolu), "LSSF = 0 non traite"
        self.rm = ReichMoore(resolu, self.awr)
        self.mf3 = {nom: mat.section_data[3, mt]["sigma"] for nom, mt in self.MT.items()}

    def zero_kelvin(self, fine=True):
        """(E, {reaction: sigma}) a 0 K sur une grille adaptee aux resonances."""
        E_res = self.rm.grille(fine=fine)
        res = self.rm.sections_efficaces(E_res)
        fond = {nom: f(E_res) for nom, f in self.mf3.items()}
        sig = {"total": res[0] + fond["total"], "elastique": res[1] + fond["elastique"],
               "fission": res[2] + fond["fission"], "capture": res[3] + fond["capture"]}
        # au-dessus de la region resolue : MF=3 tel quel
        E_haut = np.unique(np.concatenate([self.mf3[n].x for n in self.mf3]))
        E_haut = E_haut[E_haut > self.rm.EH]
        E_haut = np.unique(np.concatenate([E_haut, np.logspace(np.log10(self.rm.EH),
                                                               np.log10(E_haut.max()), 400)[1:]]))
        E = np.concatenate([E_res, E_haut])
        for nom in sig:
            sig[nom] = np.concatenate([sig[nom], self.mf3[nom](E_haut)])
        return E, sig

    def temperature(self, E, sig0, T, reactions=("total", "elastique", "fission", "capture"),
                    E_sortie=None, emin=None, emax=None):
        """Sections a la temperature T, sur la grille E_sortie (par defaut E).
        L'integrale de Doppler utilise la grille d'entree E (fine) ; seule la
        region resolue est elargie, le reste (lisse) est interpole.
        emin/emax : restreindre le calcul a une fenetre (zoom)."""
        if E_sortie is None:
            E_sortie = E
        lo = self.rm.EL if emin is None else emin
        hi = self.rm.EH * 1.0000001 if emax is None else min(emax, self.rm.EH * 1.0000001)
        m = (E_sortie >= lo) & (E_sortie <= hi)
        # le noyau a besoin des points voisins : on garde une marge en entree
        marge = (E >= lo / 1.5) & (E <= hi * 1.5) & (E <= self.rm.EH * 1.0000001)
        pile = np.array([sig0[nom][marge] for nom in reactions])
        br = np.atleast_2d(doppler(E[marge], pile, self.awr, T, E_cible=E_sortie[m]))
        out = {}
        for j, nom in enumerate(reactions):
            s = np.interp(E_sortie, E, sig0[nom])
            s[m] = br[j]
            out[nom] = s
        return out
