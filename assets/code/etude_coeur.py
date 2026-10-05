"""
etude_coeur.py -- Etude de conception d'un coeur REP.

Six questions, dans l'ordre ou un concepteur de coeur se les pose :

  1. De quoi dispose-t-on ?          k_inf des combustibles
  2. Que coutent les fuites ?        k_inf -> k_eff, radiales et axiales
  3. Comment aplatir la puissance ?  plan de chargement et F_xy
  4. Comment rendre critique ?       recherche du bore critique
  5. Comment piloter ?               poids differentiel du bore
  6. Comment arreter ?               efficacite des grappes, et leur cout

Lancer :  python etude_coeur.py
"""

import numpy as np
from coeur import (Coeur, CHARGEMENT_UNIFORME, CHARGEMENT_ZONE, GRAPPES,
                   BUCKLING_AXIAL, bibliotheque, k_infini, reactivite_pcm,
                   poids_pcm, bore_critique, afficher_carte,
                   SIGMA_A_PAR_PPM_BORE, SIGMA_A_GRAPPE)


def titre(n, texte):
    print()
    print("=" * 74)
    print(f"{n}. {texte}")
    print("=" * 74)


# ==========================================================================
def etape1_combustibles():
    titre(1, "DE QUOI DISPOSE-T-ON ? k_infini des combustibles")
    print("""
  k_inf est le facteur de multiplication d'un milieu INFINI : aucune fuite.
  C'est le potentiel intrinseque du combustible, le plafond que le coeur
  reel ne pourra jamais atteindre.

      k_inf = (nuSf1 + nuSf2 * Ss12/Sa2) / (Sa1 + Ss12)

  Le numerateur est la production : fissions rapides, plus les neutrons qui
  ralentissent (Ss12/Sa2 = probabilite de survivre jusqu'a la fission
  thermique). Le denominateur est la disparition du groupe rapide.
""")
    mats = bibliotheque()
    print(f"  {'zone':<22} {'Sa2':>8} {'nuSf2':>8} {'k_inf':>9}")
    for cle in ("zone1", "zone2", "zone3"):
        m = mats[cle]
        print(f"  {m.name:<22} {m.sigma_a[1]:>8.4f} {m.nu_sigma_f[1]:>8.4f} "
              f"{k_infini(m):>9.5f}")
    print("""
  L'enrichissement joue surtout sur nuSf2. Il augmente aussi un peu Sa2
  (l'U235 absorbe sans fissionner une fois sur six environ), c'est pourquoi
  k_inf ne croit pas proportionnellement a l'enrichissement.
""")


# ==========================================================================
def etape2_fuites():
    titre(2, "QUE COUTENT LES FUITES ? de k_inf a k_eff")
    print("""
  Un coeur reel perd des neutrons par ses frontieres :

      k_eff = k_inf * P_NL        (P_NL = probabilite de non-fuite)

  On separe les deux directions. La fuite radiale est deja dans le calcul
  2D. La fuite axiale n'y est pas : le modele est plan. On la reintroduit
  par un buckling B_z^2 = (pi/H_ex)^2, qui ajoute une pseudo-absorption
  D*B_z^2. Si le flux axial est un cosinus, c'est exact.
""")
    mats = bibliotheque()
    kinf = k_infini(mats["zone3"])

    c_sans = Coeur(plan=CHARGEMENT_UNIFORME, buckling_axial=0.0)
    k_rad = c_sans.resoudre().k_eff
    c_avec = Coeur(plan=CHARGEMENT_UNIFORME, buckling_axial=BUCKLING_AXIAL)
    k_tot = c_avec.resoudre().k_eff

    print(f"  B_z^2 = {BUCKLING_AXIAL:.4e} 1/cm2   (hauteur active 366 cm)")
    print()
    print(f"  {'k_inf (zone3, milieu infini)':<42} {kinf:>10.5f}")
    print(f"  {'k_eff, fuites radiales seules':<42} {k_rad:>10.5f}   "
          f"P_NL,rad = {k_rad / kinf:.5f}")
    print(f"  {'k_eff, fuites radiales + axiales':<42} {k_tot:>10.5f}   "
          f"P_NL,tot = {k_tot / kinf:.5f}")
    print()
    print(f"  cout des fuites radiales : {reactivite_pcm(kinf) - reactivite_pcm(k_rad):7.0f} pcm")
    print(f"  cout des fuites axiales  : {poids_pcm(k_rad, k_tot):7.0f} pcm")
    print("""
  Les fuites radiales sont faibles parce que le reflecteur en eau renvoie
  une bonne partie des neutrons dans le coeur : c'est tout son role, et
  c'est pour cela qu'il est economiquement rentable.
""")


# ==========================================================================
def etape3_aplatissement():
    titre(3, "COMMENT APLATIR LA PUISSANCE ? le plan de chargement")
    print("""
  Le facteur de point chaud F_xy = puissance max / puissance moyenne est
  une contrainte de SURETE, pas de performance. Il fixe la marge avant
  ebullition (crise d'ebullition, DNB) sur le crayon le plus chaud. Un
  coeur avec un F_xy trop eleve doit etre bride en puissance, meme si sa
  reactivite est excellente.

  Le flux naturel d'un coeur homogene est en Bessel J0 : tres pique au
  centre. La parade est le chargement zone, dit "out-in" : on met le
  combustible le plus reactif en PERIPHERIE, la ou le flux serait faible.
""")
    resultats = {}
    for plan, nom in ((CHARGEMENT_UNIFORME, "uniforme (tout zone3)"),
                      (CHARGEMENT_ZONE, "zone out-in (1/2/3)")):
        c = Coeur(plan=plan)
        c.resoudre()
        p = c.puissance_assemblages()
        resultats[nom] = (c, p)
        print(f"  {nom:<24} k_eff = {c.sol.k_eff:.5f}   "
              f"F_xy = {c.facteur_point_chaud():.4f}")

    for nom, (c, p) in resultats.items():
        afficher_carte(p[:c.n_ass, :c.n_ass], c.plan,
                       f"Puissance par assemblage -- {nom}  (moyenne = 1)")

    c_u = resultats["uniforme (tout zone3)"][0]
    c_z = resultats["zone out-in (1/2/3)"][0]
    print(f"""
  F_xy passe de {c_u.facteur_point_chaud():.3f} a {c_z.facteur_point_chaud():.3f}.
  Le prix a payer est {poids_pcm(c_u.sol.k_eff, c_z.sol.k_eff):.0f} pcm de reactivite en moins,
  puisqu'on a remplace du combustible riche par du combustible pauvre.

  C'est LE compromis central de la conception de coeur : reactivite contre
  aplatissement. On achete de la marge thermique avec de la reactivite.
""")
    return resultats


# ==========================================================================
def etape4_bore_critique():
    titre(4, "COMMENT RENDRE CRITIQUE ? le bore soluble")
    print("""
  Un coeur neuf est tres sur-critique : il faut de la reactivite en reserve
  pour tenir tout le cycle, puisque le combustible s'use. Il faut donc
  compenser cet exces au demarrage, sans deformer la puissance.

  Le bore dilue dans l'eau du circuit primaire est ideal pour cela : il est
  reparti uniformement, donc il ne cree aucun point chaud. C'est un
  absorbant quasi purement thermique (loi en 1/v), il n'agit donc que sur
  le groupe 2.
""")
    c0 = Coeur(plan=CHARGEMENT_ZONE)
    c0.resoudre()
    rho0 = reactivite_pcm(c0.sol.k_eff)
    print(f"  reactivite sans bore : {rho0:.0f} pcm  (k = {c0.sol.k_eff:.5f})")
    print(f"  Sigma_a2 ajoutee     : {SIGMA_A_PAR_PPM_BORE:.2e} 1/cm par ppm")
    print("\n  recherche du bore critique (methode de la secante) :")
    ppm = bore_critique(plan=CHARGEMENT_ZONE, verbose=True)
    print(f"\n  --> bore critique = {ppm:.0f} ppm")
    print(f"""
  A COMPARER AU REEL : un REP demarre son cycle vers 1200-1800 ppm, pas
  {ppm:.0f}. L'ecart n'est pas une erreur du modele, il est PHYSIQUE : notre
  coeur est entierement neuf, alors qu'un coeur reel n'est recharge qu'au
  tiers ou au quart, le reste etant du combustible deja use, donc moins
  reactif.

  Et c'est precisement pour cela que les POISONS CONSOMMABLES existent
  (gadolinium, crayons boures). Un bore critique aussi eleve serait
  inacceptable : au-dela de ~1800 ppm le coefficient de temperature du
  moderateur devient POSITIF, ce qui est interdit. Voir etape 5.
""")
    return ppm


# ==========================================================================
def etape5_poids_bore(ppm_crit):
    titre(5, "COMMENT PILOTER ? poids differentiel du bore")
    print("""
  Le poids differentiel drho/dppm est ce qui gouverne le pilotage : c'est
  la reactivite gagnee en diluant d'un ppm. On dilue progressivement le
  bore au long du cycle pour compenser l'usure du combustible.
""")
    c0 = Coeur(plan=CHARGEMENT_ZONE)
    k0 = c0.resoudre().k_eff
    print(f"  {'ppm':>7} {'k_eff':>10} {'rho [pcm]':>11} {'poids cumule':>13} "
          f"{'differentiel':>14}")
    liste, precedent = [], None
    for ppm in (0, 500, 1000, 1500, 2000, 2500):
        c = Coeur(plan=CHARGEMENT_ZONE, ppm_bore=ppm)
        k = c.resoudre().k_eff
        rho = reactivite_pcm(k)
        cumul = poids_pcm(k0, k)
        diff = ""
        if precedent is not None:
            d = (precedent[1] - rho) / (ppm - precedent[0])
            diff = f"{-d:>10.2f} pcm/ppm"
            liste.append(-d)
        print(f"  {ppm:>7} {k:>10.5f} {rho:>11.0f} {cumul:>13.0f} {diff:>14}")
        precedent = (ppm, rho)

    print(f"""
  Le poids differentiel vaut environ {np.mean(liste):.2f} pcm/ppm et il DECROIT
  en valeur absolue quand la concentration augmente ({liste[0]:.2f} -> {liste[-1]:.2f}).

  Pourquoi ? Auto-protection. Chaque atome de bore ajoute absorbe des
  neutrons thermiques, ce qui deprime le flux thermique. Les atomes
  suivants voient donc moins de neutrons a absorber et rapportent moins.
  Tout absorbant sature ainsi ; c'est un comportement general, pas une
  particularite du bore.

  CONSEQUENCE DE SURETE. Le bore est dans l'EAU. Si la temperature monte,
  l'eau se dilate, il y a moins d'eau dans le coeur -- donc moins de bore,
  donc un gain de reactivite. Cet effet s'oppose au coefficient de
  temperature moderateur naturellement negatif. Trop de bore et la somme
  devient POSITIVE : une hausse de temperature augmenterait la puissance,
  qui augmenterait la temperature. Emballement. C'est la vraie limite
  haute sur la concentration en bore, et la raison d'etre des poisons
  consommables.
""")


# ==========================================================================
def etape6_grappes(ppm_crit):
    titre(6, "COMMENT ARRETER ? efficacite des grappes de commande")
    print("""
  Le bore est trop lent pour piloter (il faut diluer ou borer, des dizaines
  de minutes) et incapable d'arreter le reacteur rapidement. Les grappes de
  commande, elles, tombent en quelques secondes.

  Contrairement au bore, elles sont LOCALES : elles deforment la puissance.
  On paie leur efficacite en point chaud.
""")
    ppm = round(ppm_crit)
    c_ext = Coeur(plan=CHARGEMENT_ZONE, ppm_bore=ppm)
    k_ext = c_ext.resoudre().k_eff
    c_ins = Coeur(plan=CHARGEMENT_ZONE, ppm_bore=ppm, grappes=GRAPPES)
    k_ins = c_ins.resoudre().k_eff

    print(f"  Sigma_a2 ajoutee par grappe : {SIGMA_A_GRAPPE:.4f} 1/cm")
    print(f"  bore fixe a {ppm} ppm (coeur critique grappes extraites)")
    print()
    print(f"  {'grappes extraites':<26} k = {k_ext:.5f}   rho = {reactivite_pcm(k_ext):>7.0f} pcm"
          f"   F_xy = {c_ext.facteur_point_chaud():.4f}")
    print(f"  {'grappes inserees':<26} k = {k_ins:.5f}   rho = {reactivite_pcm(k_ins):>7.0f} pcm"
          f"   F_xy = {c_ins.facteur_point_chaud():.4f}")
    print(f"\n  --> efficacite totale = {poids_pcm(k_ext, k_ins):.0f} pcm")

    # efficacite grappe par grappe, pour montrer l'effet de position
    print("\n  Efficacite individuelle, selon la position :")
    print(f"  {'position':>12} {'k_eff':>10} {'efficacite':>13}")
    individuelles = []
    for g in GRAPPES:
        c = Coeur(plan=CHARGEMENT_ZONE, ppm_bore=ppm, grappes=[g])
        k = c.resoudre().k_eff
        eff = poids_pcm(k_ext, k)
        individuelles.append(eff)
        print(f"  {str(g):>12} {k:>10.5f} {eff:>10.0f} pcm")

    somme = sum(individuelles)
    total = poids_pcm(k_ext, k_ins)
    print(f"""
  PIEGE DE SYMETRIE. La position (0,0) est au CENTRE du coeur, la ou le
  flux est maximal -- et pourtant c'est la moins efficace ({individuelles[0]:.0f} pcm).
  Ce n'est pas une erreur. Le plan de symetrie est au BORD de l'assemblage
  (0,0), pas en son milieu : en coeur complet, cette position se reflete en
  un bloc serre de 2x2 grappes accolees au centre. Elles s'ombragent
  mutuellement au maximum. Toute position en quart de coeur represente
  quatre assemblages physiques, et leur agencement compte autant que leur
  nombre.

  somme des efficacites individuelles : {somme:.0f} pcm
  efficacite du groupe entier         : {total:.0f} pcm
  ecart                               : {total - somme:+.0f} pcm ({(total / somme - 1) * 100:+.1f} %)

  LES EFFICACITES NE S'ADDITIONNENT PAS. Ici le groupe vaut PLUS que la
  somme de ses parties, ce qui surprend. L'experience ci-dessous separe les
  deux mecanismes en jeu.
""")

    # ---- experience : ombrage contre anti-ombrage --------------------
    print("  Experience : deux grappes, accolees puis eloignees.")
    print(f"  {'cas':<14} {'seule A':>9} {'seule B':>9} {'somme':>8} "
          f"{'ensemble':>10} {'ecart':>9}")

    def poids(gs):
        c = Coeur(plan=CHARGEMENT_ZONE, ppm_bore=ppm, grappes=gs)
        return poids_pcm(k_ext, c.resoudre().k_eff)

    for nom, a, b in (("accolees", (3, 3), (3, 4)),
                      ("eloignees", (1, 1), (5, 4)),
                      ("eloignees 2", (0, 6), (6, 0))):
        wa, wb, wab = poids([a]), poids([b]), poids([a, b])
        print(f"  {nom:<14} {wa:>9.0f} {wb:>9.0f} {wa + wb:>8.0f} "
              f"{wab:>10.0f} {100 * (wab / (wa + wb) - 1):>+8.1f} %")

    print(f"""
  Deux effets opposes, et c'est la DISTANCE qui decide lequel gagne :

  * ACCOLEES : -14 %. C'est l'ombrage (shadowing) classique. Chaque grappe
    creuse le flux thermique autour d'elle ; sa voisine baigne dans ce
    creux, voit moins de neutrons, et en absorbe donc moins. Deux grappes
    collees se genent.

  * ELOIGNEES : +7 a +18 %. C'est l'anti-ombrage, et le mecanisme est
    different. Une grappe SEULE fait basculer la nappe de flux a l'oppose
    d'elle-meme : le coeur "fuit" lateralement la zone empoisonnee, ce qui
    limite l'efficacite de la grappe. Avec des grappes reparties
    symetriquement, le flux n'a plus nulle part ou basculer. Chacune
    conserve alors la totalite de son efficacite au lieu de se la faire
    rogner par sa propre deformation de flux.

  Le groupe complet de {len(GRAPPES)} grappes est reparti, donc l'anti-ombrage l'emporte
  et le total depasse la somme de {(total / somme - 1) * 100:+.0f} %.

  CONSEQUENCE PRATIQUE. L'efficacite d'un groupe de grappes ne se deduit
  JAMAIS de mesures individuelles : ni par addition, ni par un facteur
  correctif universel, puisque le signe de l'ecart change avec la
  geometrie. Elle se calcule, ou se mesure groupe par groupe. C'est
  exactement ce que font les essais physiques au demarrage d'un coeur
  neuf : on mesure l'efficacite de chaque GROUPE de grappes tel qu'il sera
  reellement manoeuvre.
""")
    return c_ext, c_ins


# ==========================================================================
def graphiques(res_chargement, c_ext, c_ins):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  (graphiques ignores : {e})")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 11))
    cas = [
        (res_chargement["uniforme (tout zone3)"], "Chargement uniforme"),
        (res_chargement["zone out-in (1/2/3)"], "Chargement zone out-in"),
    ]
    for ax, ((c, p), nom) in zip(axes[0], cas):
        n = c.n_ass
        pm = np.ma.masked_where(p[:n, :n] == 0, p[:n, :n])
        im = ax.imshow(pm, origin="lower", cmap="inferno", vmin=0, vmax=2.2)
        ax.set_title(f"{nom}\nF_xy = {c.facteur_point_chaud():.3f}")
        ax.set_xlabel("assemblage i")
        ax.set_ylabel("assemblage j")
        for j in range(n):
            for i in range(n):
                if p[j, i] > 0:
                    ax.text(i, j, f"{p[j, i]:.2f}", ha="center", va="center",
                            fontsize=6.5,
                            color="white" if p[j, i] < 1.4 else "black")
        fig.colorbar(im, ax=ax, label="puissance relative")

    for ax, c, nom in ((axes[1][0], c_ext, "Grappes extraites"),
                       (axes[1][1], c_ins, "Grappes inserees")):
        p = c.puissance_assemblages()
        n = c.n_ass
        pm = np.ma.masked_where(p[:n, :n] == 0, p[:n, :n])
        im = ax.imshow(pm, origin="lower", cmap="inferno", vmin=0, vmax=2.2)
        ax.set_title(f"{nom}\nk = {c.sol.k_eff:.5f}, F_xy = {c.facteur_point_chaud():.3f}")
        ax.set_xlabel("assemblage i")
        ax.set_ylabel("assemblage j")
        for (gj, gi) in GRAPPES:
            if nom.startswith("Grappes ins"):
                ax.plot(gi, gj, "c+", ms=13, mew=2)
        fig.colorbar(im, ax=ax, label="puissance relative")

    fig.suptitle("Distribution de puissance par assemblage (quart de coeur, "
                 "centre en bas a gauche)", fontsize=13)
    fig.tight_layout()
    fig.savefig("cartes_puissance.png", dpi=130)
    print("  --> cartes_puissance.png ecrit")


# ==========================================================================
def main():
    etape1_combustibles()
    etape2_fuites()
    res = etape3_aplatissement()
    ppm = etape4_bore_critique()
    etape5_poids_bore(ppm)
    c_ext, c_ins = etape6_grappes(ppm)
    print("=" * 74)
    graphiques(res, c_ext, c_ins)


if __name__ == "__main__":
    main()
