---
title: "Neutronique du cœur : comprendre et concevoir"
date: 2026-10-02 09:00:00 +0200
categories: ["III · Du noyau au réacteur", "Conception de cœur"]
tags: [REP, bore, grappes, facteur de point chaud]
description: "Plan de chargement, bore critique, efficacité des grappes et anti-ombrage, sur un cœur REP modélisé."
image:
  path: /assets/img/nucleaire/cartes_puissance.png
  alt: "Neutronique du cœur : comprendre et concevoir"
lecon: 8
partie: "III"
objectifs:
  - 'Passer de \(k_\infty\) à \(k_{\text{eff}}\) : le rôle des fuites et du réflecteur.'
  - "Aplatir la puissance avec un plan de chargement, et en mesurer le prix en réactivité."
  - "Comparer le bore soluble et les grappes de commande ; comprendre ombrage et anti-ombrage."
  - "Savoir ce qu'un calcul de neutronique statique ne contient pas."
prerequis: [7]
code: [coeur.py, etude_coeur.py, diffusion.py]
sources: >-
  Tous les chiffres (k, pcm, F_xy, ppm) : `coeur.py` et `etude_coeur.py`, avec les constantes à deux groupes de `coeur.py`, dont l'origine est discutée à la leçon 7. Ordres de grandeur d'exploitation (xénon, β_eff) : valeurs usuelles des manuels cités.
bibliographie:
  - cle: bussac1985
    note: "La référence francophone."
  - cle: reuss2003
    note: "Plus accessible, excellent sur les deux groupes et le contrôle de réactivité."
  - cle: duderstadt1976
    note: "Le standard anglo-saxon."
  - cle: stacey2007
    note: "Méthodes nodales et facteurs de discontinuité."
  - cle: smith1986
    note: "L'homogénéisation et les facteurs de discontinuité."
  - cle: pusa2010
    note: "La méthode CRAM pour l'épuisement."
  - cle: romano2015
    note: "OpenMC, pour confronter ces résultats à un Monte-Carlo."
---

{% include cours-entete.html %}

Document de synthèse accompagnant [`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }}) et [`etude_coeur.py`]({{ '/assets/code/etude_coeur.py' | relative_url }}).
Tous les chiffres cités sont produits par le code, pas estimés.


![Cartes de puissance : chargement uniforme ou zoné, grappes extraites ou insérées]({{ '/assets/img/nucleaire/cartes_puissance.png' | relative_url }})
_Cartes de puissance : chargement uniforme ou zoné, grappes extraites ou insérées._

---

## 1. Le cycle de vie du neutron

Tout part d'une question comptable : sur 100 neutrons d'une génération,
combien en produit-on à la suivante ? Ce rapport, c'est $$k$$.

Le parcours d'un neutron dans un REP :

1. **Naissance** par fission, à environ 2 MeV (groupe rapide).
2. **Ralentissement** par chocs sur l'hydrogène de l'eau. Il faut ~20 chocs
   pour tomber à 0,025 eV. C'est le rôle du modérateur.
3. **Traversée des résonances** de l'U238 entre 6 eV et 200 eV, où il peut
   être capturé stérilement. C'est le principal risque du trajet.
4. **Thermalisation**, puis diffusion jusqu'à absorption.
5. **Absorption** : soit une fission (on repart en 1), soit une capture
   stérile (dans l'U238, le bore, les structures, une grappe).

Le modèle à deux groupes comprime ce parcours en trois nombres par matériau :
$$\Sigma_{1\to2}$$ (le ralentissement), $$\Sigma_{a2}$$ (l'absorption thermique),
$$\nu\Sigma_{f2}$$ (la production thermique). D'où :

$$
k_\infty = \frac{\nu\Sigma_{f1} + \nu\Sigma_{f2}\, \Sigma_{1\to2}/\Sigma_{a2}}{\Sigma_{a1} + \Sigma_{1\to2}}
$$

Le terme $$\Sigma_{1\to2}/\Sigma_{a2}$$ est la probabilité de survivre au ralentissement plutôt
que d'être absorbé en route.

**Valeurs du modèle :** zone1 = 1,21644 · zone2 = 1,25194 · zone3 = 1,29054.

## 2. De $$k_\infty$$ à $$k_{\text{eff}}$$ : le prix des fuites

$$k_\infty$$ suppose un milieu infini. Un cœur réel fuit :

$$
k_{\text{eff}} = k_\infty \times P_{\text{NL}}
$$

Mesuré par le code, pour le chargement uniforme :

| | k | P_NL |
|---|---|---|
| $$k_\infty$$ (zone3, milieu infini) | 1,29054 | — |
| fuites radiales seules | 1,27841 | 0,99060 |
| fuites radiales + axiales | 1,27317 | 0,98654 |

**735 pcm** de fuite radiale, **322 pcm** de fuite axiale. C'est peu, et
c'est tout le mérite du réflecteur : les deux assemblages d'eau autour du
cœur renvoient l'essentiel des neutrons échappés. Un cœur nu perdrait
plusieurs milliers de pcm.

Noter le procédé pour la fuite axiale : le modèle est plan, il ne voit pas
la troisième dimension. On la réintroduit par un **buckling** $$B_z^2 = (\pi/H_{\text{ex}})^2 \approx 7{,}3\cdot10^{-5}~\text{cm}^{-2}$$, qui ajoute une pseudo-absorption $$D B_z^2$$. Si le
flux axial est un cosinus — ce qui est le cas pour un cœur axialement nu —
c'est exact, pas approché.

## 3. Aplatir la puissance : le vrai métier

Le flux naturel d'un cœur homogène est en Bessel J₀ : très piqué au centre.
Or ce qui limite un réacteur n'est pas sa réactivité, c'est le **crayon le
plus chaud**. Le facteur de point chaud radial

$$
F_{xy} = \frac{\text{puissance de l'assemblage le plus chargé}}{\text{puissance moyenne d'un assemblage}}
$$

fixe la marge avant crise d'ébullition (DNB). Un cœur au F_xy trop élevé
doit être bridé en puissance, quelle que soit sa réactivité.

Le résultat le plus parlant du code :

| chargement | k_eff | F_xy |
|---|---|---|
| uniforme (tout zone3) | 1,27317 | **2,123** |
| zoné « out-in » (1/2/3) | 1,23639 | **1,468** |

On place le combustible le **plus** réactif en **périphérie**, là où le flux
serait naturellement faible, et le moins réactif au centre. F_xy chute de
2,12 à 1,47.

Le prix : **3 020 pcm** de réactivité en moins, puisqu'on a remplacé du
combustible riche par du pauvre.

> **C'est le compromis central de la conception de cœur : on achète de la
> marge thermique avec de la réactivité.** Tout le reste en découle.

## 4. Contrôler la réactivité : trois outils, trois échelles de temps

Un cœur neuf a ~19 100 pcm d'excès. Il le faut : le combustible s'use, et il
faut de la réserve pour tenir le cycle. Mais il faut être critique
aujourd'hui.

### Le bore soluble — lent, uniforme

Dilué dans l'eau du primaire, donc **réparti uniformément** : il ne crée
aucun point chaud. C'est sa qualité maîtresse. Absorbant en 1/v, il n'agit
pratiquement que sur le groupe thermique.

Poids différentiel calculé, $$d\rho/d\text{ppm}$$ :

| ppm | k_eff | ρ (pcm) | différentiel |
|---|---|---|---|
| 0 | 1,23639 | 19 119 | — |
| 500 | 1,18062 | 15 299 | −7,64 pcm/ppm |
| 1 000 | 1,13061 | 11 552 | −7,49 |
| 1 500 | 1,08547 | 7 874 | −7,36 |
| 2 000 | 1,04452 | 4 263 | −7,22 |
| 2 500 | 1,00721 | 716 | −7,09 |

**Bore critique : 1 301 ppm.**

Le poids différentiel **décroît** en valeur absolue (−7,64 → −7,09). C'est
l'**auto-protection** : chaque atome de bore creuse le flux thermique, donc
les suivants voient moins de neutrons à absorber. Tout absorbant sature
ainsi — ce n'est pas propre au bore.

### La limite de sûreté du bore

Le bore est **dans l'eau**. Si la température monte, l'eau se dilate, il y a
moins d'eau — donc moins de bore — donc **gain** de réactivité. Cet effet
s'oppose au coefficient de température modérateur naturellement négatif.

Au-delà d'environ 1 800 ppm, la somme devient **positive** : une hausse de
température augmenterait la puissance, qui augmenterait la température.
Emballement. C'est la vraie limite haute sur la concentration en bore.

C'est exactement pourquoi les **poisons consommables** (gadolinium intégré
au combustible) existent : ils écrêtent l'excès de réactivité en début de
cycle sans passer par le bore, et disparaissent en s'usant.

### Les grappes de commande — rapides, locales

Elles tombent en secondes. Mais elles sont **locales**, donc elles
déforment la puissance : $$F_{xy}$$ passe de 1,485 à 1,592 à leur insertion.
On paie l'efficacité en point chaud.

Efficacité du groupe de 7 positions : **2 264 pcm**.

## 5. Ombrage et anti-ombrage : le piège des efficacités

Les efficacités **ne s'additionnent pas**. Somme des individuelles :
1 581 pcm. Groupe entier : 2 264 pcm. Écart **+43 %**.

J'avais d'abord écrit dans le code que c'était l'ombrage classique et que le
groupe valait donc *moins*. Le calcul disait l'inverse. L'expérience qui
tranche — deux grappes seules, puis ensemble :

| cas | seule A | seule B | somme | ensemble | écart |
|---|---|---|---|---|---|
| accolées (3,3)+(3,4) | 440 | 446 | 886 | 762 | **−14,0 %** |
| éloignées (1,1)+(5,4) | 159 | 171 | 329 | 351 | **+6,6 %** |
| éloignées (0,6)+(6,0) | 183 | 124 | 307 | 361 | **+17,6 %** |

Deux mécanismes opposés, et c'est la **distance** qui décide :

- **Accolées → ombrage (−14 %).** Chaque grappe creuse le flux thermique
  autour d'elle ; sa voisine baigne dans ce creux, voit moins de neutrons,
  en absorbe moins. Deux grappes collées se gênent.

- **Éloignées → anti-ombrage (+7 à +18 %).** Mécanisme différent. Une grappe
  **seule** fait basculer la nappe de flux à l'opposé d'elle-même : le cœur
  fuit latéralement la zone empoisonnée, ce qui limite sa propre efficacité.
  Avec des grappes réparties symétriquement, le flux n'a plus nulle part où
  basculer. Chacune conserve alors toute son efficacité.

Le groupe de 7 est réparti, donc l'anti-ombrage l'emporte.

**Conséquence pratique :** l'efficacité d'un groupe ne se déduit jamais de
mesures individuelles — ni par addition, ni par un facteur correctif
universel, puisque le **signe** de l'écart change avec la géométrie. C'est
pourquoi les essais physiques au démarrage mesurent l'efficacité **groupe
par groupe**, tel qu'ils seront réellement manœuvrés.

### Un piège de symétrie

La position (0,0) est au centre du cœur, flux maximal — et c'est la **moins
efficace** (82 pcm). Ce n'est pas un bug. Le plan de symétrie est au *bord*
de l'assemblage (0,0), pas en son milieu : en cœur complet, cette position
se reflète en un bloc serré de **2×2 grappes accolées**. Auto-ombrage
maximal.

Toute position en quart de cœur représente quatre assemblages physiques, et
**leur agencement compte autant que leur nombre**. À vérifier
systématiquement quand on exploite une symétrie.

---

## 6. Ce que le modèle ne contient PAS

Point le plus important du document. Le code calcule une neutronique
statique correcte, mais un cœur réel est gouverné par des effets absents
d'ici. Il ne faut pas lui faire dire ce qu'il ne sait pas.

### Les contre-réactions thermohydrauliques

Aucune. Les sections efficaces sont figées, alors qu'elles dépendent de la
température et de la densité de l'eau.

- **Effet Doppler** (coefficient combustible ; voir les résonances élargies de
  la [leçon 2]({{ '/posts/la-fission-de-l-uranium/' | relative_url }})). La température du combustible
  monte, les résonances de l'U238 s'élargissent, les captures stériles
  augmentent, la réactivité baisse. Toujours négatif, **instantané**.
  C'est la première barrière de sûreté d'un REP : c'est lui qui arrête une
  excursion de puissance avant toute action humaine ou mécanique.
- **Coefficient modérateur**. L'eau se dilate → moins de modération et moins
  de bore. Doit rester négatif, d'où la limite en bore de la section 4.

Sans ces effets, le modèle ne peut simuler ni un transitoire, ni la
stabilité du cœur.

### Le xénon 135

Absent, et c'est majeur pour l'exploitation. Produit de fission avec une
section de capture thermique gigantesque (~2,6 millions de barns), il
apparaît en partie par décroissance de l'iode 135 (période ~6,6 h).

Conséquences pratiques :
- **Équilibre xénon** : environ −2 800 pcm en marche stable.
- **Pic xénon après arrêt** : l'iode continue à décroître en xénon alors que
  le flux ne le détruit plus. Le xénon culmine ~10 h après l'arrêt, et le
  cœur peut devenir impossible à redémarrer pendant ~20 h. C'est un fait
  d'exploitation quotidien.
- **Oscillations xénon** : instabilité spatiale de période ~24 h dans les
  grands cœurs, qu'il faut piloter activement.

### L'épuisement (burnup)

Le combustible est neuf et le reste. Un cycle réel dure 12 à 18 mois,
pendant lesquels l'U235 se consomme, le plutonium se crée (et fournit
jusqu'à la moitié des fissions en fin de cycle), les produits de fission
s'accumulent. C'est ce qui explique l'écart section 4 : notre bore critique
de 1 301 ppm est celui d'un cœur **entièrement neuf**, alors qu'un REP réel
n'est rechargé qu'au tiers ou au quart.

### Les limites de la diffusion elle-même

Voir l'[annexe sur le solveur de diffusion]({{ '/posts/solveur-de-diffusion/' | relative_url }}), §1.2. La théorie de la diffusion est mauvaise près des
absorbants forts — c'est-à-dire précisément près des **grappes**. Les
efficacités calculées ici sont qualitativement justes et quantitativement
approximatives. Un calcul industriel utilise des facteurs de discontinuité
calibrés sur du transport.

---

## À retenir

- $$k_{\text{eff}} = k_\infty \times P_{\text{NL}}$$ : grâce au réflecteur, les fuites ne coûtent ici
  qu'environ 1 000 pcm.
- Ce qui limite un réacteur est le **crayon le plus chaud** : un chargement
  « out-in » fait passer $$F_{xy}$$ de 2,12 à 1,47, au prix de 3 020 pcm.
- Le **bore** contrôle lentement et uniformément (bore critique : 1 301 ppm
  pour ce cœur neuf) ; sa concentration est bornée par le signe du
  coefficient modérateur.
- Les **grappes** agissent vite mais localement ; leurs efficacités ne
  s'additionnent pas (ombrage si accolées, anti-ombrage si éloignées).
- Sans contre-réactions, xénon ni épuisement, ce modèle statique ne simule ni
  transitoire ni cycle.

---

## 7. Suite logique du développement

Par ordre de rapport pédagogique/effort :

1. **Contre-réactions Doppler et modérateur.** $$\Sigma_a(T)$$ et $$D(\rho)$$,
   bouclés par itération de Picard avec un modèle thermohydraulique
   canal-par-canal simple. C'est le pas qui transforme un calcul de
   neutronique en modèle de réacteur.
2. **Xénon et samarium.** Équations d'évolution I-135 / Xe-135 couplées au
   flux. Permet de reproduire le pic xénon — spectaculaire et très
   instructif.
3. **Épuisement.** Équations de Bateman (méthode CRAM pour l'exponentielle
   matricielle), pour suivre un cycle complet et voir la dilution du bore.
4. **Benchmarks publiés.** IAEA 2D, BIBLIS : géométries hétérogènes avec
   valeurs de référence, pour valider vraiment.
5. **Cinétique point.** Neutrons retardés, β_eff ≈ 650 pcm, équation de
   Nordheim. Indispensable pour comprendre pourquoi on raisonne en pcm et
   ce que signifie la criticité prompte.

---

## Bibliographie

{% include bibliographie.html %}

{% include cours-pied.html %}
