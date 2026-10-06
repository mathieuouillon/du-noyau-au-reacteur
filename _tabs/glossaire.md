---
title: Glossaire
icon: fas fa-book
order: 1
---

{% assign l1 = '/posts/energies-de-liaison/' | relative_url %}
{% assign l2 = '/posts/la-fission-de-l-uranium/' | relative_url %}
{% assign l3 = '/posts/au-dela-de-la-goutte-liquide/' | relative_url %}
{% assign l4 = '/posts/vers-frdm/' | relative_url %}
{% assign l5 = '/posts/orbitales-nucleaires/' | relative_url %}
{% assign l6 = '/posts/barrieres-de-fission/' | relative_url %}
{% assign l7 = '/posts/ralentir-les-neutrons/' | relative_url %}
{% assign l8 = '/posts/neutronique-du-coeur/' | relative_url %}
{% assign la = '/posts/solveur-de-diffusion/' | relative_url %}

Les notions clés du cours, regroupées par thème. Le lien en fin de définition
renvoie à la leçon où la notion est introduite.

## Masses et énergie de liaison

Énergie de liaison $$B$$
: Énergie qu'il faudrait fournir pour séparer un noyau en nucléons libres.
  Elle se calcule à partir des masses mesurées :
  $$B = \bigl[Z\,m({}^{1}\mathrm{H}) + N\,m_n - M(A,Z)\bigr]\,c^2$$. On la compare d'un noyau à l'autre
  par nucléon, $$B/A$$, maximal pour Ni-62. ([leçon 1]({{ l1 }}))

AME2020
: *Atomic Mass Evaluation* 2020, la référence mondiale des masses atomiques :
  un ajustement global de toutes les mesures publiées. Une valeur suivie de
  `#` est une estimation, pas une mesure. ([leçon 1]({{ l1 }}))

Goutte liquide (formule semi-empirique de masse)
: Modèle qui traite le noyau comme une goutte de liquide chargée, avec cinq
  termes : volume, surface, Coulomb, asymétrie et appariement. Excellent en
  moyenne, aveugle aux couches. ([leçon 1]({{ l1 }}), [leçon 3]({{ l3 }}))

Énergie de séparation $$S_n$$, $$S_{2n}$$
: Énergie nécessaire pour arracher un neutron (ou deux) à un noyau. $$S_n$$ du
  noyau composé décide de la fissilité ; la chute de $$S_{2n}$$ le long d'une
  chaîne isotopique révèle les nombres magiques. ([leçon 1]({{ l1 }}))

Nombres magiques
: 2, 8, 20, 28, 50, 82, 126 : nombres de protons ou de neutrons qui ferment
  une couche et rendent le noyau particulièrement lié. Ils se lisent dans les
  masses mesurées ; le spin-orbite les explique. ([leçon 1]({{ l1 }}), [leçon 4]({{ l4 }}))

Appariement
: Surliaison des nucléons de même type groupés par paires de moments opposés.
  Elle crée l'effet pair-impair des masses et explique la fissilité de
  l'U-235. ([leçon 2]({{ l2 }}), [leçon 4]({{ l4 }}))

Peau de neutrons
: Excès de neutrons logé à la surface des noyaux riches en neutrons, où il
  coûte moins cher qu'au cœur. Elle se traduit par le terme de symétrie de
  surface des modèles de masse. ([leçon 3]({{ l3 }}))

Extrapolation
: Prédiction de noyaux éloignés de ceux qui ont servi à ajuster un modèle,
  vers la limite d'existence des noyaux (*drip line*). C'est l'épreuve qui
  sépare un modèle physique d'une correction statistique. ([leçon 3]({{ l3 }}))

## Fission

Fissile, fissionnable
: Un noyau **fissile** (U-235, Pu-239, U-233) fissionne après capture d'un
  neutron lent ; un noyau seulement **fissionnable** (U-238) a besoin d'un
  neutron rapide, au-dessus d'environ 1 MeV. ([leçon 2]({{ l2 }}))

Barrière de fission
: Hauteur du point selle au-dessus de l'état fondamental : l'énergie qu'un
  noyau doit acquérir pour se déformer jusqu'à la scission. Quelques MeV pour
  les actinides (5 à 6 MeV pour U-236). Pour eux, elle a deux bosses, effet
  des couches. ([leçon 2]({{ l2 }}), [leçon 6]({{ l6 }}))

Point selle
: Col du paysage d'énergie potentielle : minimum dans toutes les directions
  sauf celle du chemin de fission, où il est maximum. ([leçon 6]({{ l6 }}))

Double bosse, isomère de fission
: La correction de couches, qui oscille avec la déformation, découpe la
  barrière des actinides en deux bosses séparées par un second puits
  superdéformé. Un noyau piégé dans ce puits est un **isomère de fission** :
  il fissionne bien plus vite que l'état fondamental (Am-242, 14 ms,
  Polikanov 1962). ([leçon 6]({{ l6 }}))

Scission
: Rupture du col qui relie les deux fragments naissants : la fin du chemin de
  fission. ([leçon 6]({{ l6 }}))

Pénétrabilité (Hill-Wheeler)
: Probabilité de traverser une barrière par effet tunnel :
  $$T = 1/\bigl(1 + \exp[2\pi(B_f - E)/\hbar\omega]\bigr)$$ pour une barrière parabolique.
  ([leçon 6]({{ l6 }}))

Paramètre de fissilité $$x$$
: $$x = E_c / 2E_s$$, rapport entre énergie coulombienne et énergie de surface
  d'une goutte sphérique. Pour $$x \ge 1$$, la sphère est instable et la goutte
  fissionne spontanément. Pour U-236, x vaut 0,71 à 0,84 selon la goutte
  utilisée, et la barrière en dépend au cube. ([leçon 4]({{ l4 }}), [leçon 6]({{ l6 }}))

Q de fission
: Énergie libérée par une fission, différence entre les énergies de liaison
  des fragments et celle du noyau qui fissionne : environ 200 MeV pour
  l'uranium. ([leçon 2]({{ l2 }}))

Fission asymétrique
: L'uranium se casse le plus souvent en deux fragments inégaux (A ≈ 95 et
  A ≈ 139), alors que l'énergie disponible culmine plus près du partage
  symétrique, pour un fragment lourd ¹³²Sn.
  La goutte liquide prédit une fission symétrique : l'asymétrie vient des
  couches des fragments, qui façonnent le chemin vers la scission.
  ([leçon 2]({{ l2 }}), [leçon 6]({{ l6 }}))

Section efficace
: Aire apparente qu'un noyau présente à un neutron pour une réaction donnée :
  le taux de réaction vaut $$R = N\sigma\varphi$$. Unité : le barn,
  $$10^{-24}~\text{cm}^2$$. ([leçon 2]({{ l2 }}))

Noyau composé
: Noyau formé par l'absorption d'un neutron, excité de
  $$E^* = S_n + E_n A/(A+1)$$, qui « oublie » comment il a été formé avant de
  fissionner, d'émettre un gamma ou de réémettre un neutron (Bohr, 1936).
  ([leçon 2]({{ l2 }}))

Résonance, formule de Breit-Wigner
: Pic de section efficace quand l'énergie du neutron tombe sur un état du noyau
  composé ; sa forme est donnée par la formule de Breit et Wigner (1936), avec
  des largeurs partielles $$\Gamma_n$$, $$\Gamma_\gamma$$, $$\Gamma_f$$.
  ([leçon 2]({{ l2 }}))

Loi en 1/v
: Loin sous les résonances, la section efficace d'absorption décroît comme
  l'inverse de la vitesse du neutron : les neutrons lents sont absorbés bien
  plus souvent que les rapides. ([leçon 2]({{ l2 }}))

Évaluation (ENDF)
: Fichier qui rassemble, pour un noyau, toutes les sections efficaces
  recommandées, issues des mesures et des modèles (ENDF/B-VIII.0, JENDL,
  JEFF). Dans la région des résonances, il donne leurs paramètres, et la
  courbe se reconstruit (formalisme de Reich-Moore). ([leçon 2]({{ l2 }}))

Effet Doppler
: Élargissement des résonances par l'agitation thermique des noyaux : le pic
  s'abaisse, s'élargit, l'aire est conservée. Dans le combustible, il augmente
  la capture de l'U-238 quand la température monte : coefficient de
  température négatif. ([leçon 2]({{ l2 }}), [leçon 8]({{ l8 }}))

Facteur de reproduction $$\eta$$
: Nombre moyen de neutrons émis par neutron absorbé dans un noyau fissile,
  $$\eta = \nu\,\sigma_f/(\sigma_f + \sigma_\gamma)$$ : 2,08 pour l'U-235
  thermique, plus de 2,8 pour le Pu-239 rapide. ([leçon 2]({{ l2 }}), [leçon 7]({{ l7 }}))

## Couches et formes des noyaux

Déformation $$\delta$$, $$\beta_2$$
: Mesures de l'écart à la sphère. $$\delta > 0$$ : noyau **allongé** (ballon de
  rugby, *prolate*) ; $$\delta < 0$$ : **aplati** (galette, *oblate*). $$\beta_2$$ est la
  convention expérimentale, déduite des transitions E2. ([leçon 4]({{ l4 }}))

Modèle de Nilsson
: Niveaux d'énergie d'un nucléon dans un oscillateur harmonique déformé, avec
  un terme de spin-orbite et un terme en $$\boldsymbol{l}^2$$. Son diagramme montre comment
  les couches se brisent quand le noyau se déforme. ([leçon 4]({{ l4 }}))

Spin-orbite
: Couplage entre le moment orbital $$\boldsymbol{l}$$ et le spin $$\boldsymbol{s}$$ d'un nucléon, de
  l'ordre du MeV dans le noyau. C'est lui qui produit les nombres magiques
  28, 50, 82, 126 (Goeppert Mayer, Haxel, Jensen, Suess, 1949).
  ([leçon 4]({{ l4 }}), [leçon 5]({{ l5 }}))

$$\Omega$$
: Projection du moment angulaire d'un nucléon sur l'axe de symétrie d'un
  noyau déformé, le seul bon nombre quantique qui reste. Dans un noyau
  allongé, les orbitales de petit $$\Omega$$ sont les plus basses. ([leçon 5]({{ l5 }}))

Correction de couches (Strutinsky)
: Différence entre la somme des énergies des niveaux occupés et la même somme
  lissée : la seule partie des niveaux qui fluctue d'une couche à l'autre. On
  l'ajoute à la goutte liquide. ([leçon 4]({{ l4 }}))

BCS
: Théorie de l'appariement (Bardeen, Cooper, Schrieffer), transposée des
  supraconducteurs aux noyaux. Le **gap** $$\Delta$$ mesure l'appariement ; un
  nucléon célibataire **bloque** son niveau. ([leçon 4]({{ l4 }}))

Modèle macroscopique-microscopique, FRDM
: Goutte déformable plus corrections de couches et d'appariement calculées à
  partir des niveaux individuels. FRDM (*Finite-Range Droplet Model*, Los
  Alamos) en est la version de référence. ([leçon 4]({{ l4 }}))

Orbitale nucléaire
: Densité de probabilité $$\lvert\psi\rvert^2$$ d'un nucléon dans un état de Nilsson, tracée
  comme une orbitale électronique. La somme des orbitales occupées donne la
  forme du noyau. ([leçon 5]({{ l5 }}))

## Réacteur

Modérateur, rapport de modération
: Matériau qui ralentit les neutrons par chocs (eau, eau lourde, graphite).
  Il doit ralentir vite ($$\xi\Sigma_s$$ grand) **et** peu absorber ($$\xi\Sigma_s/\Sigma_a$$ grand).
  ([leçon 7]({{ l7 }}))

Formule des quatre facteurs
: $$k_\infty = \eta\, \varepsilon\, p\, f$$ : facteur de reproduction, facteur de fission
  rapide, probabilité d'échapper aux résonances, facteur d'utilisation
  thermique. ([leçon 7]({{ l7 }}))

Sous-modération
: Réseau qui contient moins d'eau que l'optimum de réactivité. Un REP est
  volontairement sous-modéré pour que son coefficient de température
  modérateur soit négatif. ([leçon 7]({{ l7 }}))

$$k_{\text{eff}}$$, pcm
: Facteur de multiplication effectif, fuites comprises : le cœur est critique
  pour $$k_{\text{eff}} = 1$$. La réactivité $$\rho = (k-1)/k$$ se compte en pcm
  (pour cent mille). ([leçon 8]({{ l8 }}))

Facteur de point chaud $$F_{xy}$$
: Rapport entre la puissance de l'assemblage le plus chargé et la puissance
  moyenne. Il fixe la marge thermique d'un cœur. ([leçon 8]({{ l8 }}))

Ombrage, anti-ombrage
: Deux grappes accolées se gênent et absorbent moins que la somme de leurs
  efficacités (ombrage) ; éloignées, elles peuvent absorber plus (anti-ombrage).
  ([leçon 8]({{ l8 }}))

Flux scalaire, théorie de la diffusion
: Le flux $$\varphi = n v$$ est la longueur de trajectoire des neutrons par unité de
  volume et de temps ; les taux de réaction valent $$\Sigma\varphi$$. La diffusion
  approche le transport par la loi de Fick $$J = -D\, \nabla\varphi$$. ([annexe]({{ la }}))
