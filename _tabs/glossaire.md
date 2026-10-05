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
{% assign l6 = '/posts/ralentir-les-neutrons/' | relative_url %}
{% assign l7 = '/posts/neutronique-du-coeur/' | relative_url %}
{% assign la = '/posts/solveur-de-diffusion/' | relative_url %}

Les notions clés du cours, regroupées par thème. Le lien en fin de définition
renvoie à la leçon où la notion est introduite.

## Masses et énergie de liaison

Énergie de liaison `B`
: Énergie qu'il faudrait fournir pour séparer un noyau en nucléons libres.
  Elle se calcule à partir des masses mesurées :
  `B = [Z·m(¹H) + N·m(n) − M(A,Z)]·c²`. On la compare d'un noyau à l'autre
  par nucléon, `B/A`, maximal pour Ni-62. ([leçon 1]({{ l1 }}))

AME2020
: *Atomic Mass Evaluation* 2020, la référence mondiale des masses atomiques :
  un ajustement global de toutes les mesures publiées. Une valeur suivie de
  `#` est une estimation, pas une mesure. ([leçon 1]({{ l1 }}))

Goutte liquide (formule semi-empirique de masse)
: Modèle qui traite le noyau comme une goutte de liquide chargée, avec cinq
  termes : volume, surface, Coulomb, asymétrie et appariement. Excellent en
  moyenne, aveugle aux couches. ([leçon 1]({{ l1 }}), [leçon 3]({{ l3 }}))

Énergie de séparation `S_n`, `S₂ₙ`
: Énergie nécessaire pour arracher un neutron (ou deux) à un noyau. `S_n` du
  noyau composé décide de la fissilité ; la chute de `S₂ₙ` le long d'une
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
: Énergie qu'un noyau doit franchir pour se déformer jusqu'à la scission,
  quelques MeV pour les actinides (environ 6 MeV pour U-236). Pour les
  actinides, elle présente en réalité deux bosses, effet des couches.
  ([leçon 2]({{ l2 }}), [leçon 4]({{ l4 }}))

Paramètre de fissilité `x`
: `x = E_c / 2E_s`, rapport entre énergie coulombienne et énergie de surface
  d'une goutte sphérique. Pour `x ≥ 1`, la sphère est instable et la goutte
  fissionne spontanément ; U-236 est à 0,71. ([leçon 4]({{ l4 }}))

Q de fission
: Énergie libérée par une fission, différence entre les énergies de liaison
  des fragments et celle du noyau qui fissionne : environ 200 MeV pour
  l'uranium. ([leçon 2]({{ l2 }}))

Fission asymétrique
: L'uranium se casse le plus souvent en deux fragments inégaux (A ≈ 95 et
  A ≈ 139), alors que le partage symétrique libérerait plus d'énergie.
  L'asymétrie vient des couches Z = 50 et N = 82 des fragments, qui façonnent
  le chemin vers la scission. ([leçon 2]({{ l2 }}))

## Couches et formes des noyaux

Déformation `δ`, `β₂`
: Mesures de l'écart à la sphère. `δ > 0` : noyau **allongé** (ballon de
  rugby, *prolate*) ; `δ < 0` : **aplati** (galette, *oblate*). `β₂` est la
  convention expérimentale, déduite des transitions E2. ([leçon 4]({{ l4 }}))

Modèle de Nilsson
: Niveaux d'énergie d'un nucléon dans un oscillateur harmonique déformé, avec
  un terme de spin-orbite et un terme en `l²`. Son diagramme montre comment
  les couches se brisent quand le noyau se déforme. ([leçon 4]({{ l4 }}))

Spin-orbite
: Couplage entre le moment orbital `l` et le spin `s` d'un nucléon, de
  l'ordre du MeV dans le noyau. C'est lui qui produit les nombres magiques
  28, 50, 82, 126 (Goeppert Mayer, Haxel, Jensen, Suess, 1949).
  ([leçon 4]({{ l4 }}), [leçon 5]({{ l5 }}))

`Ω`
: Projection du moment angulaire d'un nucléon sur l'axe de symétrie d'un
  noyau déformé, le seul bon nombre quantique qui reste. Dans un noyau
  allongé, les orbitales de petit `Ω` sont les plus basses. ([leçon 5]({{ l5 }}))

Correction de couches (Strutinsky)
: Différence entre la somme des énergies des niveaux occupés et la même somme
  lissée : la seule partie des niveaux qui fluctue d'une couche à l'autre. On
  l'ajoute à la goutte liquide. ([leçon 4]({{ l4 }}))

BCS
: Théorie de l'appariement (Bardeen, Cooper, Schrieffer), transposée des
  supraconducteurs aux noyaux. Le **gap** `Δ` mesure l'appariement ; un
  nucléon célibataire **bloque** son niveau. ([leçon 4]({{ l4 }}))

Modèle macroscopique-microscopique, FRDM
: Goutte déformable plus corrections de couches et d'appariement calculées à
  partir des niveaux individuels. FRDM (*Finite-Range Droplet Model*, Los
  Alamos) en est la version de référence. ([leçon 4]({{ l4 }}))

Orbitale nucléaire
: Densité de probabilité `|ψ|²` d'un nucléon dans un état de Nilsson, tracée
  comme une orbitale électronique. La somme des orbitales occupées donne la
  forme du noyau. ([leçon 5]({{ l5 }}))

## Réacteur

Modérateur, rapport de modération
: Matériau qui ralentit les neutrons par chocs (eau, eau lourde, graphite).
  Il doit ralentir vite (`ξΣs` grand) **et** peu absorber (`ξΣs/Σa` grand).
  ([leçon 6]({{ l6 }}))

Formule des quatre facteurs
: `k_inf = η · ε · p · f` : facteur de reproduction, facteur de fission
  rapide, probabilité d'échapper aux résonances, facteur d'utilisation
  thermique. ([leçon 6]({{ l6 }}))

Sous-modération
: Réseau qui contient moins d'eau que l'optimum de réactivité. Un REP est
  volontairement sous-modéré pour que son coefficient de température
  modérateur soit négatif. ([leçon 6]({{ l6 }}))

`k_eff`, pcm
: Facteur de multiplication effectif, fuites comprises : le cœur est critique
  pour `k_eff = 1`. La réactivité `ρ = (k − 1)/k` se compte en pcm
  (pour cent mille). ([leçon 7]({{ l7 }}))

Facteur de point chaud `F_xy`
: Rapport entre la puissance de l'assemblage le plus chargé et la puissance
  moyenne. Il fixe la marge thermique d'un cœur. ([leçon 7]({{ l7 }}))

Ombrage, anti-ombrage
: Deux grappes accolées se gênent et absorbent moins que la somme de leurs
  efficacités (ombrage) ; éloignées, elles peuvent absorber plus (anti-ombrage).
  ([leçon 7]({{ l7 }}))

Flux scalaire, théorie de la diffusion
: Le flux `φ = n·v` est la longueur de trajectoire des neutrons par unité de
  volume et de temps ; les taux de réaction valent `Σ·φ`. La diffusion
  approche le transport par la loi de Fick `J = −D ∇φ`. ([annexe]({{ la }}))
