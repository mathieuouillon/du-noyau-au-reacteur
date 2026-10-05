---
title: "La fission de l'uranium : énergie, fissilité, fragments"
date: 2026-09-27 09:00:00 +0200
categories: ["I · Énergie du noyau et fission", "Fission"]
tags: [fission, uranium, appariement, Q de fission]
description: "Pourquoi un noyau lourd libère de l'énergie en se cassant, pourquoi l'U-235 fissionne avec des neutrons lents et pas l'U-238, ce que font les sections efficaces et les résonances, et où passent les 200 MeV."
image:
  path: /assets/img/nucleaire/fission_bilan.png
  alt: "La fission de l'uranium"
redirect_from:
  - /posts/de-l-uranium-au-reacteur/
lecon: 2
partie: "I"
objectifs:
  - "Expliquer, à partir de la courbe B/A, pourquoi un noyau lourd libère de l'énergie en se cassant."
  - 'Distinguer noyau fissile et fissionnable : noyau composé, appariement, \(S_n\) face à la barrière.'
  - "Comprendre pourquoi les neutrons lents sont si efficaces : taille quantique, résonances de Breit-Wigner, loi en 1/v."
  - "Comparer neutrons thermiques et rapides : sections efficaces et nombre de neutrons η par neutron absorbé."
  - "Relier la fission sous la barrière à l'effet tunnel, et le seuil de l'U-238 à la hauteur de la barrière."
  - "Calculer le Q d'une fission selon le partage, et l'origine coulombienne des 168 MeV des fragments."
prerequis: [1]
code: [fission.py, etude_fission.py, sections_efficaces.py]
sources: >-
  Énergies de liaison et Q : formule semi-empirique de `fission.py` (coefficients de manuel) et masses
  mesurées d'AME2020. Sections efficaces évaluées (fission, capture, ν, facteurs de Westcott) : JENDL-4.0
  à 300 K, lues dans les tables publiées par la JAEA. Résonance de l'U-238 à 6,67 eV : Mughabghab (2018).
  Barrières : RIPL-3, via Ryssens et al. (2023). Calculs : `sections_efficaces.py`. Les barrières typiques
  (6,2 et 6,6 MeV) et le bilan des ~200 MeV sont des valeurs classiques saisies dans `fission.py` et
  `etude_fission.py` sans référence précise.
bibliographie:
  - cle: fermi1935
    note: "La découverte de l'efficacité des neutrons lents."
  - cle: bohr1936
    note: "Le noyau composé."
  - cle: breit1936
    note: "La formule des résonances."
  - cle: shibata2011
    note: "JENDL-4.0 : les sections efficaces évaluées utilisées ici."
  - cle: jendl_tables
    note: "Les tables de valeurs thermiques et moyennes de la JAEA."
  - cle: mughabghab2018
    note: "Les paramètres de la résonance de l'U-238 à 6,67 eV."
  - cle: hill1953
    note: "La pénétrabilité d'une barrière parabolique."
  - cle: ryssens2023
    note: "Barrières de référence (RIPL-3)."
  - cle: hahn1939
    note: "La découverte : du baryum dans l'uranium irradié."
  - cle: meitner1939
    note: "L'interprétation : le noyau se casse en deux."
  - cle: bohr1939
    note: "La théorie de la fission par la goutte liquide."
  - cle: wang2021
    note: "Les masses mesurées, pour les énergies de séparation et les Q."
  - cle: madland2006
    note: "Le bilan détaillé de l'énergie libérée, mesuré et évalué."
  - cle: capote2009
    note: "Barrières de fission recommandées."
  - cle: lamarsh1966
    note: "Sections efficaces et résonances, très pédagogique."
  - cle: moller2001
    note: "Pourquoi la fission de l'uranium est asymétrique."
  - cle: scamps2018
    note: "Le rôle des fragments déformés en poire."
  - cle: vandenbosch1973
    note: "Le manuel classique sur la fission."
  - cle: wagemans1991
    note: "Un ouvrage collectif plus récent, très complet."
  - cle: krappe2012
    note: "La théorie moderne de la fission."
  - cle: reuss2003
    note: "La fission vue du neutronicien, en français."
  - cle: nudat
    note: "Pour vérifier énergies de liaison et sections efficaces."
---

{% include cours-entete.html %}

Les sections 1 à 3 de [`fission.py`]({{ '/assets/code/fission.py' | relative_url }}) et de
[`etude_fission.py`]({{ '/assets/code/etude_fission.py' | relative_url }}), complétées par
[`sections_efficaces.py`]({{ '/assets/code/sections_efficaces.py' | relative_url }}) pour les sections efficaces, le
passage de la barrière et le partage des fragments. Tous les chiffres sont produits par le code.
La suite de ces deux scripts, qui ralentit les neutrons et assemble un réseau de
réacteur, fait l'objet de la [leçon 7]({{ '/posts/ralentir-les-neutrons/' | relative_url }}).

![Courbe de liaison, ralentissement, courbe de modération et quatre facteurs]({{ '/assets/img/nucleaire/fission_bilan.png' | relative_url }})
_Figure produite par `etude_fission.py`. Le panneau en haut à gauche (courbe de liaison) concerne cette leçon ; les trois autres, la [leçon 7]({{ '/posts/ralentir-les-neutrons/' | relative_url }})._

---

## 1. Pourquoi un noyau lourd libère de l'énergie en se cassant

Une seule courbe explique la fission **et** la fusion : l'énergie de liaison
par nucléon $$B/A$$. Elle monte de l'hydrogène jusqu'à un maximum vers A ≈ 60
(**Ni-62, 8,7946 MeV** — pas Fe-56, voir [Énergies de liaison : les données mesurées]({{ '/posts/energies-de-liaison/' | relative_url }})), puis
redescend. Un noyau lourd qui se casse en deux morceaux plus proches du fer
devient plus lié ; la différence est libérée.

On la calcule par la formule semi-empirique de masse, qui traite le noyau
comme une goutte de liquide chargée.

| noyau | B/A calculé | B/A mesuré |
|---|---|---|
| He-4 | 5,486 | 7,074 |
| Fe-56 | 8,760 | 8,790 |
| U-235 | **7,600** | **7,591** |

Excellent sur les noyaux lourds, mauvais sur l'hélium 4 : le modèle ignore
les couches nucléaires, et He-4 est doublement magique donc anormalement
lié. Retenir cette limite, elle resservira.

**Le terme coupable** est le terme coulombien, en $$Z^2/A^{1/3}$$. Les protons
se repoussent **tous** mutuellement (en $$Z^2$$), alors que la force nucléaire
ne lie qu'aux voisins immédiats (en $$A$$). Passé le fer, la répulsion
l'emporte. L'uranium est déjà presque instable ; il suffit de le déformer.

## 2. Pourquoi l'U-235 et pas l'U-238

Les deux sont *fissionnables*. Un seul est **fissile**, c'est-à-dire cassable
par un neutron lent. La différence tient à l'appariement des nucléons, et le
modèle de la goutte liquide la reproduit :

| réaction | noyau composé | S_n calculé | barrière | fissile ? |
|---|---|---|---|---|
| U-235 + n → U-236* | U-236 | **6,70 MeV** | 6,2 MeV | **OUI** |
| U-238 + n → U-239* | U-239 | **4,87 MeV** | 6,6 MeV | non |

(mesures : 6,545 et 4,806 MeV — l'accord est excellent)

> Les barrières de ce tableau sont des valeurs typiques. La
> [leçon 6]({{ '/posts/barrieres-de-fission/' | relative_url }}) les calcule, et donne les valeurs recommandées par
> RIPL-3 pour U-236 : 5,67 MeV pour la bosse la plus haute, 5,00 MeV pour
> l'autre. La conclusion ne change pas : $$S_n$$ = 6,55 MeV passe au-dessus.
{: .prompt-info }

Le mécanisme en une phrase : **l'U-235 a un nombre impair de neutrons (143)**.
Le neutron incident vient compléter une paire, et l'énergie d'appariement
ainsi libérée suffit à franchir la barrière. L'U-238 a un nombre pair (146) ;
le neutron incident reste célibataire, ne reçoit pas ce bonus, et
l'excitation est insuffisante.

Un neutron thermique apporte 0,025 eV, soit rien du tout. **Tout vient de
S_n.** La fissilité ne dépend donc pas de la vitesse du neutron mais de la
**parité du noyau cible**. Même règle pour Pu-239 (N=145) et U-233 (N=141),
sans exception.

L'U-238 fissionne quand même si on lui apporte la différence en énergie
cinétique — seuil ~1 MeV. C'est le facteur $$\varepsilon \approx 1{,}03$$ de la [leçon 7]({{ '/posts/ralentir-les-neutrons/' | relative_url }}) : 3 % des fissions
d'un REP ont lieu sur l'U-238, par des neutrons encore rapides.

### Sept noyaux, une seule règle

La règle se vérifie sur tous les actinides courants. Le neutron absorbé forme
un **noyau composé** (Bohr, 1936) : il partage aussitôt son énergie entre tous
les nucléons, et le noyau « oublie » comment il a été formé. Son excitation
vaut

$$
E^* = S_n(\text{noyau composé}) + E_n\,\frac{A}{A+1}
$$

où le second terme est l'énergie cinétique du neutron dans le centre de masse.
Pour un neutron thermique, il est négligeable : tout vient de $$S_n$$. Avec les
masses mesurées (AME2020), les barrières recommandées (RIPL-3, quand le noyau
composé est pair-pair) et les sections efficaces évaluées (JENDL-4.0) :

| cible | N de la cible | noyau composé | $$S_n$$ (MeV) | barrière (MeV) | marge | fission à 0,0253 eV |
|---|---|---|---|---|---|---|
| Th-232 | 142 (pair) | Th-233 | 4,786 | — | — | 53,7 μb |
| U-238 | 146 (pair) | U-239 | 4,806 | — | — | 16,8 μb |
| Pu-240 | 146 (pair) | Pu-241 | 5,241 | — | — | 36,2 mb |
| Pu-241 | 147 (impair) | Pu-242 | 6,310 | 5,85 | +0,46 | 1 012 b |
| Pu-239 | 145 (impair) | Pu-240 | 6,534 | 6,05 | +0,48 | 747,4 b |
| U-235 | 143 (impair) | U-236 | 6,546 | 5,67 | +0,88 | 585,1 b |
| U-233 | 141 (impair) | U-234 | 6,845 | 5,50 | +1,35 | 531,3 b |

![Section efficace de fission thermique en fonction de l'énergie apportée par le neutron, pour sept actinides]({{ '/assets/img/nucleaire/fission_seuil.png' | relative_url }})
_Figure produite par [`sections_efficaces.py`]({{ '/assets/code/sections_efficaces.py' | relative_url }})._

Entre 5,2 et 6,3 MeV d'énergie apportée, la fission par neutron lent gagne
**sept ordres de grandeur**. Les cibles à N impair apportent 6,3 à 6,8 MeV,
au-dessus de la barrière : elles sont fissiles. Les cibles à N pair n'apportent
que 4,8 à 5,2 MeV : elles fissionnent à peine, et d'autant moins qu'elles sont
loin de la barrière.

**La signature de l'appariement dans les masses.** L'énergie de séparation du
neutron le long des isotopes de l'uranium suit une dent de scie :

| noyau | U-233 | U-234 | U-235 | U-236 | U-237 | U-238 | U-239 | U-240 |
|---|---|---|---|---|---|---|---|---|
| N | 141 | **142** | 143 | **144** | 145 | **146** | 147 | **148** |
| $$S_n$$ (MeV) | 5,762 | **6,845** | 5,298 | **6,546** | 5,126 | **6,154** | 4,806 | **5,929** |

Les noyaux à N pair lient leur dernier neutron 1,44 MeV de plus, en moyenne,
que leurs voisins à N impair : c'est l'énergie gagnée en complétant une paire,
environ deux fois le gap d'appariement $$\Delta \approx 12/\sqrt{A}$$ = 0,78 MeV.
C'est exactement ce que la goutte liquide encode dans son terme d'appariement,
et ce qui sépare un noyau fissile d'un noyau qui ne l'est pas.

---

## 3. Pourquoi des neutrons lents : la section efficace

L'énergie apportée décide **si** le noyau peut fissionner. Mais un réacteur a
aussi besoin que la capture du neutron soit **probable**. C'est la section
efficace qui le mesure, et c'est elle qui explique pourquoi les réacteurs
ralentissent leurs neutrons (Fermi, 1934-1935).

### Ce qu'est une section efficace

Un faisceau de neutrons de densité $$n$$ et de vitesse $$v$$ traverse une cible
de $$N$$ noyaux par unité de volume. Le nombre de réactions par unité de volume
et de temps est

$$
R = N\,\sigma\,n v = N\,\sigma\,\varphi
$$

La section efficace $$\sigma$$ a la dimension d'une **aire** : c'est la surface
apparente que le noyau présente au neutron pour une réaction donnée. Son unité,
le barn, vaut $$10^{-24}~\text{cm}^2 = 100~\text{fm}^2$$, à peu près l'aire
géométrique d'un noyau d'uranium.

### Un neutron lent est une onde étendue

Un neutron est aussi une onde, de longueur d'onde réduite

$$
ƛ = \frac{\hbar}{p} = \frac{\hbar c}{\sqrt{2 m_n c^2 E}}\cdot\frac{A+1}{A}
$$

(le facteur $$(A+1)/A$$ passe dans le centre de masse). Quand $$ƛ$$ dépasse la
taille du noyau, la notion de cible « géométrique » n'a plus de sens : l'aire
qui compte devient $$\pi ƛ^2$$.

| énergie du neutron | $$ƛ$$ | $$\pi ƛ^2$$ | rapport à $$\pi R^2$$ |
|---|---|---|---|
| 0,0253 eV (thermique) | 28 740 fm | 2,6 × 10⁷ b | 1,5 × 10⁷ |
| 1 eV | 4 571 fm | 6,6 × 10⁵ b | 3,8 × 10⁵ |
| 1 keV | 145 fm | 657 b | 381 |
| 1 MeV | 4,6 fm | 0,66 b | 0,38 |
| 2 MeV (fission) | 3,2 fm | 0,33 b | 0,19 |

(U-235 : $$R = 1{,}2\,A^{1/3}$$ = 7,41 fm, $$\pi R^2$$ = 1,72 b.)

![Taille quantique du neutron comparée à la taille du noyau d'U-235]({{ '/assets/img/nucleaire/sections_taille.png' | relative_url }})

Un neutron de fission (2 MeV) est plus petit que le noyau : sa section efficace
ne peut guère dépasser l'aire géométrique, et la fission de l'U-235 moyennée
sur le spectre de fission ne vaut que **1,218 b**. Un neutron
thermique « voit » une aire quinze millions de fois plus grande que le noyau.
La section efficace réelle (585,1 b) est loin de ce plafond, parce qu'il faut
encore que le neutron **entre** dans le noyau : c'est le rôle des résonances.

### Résonances : la formule de Breit et Wigner

Le noyau composé a des états excités bien définis. Quand l'énergie du neutron
tombe juste sur l'un d'eux, la capture devient énorme : c'est une
**résonance**. Pour une résonance isolée d'énergie $$E_0$$, Breit et Wigner
(1936) ont montré que

$$
\sigma_{n,x}(E) = \pi ƛ^2\, g\, \frac{\Gamma_n(E)\,\Gamma_x}{(E - E_0)^2 + \Gamma^2/4},
\qquad \Gamma = \Gamma_n + \Gamma_\gamma + \Gamma_f
$$

- $$\Gamma_n$$, $$\Gamma_\gamma$$, $$\Gamma_f$$ sont les **largeurs partielles** :
  les probabilités, par unité de temps, que le noyau composé réémette le
  neutron, émette un gamma, ou fissionne ;
- $$g = (2J+1)/\bigl(2(2I+1)\bigr)$$ est un facteur de spin, qui vaut 1 pour
  une cible de spin nul comme l'U-238 ;
- pour un neutron de moment orbital nul (onde s), $$\Gamma_n(E) \propto \sqrt{E}$$ :
  le neutron entre d'autant plus facilement qu'il a plus d'impulsion.

**Un exemple : la résonance de l'U-238 à 6,67 eV.** Avec $$E_0$$ = 6,674 eV,
$$\Gamma_n$$ = 1,48 meV et $$\Gamma_\gamma$$ = 23,0 meV (Mughabghab), la
capture atteint **22 349 b** au sommet, treize mille fois l'aire géométrique.
C'est cette résonance, et ses voisines, que les neutrons doivent traverser
pendant leur ralentissement sans être capturés
([leçon 7]({{ '/posts/ralentir-les-neutrons/' | relative_url }})). Prolongée jusqu'aux énergies thermiques, elle
donne une capture de 1,230 b à 0,0253 eV, soit **46 %** de la valeur évaluée
(2,683 b) : le reste vient des autres résonances.

![Résonance de l'U-238 et loi en 1/v de la fission de l'U-235]({{ '/assets/img/nucleaire/sections_resonance.png' | relative_url }})

### La loi en 1/v

Loin sous une résonance ($$E \ll E_0$$), le dénominateur $$(E - E_0)^2$$ ne
dépend presque plus de $$E$$. Il reste $$\pi ƛ^2 \propto 1/E$$ et
$$\Gamma_n \propto \sqrt{E}$$, d'où

$$
\sigma(E) \propto \frac{1}{E}\cdot\sqrt{E} = \frac{1}{\sqrt{E}} \propto \frac{1}{v}
$$

C'est la **loi en 1/v** : plus le neutron est lent, plus il passe de temps près
du noyau, et plus il a de chances d'être absorbé. Le modèle la vérifie
numériquement : $$\sigma\sqrt{E}$$ reste constant à 3 % près entre 0,001 et
0,1 eV. Le taux de réaction par neutron, $$\sigma v$$, est alors **indépendant
de la vitesse** : le nombre de réactions ne dépend que du nombre de neutrons
présents, pas de leur répartition en vitesse.

Les évaluations le mesurent par le facteur de Westcott $$g$$, rapport de la
moyenne maxwellienne à la valeur à 0,0253 eV, qui vaut 1 pour une loi en 1/v
exacte : 0,977 pour la fission de l'U-235 et 0,997 pour l'U-233, mais **1,058**
pour le Pu-239. Sa grande résonance de fission, à 0,3 eV, est déjà dans la
queue du spectre thermique : quand le modérateur chauffe, le spectre glisse
vers elle, et la fission comme la capture du Pu-239 augmentent.

### Lent contre rapide

| cible | fission à 0,0253 eV | fission rapide | rapport | $$\nu$$ | $$\eta$$ thermique | $$\eta$$ rapide |
|---|---|---|---|---|---|---|
| U-233 | 531,3 b | 1,908 b | 278 | 2,485 | **2,290** | ≥ 2,41 |
| U-235 | 585,1 b | 1,218 b | 480 | 2,436 | **2,084** | ≥ 2,27 |
| Pu-239 | 747,4 b | 1,802 b | 415 | 2,878 | **2,111** | ≥ 2,80 |
| Pu-241 | 1 012 b | 1,626 b | 622 | 2,947 | **2,169** | ≥ 2,82 |
| U-238 | 16,8 μb | 0,306 b | — | — | — | — |
| Th-232 | 53,7 μb | 0,080 b | — | — | — | — |
| Pu-240 | 36,2 mb | 1,328 b | — | — | — | — |

(JENDL-4.0. « Fission rapide » : moyenne sur le spectre des neutrons de fission
de l'U-235, d'énergie moyenne voisine de 2 MeV. $$\nu$$ : neutrons prompts et
retardés par fission. $$\eta = \nu\,\sigma_f/(\sigma_f + \sigma_\gamma)$$ :
neutrons émis par neutron absorbé. Le $$\eta$$ rapide utilise les sections
moyennées sur le spectre de fission et le $$\nu$$ thermique ; comme $$\nu$$ croît
avec l'énergie du neutron, c'est une borne basse.)

Trois conséquences :

- **Ralentir multiplie la fission de l'U-235 par 480.** C'est la raison d'être
  du modérateur. À composition égale, un cœur thermique a besoin de beaucoup
  moins de matière fissile qu'un cœur rapide.
- **En neutrons rapides, tout fissionne un peu**, y compris l'U-238 (0,306 b)
  et le Pu-240 (1,328 b), qui ne fissionnent presque pas en neutrons lents.
- **$$\eta$$ doit dépasser 2** pour qu'un combustible puisse à la fois entretenir
  la réaction et en régénérer un autre. En neutrons thermiques, l'U-233 est le
  meilleur ($$\eta$$ = 2,29) : c'est l'argument du cycle du thorium. Le Pu-239
  a un $$\eta$$ thermique médiocre (2,11), parce qu'il capture beaucoup. En
  neutrons rapides, sa capture s'effondre et son $$\eta$$ dépasse 2,8, loin
  devant l'U-235 : c'est le principe des surgénérateurs à neutrons rapides.

---

## 4. Franchir la barrière : effet tunnel et seuil

Une barrière n'est pas un mur. Un noyau composé un peu sous la barrière peut
fissionner par **effet tunnel**, avec une probabilité donnée par la formule de
Hill et Wheeler ([leçon 6]({{ '/posts/barrieres-de-fission/' | relative_url }})) :

$$
T = \frac{1}{1 + \exp\bigl[2\pi\,(B_f - E^*)/\hbar\omega\bigr]}
$$

Le rapport $$\sigma_f/\sigma_\gamma$$ mesure directement la compétition, dans le
noyau composé, entre fission et émission d'un gamma :

| cible | $$S_n$$ du composé (MeV) | $$\sigma_f/\sigma_\gamma$$ à 0,0253 eV |
|---|---|---|
| Th-232 | 4,786 | 7,3 × 10⁻⁶ |
| U-238 | 4,806 | 6,3 × 10⁻⁶ |
| Pu-240 | 5,241 | 1,3 × 10⁻⁴ |
| Pu-239 | 6,534 | 2,75 |
| Pu-241 | 6,310 | 2,79 |
| U-235 | 6,546 | 5,93 |
| U-233 | 6,845 | 11,7 |

Pour l'U-239, noyau composé de l'U-238, nous n'avons pas de barrière de
référence. Avec une barrière de 5,5 à 6,3 MeV (celles des actinides pair-pair
voisins) et une courbure $$\hbar\omega$$ de 0,5 à 1 MeV, la pénétrabilité à
$$E^* = S_n$$ = 4,806 MeV va de $$10^{-8}$$ à $$10^{-2}$$ : trop large pour en tirer
une valeur précise, mais le rapport mesuré, $$6 \times 10^{-6}$$, a bien l'ordre
de grandeur d'un passage par effet tunnel. Pour passer **par-dessus**, il faut
apporter $$B_f - S_n$$ en énergie cinétique, soit un neutron de **0,7 à 1,5 MeV**
selon la barrière : c'est le seuil observé, autour de 1 MeV.

Même au-dessus de la barrière, la fission n'est pas certaine : pour l'U-235,
un neutron thermique absorbé sur sept est simplement capturé (rapport 5,93),
et un sur quatre pour le Pu-239. Cette capture stérile est ce qui limite
$$\eta$$, et elle produit des noyaux (U-236, Pu-240) qui ne sont pas fissiles.

---

## 5. Combien d'énergie

$$Q = B(\text{fragment 1}) + B(\text{fragment 2}) - B(\text{U-236})$$. Les neutrons libres ont une
énergie de liaison nulle et ne comptent pas.

| partition de U-236 | Q (goutte liquide) |
|---|---|
| 118/46 + 118/46 (symétrique) | **185,0 MeV** |
| 133/52 + 100/40 + 3n | 164,6 |
| 141/56 + 92/36 + 3n | 156,4 |

**Quel partage ?** En balayant tous les partages de l'U-236 en deux
fragments, avec les masses mesurées et, pour chaque masse, la charge la plus
favorable, l'énergie disponible ne culmine pas au partage symétrique :

| fragment léger | fragment lourd | charges | Q (MeV) | Coulomb au contact (MeV) |
|---|---|---|---|---|
| A = 80 | A = 156 | 32 / 60 | 172,2 | 237,7 |
| A = 90 | A = 146 | 36 / 56 | 182,3 | 248,2 |
| A = 95 | A = 141 | 38 / 54 | 185,8 | 252,1 |
| A = 100 | A = 136 | 40 / 52 | 193,2 | 255,1 |
| **A = 104** | **A = 132** | **42 / 50** | **199,3** | 257,3 |
| A = 110 | A = 126 | 44 / 48 | 194,8 | 258,5 |
| A = 118 | A = 118 | 46 / 46 | 193,2 | 258,8 |

![Énergie disponible Q et répulsion coulombienne selon le partage de l'U-236]({{ '/assets/img/nucleaire/fission_partage.png' | relative_url }})
_Figure produite par [`sections_efficaces.py`]({{ '/assets/code/sections_efficaces.py' | relative_url }})._

Le maximum, 199,3 MeV, est atteint quand le fragment lourd est **¹³²Sn**,
doublement magique (Z = 50, N = 82) : les couches se lisent déjà dans le bilan
d'énergie. Mais le pic observé du fragment lourd est à A ≈ 139, pas à 132 :
le noyau ne choisit pas le partage le plus exothermique.

> **Correction.** Les versions précédentes de cette leçon et de la leçon 1
> affirmaient que « le partage symétrique libère le plus d'énergie », en
> comparant seulement le partage 118/118 (193 MeV) à Ba-141 + Kr-92 + 3n
> (167 MeV). Le balayage complet montre que le maximum est près du partage
> symétrique, mais au fragment ¹³²Sn. La conclusion ne change pas :
> l'asymétrie observée n'est pas un effet de bilan, elle vient du chemin vers
> la scission ([leçon 6]({{ '/posts/barrieres-de-fission/' | relative_url }})).
{: .prompt-warning }

**D'où viennent les 168 MeV.** L'énergie cinétique des fragments est leur
**répulsion coulombienne**, libérée après la rupture du col. Deux sphères au
contact se repousseraient de $$Z_1 Z_2 e^2 / d$$ avec
$$d = r_0\,(A_1^{1/3} + A_2^{1/3})$$ : pour Sr-95 + Xe-141, $$d$$ = 11,7 fm et
252 MeV, plus que le $$Q$$ disponible. Pour obtenir les 168 MeV observés, il faut
$$d$$ = 17,6 fm : à la scission, les centres des fragments sont ~6 fm plus
éloignés qu'au contact. Les fragments naissent donc **allongés**, reliés par un
col. Leur énergie de déformation, rendue après la rupture, les chauffe ; ils
l'évacuent en émettant les neutrons prompts (2,4 en moyenne pour l'U-235) et
des gammas.

Pour Ba-141 + Kr-92 + 3n : 156 MeV calculés contre 167 mesurés. 6 % d'erreur
pour un modèle à cinq paramètres.

**Les ~200 MeV, en détail :**

| forme | MeV | récupérable |
|---|---|---|
| énergie cinétique des fragments | 168 | oui |
| neutrons de fission | 5 | oui |
| gammas prompts | 7 | oui |
| bêta des produits de fission | 8 | oui |
| gammas des produits de fission | 7 | oui |
| **antineutrinos** | **12** | **NON** |
| captures (n,γ) hors fission | 10 | oui |

207 MeV libérés par fission, dont 12 partent en antineutrinos qui traversent
la Terre sans interagir. D'où les ~200 MeV récupérables conventionnels.

Ordre de grandeur : **1 g d'U-235 fissionné ≈ 1 MW·jour ≈ 2 à 3 tonnes de
charbon.** Le facteur million est simplement le rapport entre énergies
nucléaires (MeV) et chimiques (eV).

## À retenir

- La courbe $$B/A$$ monte jusqu'à Ni-62 puis redescend : couper un noyau lourd
  en deux fragments de masse moyenne libère environ 0,8 MeV par nucléon.
- Le responsable est le terme **coulombien** en $$Z^2/A^{1/3}$$ : la répulsion
  entre protons croît plus vite que la liaison nucléaire.
- Le neutron forme un **noyau composé** excité de $$E^* = S_n + E_n A/(A+1)$$.
  Cible à N impair (U-233, U-235, Pu-239, Pu-241) : $$S_n$$ = 6,3 à 6,8 MeV,
  au-dessus de la barrière, **fissile**. Cible à N pair (Th-232, U-238,
  Pu-240) : 4,8 à 5,2 MeV, en dessous. L'écart, 1,4 MeV, est l'appariement.
- Un neutron lent est une onde de 29 000 fm : il « voit » une aire quinze
  millions de fois celle du noyau. Sous les **résonances** (Breit-Wigner), sa
  section efficace suit la **loi en 1/v** : ralentir multiplie la fission de
  l'U-235 par 480.
- En neutrons rapides, la capture s'effondre : le $$\eta$$ du Pu-239 passe de
  2,11 à plus de 2,8, d'où les surgénérateurs. En thermique, l'U-233 est le
  meilleur ($$\eta$$ = 2,29).
- Sous la barrière, la fission passe par **effet tunnel** (U-238 : $$6 \times 10^{-6}$$
  fission par capture) ; par-dessus, il faut un neutron de ~1 MeV.
- Une fission libère ~207 MeV, dont ~200 MeV récupérables. Les 168 MeV des
  fragments sont leur répulsion coulombienne, libérée à ~18 fm : les fragments
  naissent allongés.
- L'énergie disponible culmine pour le fragment lourd **¹³²Sn**, doublement
  magique, mais la fission de l'uranium donne un fragment lourd vers A ≈ 139 :
  le noyau suit le chemin de la surface d'énergie potentielle, pas le bilan
  ([leçon 6]({{ '/posts/barrieres-de-fission/' | relative_url }})).

---

## Bibliographie

{% include bibliographie.html %}

{% include cours-pied.html %}
