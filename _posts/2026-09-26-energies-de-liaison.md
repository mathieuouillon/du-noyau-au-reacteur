---
title: "Énergies de liaison : les données mesurées"
date: 2026-09-26 09:00:00 +0200
categories: ["I · Énergie du noyau et fission", "Énergie de liaison"]
tags: [AME2020, énergie de liaison, nombres magiques, goutte liquide]
description: "2825 nucléides de l'évaluation AME2020, leur provenance, et ce qu'ils révèlent : le maximum en Ni-62, la goutte liquide et les couches nucléaires."
image:
  path: /assets/img/nucleaire/energie_liaison_modeles.png
  alt: "Énergies de liaison : les données mesurées"
lecon: 1
partie: "I"
objectifs:
  - "Calculer une énergie de liaison à partir des masses atomiques mesurées (AME2020)."
  - "Lire la courbe B/A : son maximum en Ni-62, et pourquoi fission et fusion libèrent de l'énergie."
  - "Interpréter les termes de la goutte liquide et repérer ce qu'elle ignore."
  - "Mettre en évidence les nombres magiques directement dans les masses, avec S₂ₙ."
prerequis: []
code: [donnees_liaison.py, analyse_liaison.py, trace_liaison.py, energies_liaison.csv]
sources: >-
  Masses : évaluation AME2020 (fichier `massround.mas20.txt`, lu par le paquet `periodictable` 2.1.0) ; constantes : CODATA 2018 ; abondances : CIAAW, via le même paquet (uranium : Meija et al. 2016). Calculs : `donnees_liaison.py`, `analyse_liaison.py`, `trace_liaison.py`.
bibliographie:
  - cle: wang2021
    note: "La table de masses utilisée dans tout le cours."
  - cle: huang2021
    note: "Comment les mesures sont combinées en une évaluation."
  - cle: kondev2021
    note: "NUBASE2020 : demi-vies, spins, modes de désintégration."
  - cle: tiesinga2021
    note: "Les constantes (unité de masse atomique, masse du neutron)."
  - cle: meija2016
    note: "Compositions isotopiques recommandées."
  - cle: periodictable
    note: "Le paquet Python qui embarque AME2020."
  - cle: weizsacker1935
    note: "La formule semi-empirique de masse."
  - cle: bethe1936
    note: "La goutte liquide, présentée dans la « bible » de Bethe."
  - cle: mayer1949
    note: "Les nombres magiques expliqués par le spin-orbite."
  - cle: haxel1949
    note: "La même découverte, indépendante."
  - cle: nudat
    note: "Pour consulter un nucléide isolé."
  - cle: krane1988
    note: "Le manuel de référence : énergies de liaison, chapitre 3."
  - cle: valentin1982
    note: "Un manuel en français."
  - cle: basdevant2005
    note: "Un manuel moderne, de la structure du noyau à l'astrophysique."
---

{% include cours-entete.html %}

Tout le cours repose sur une seule grandeur : l'**énergie de liaison** `B`
d'un noyau, l'énergie qu'il faudrait fournir pour le séparer en nucléons
libres. Cette première leçon part des **masses mesurées**, avant tout modèle.

Jeu de données complet : **[`energies_liaison.csv`]({{ '/assets/code/energies_liaison.csv' | relative_url }})**, 2825 nuclides, Z de 1 à 100,
issus de l'évaluation **AME2020** — dont **2450 masses mesurées** et **375 masses
estimées** (voir « D'où viennent les données »).

| colonne | contenu |
|---|---|
| `Z`, `N`, `A` | protons, neutrons, masse |
| `symbole` | symbole chimique |
| `masse_atomique_u` | masse atomique (u) |
| `energie_liaison_MeV` | B totale |
| `B_sur_A_MeV` | **B/A** |
| `abondance_pct` | abondance naturelle (%) |
| `origine` | `mesuree` ou `estimee` (marqueur `#` d'AME2020) |

## D'où viennent les données

### La chaîne complète

```
mesures expérimentales de masses
  (spectromètres à pièges de Penning, anneaux de stockage,
   bilans Q de réactions et de désintégrations)
        │
        ▼
AME2020 — Atomic Mass Evaluation 2020
  ajustement global par moindres carrés de toutes les mesures publiées
        │
        ▼
fichier massround.mas20.txt  (AMDC, hébergé par l'AIEA)
        │  téléchargé le 2023-07-06
        ▼
paquet Python periodictable 2.1.0  (domaine public, P. Kienzle)
        │
        ▼
donnees_liaison.py  →  B, B/A, S_n, Q   →  energies_liaison.csv
```

### La source : AME2020

> M. Wang, W.J. Huang, F.G. Kondev, G. Audi, S. Naimi,
> *The AME 2020 atomic mass evaluation (II). Tables, graphs and references*,
> **Chinese Physics C 45** (2021) 030003.
> DOI : [10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf)

L'**Atomic Mass Evaluation** est la référence mondiale des masses atomiques,
rééditée régulièrement depuis les années 1950 (lignée initiée par
A.H. Wapstra, longtemps portée par G. Audi à Orsay). Elle est aujourd'hui
produite par l'Atomic Mass Data Center (AMDC), hébergé à l'Institute of Modern
Physics de Lanzhou ; les auteurs de l'édition 2020 sont rattachés notamment à
Lanzhou, Argonne, Orsay et RIKEN.

Ce n'est pas une compilation : c'est une **évaluation**. Plusieurs milliers de
mesures indépendantes — masses directes, mais aussi énergies de réactions
nucléaires et de désintégrations, qui relient les masses entre elles — sont
ajustées **ensemble** par moindres carrés. Chaque masse publiée bénéficie ainsi
de toutes les mesures qui la contraignent, directement ou par chaînage.

Le fichier exact est `massround.mas20.txt`, la version de la table publiée avec
l'article, aux incertitudes arrondies. Le paquet `periodictable` l'a téléchargé
le 6 juillet 2023 depuis le serveur de l'AMDC à l'AIEA
(`www-nds.iaea.org/amdc/ame2020/`) ; c'est ce que dit sa documentation, et je
l'ai vérifié en relisant le fichier embarqué.

### Mesurées ou estimées : la colonne `origine`

Point à ne pas manquer. Dans AME, **une masse suivie de `#` n'est pas purement
expérimentale** : pour les noyaux trop exotiques pour avoir été mesurés,
l'évaluation fournit une estimation tirée des tendances des noyaux voisins
(« systématique »). Ces valeurs sont utiles pour combler la carte, mais ce ne
sont pas des mesures.

| | nombre |
|---|---|
| nuclides retenus (Z ≤ 100) | 2825 |
| masses **mesurées** | **2450** |
| masses **estimées** (`#`) | **375** |

Le marqueur est conservé dans la colonne `origine` du CSV, et
`table_liaison(mesurees_seulement=True)` écarte les estimations.

Contrôles effectués :

- les **289 isotopes naturels** sont tous mesurés ;
- **tous les nuclides cités dans ce document** sont mesurés (U, Pu, Ni-62,
  Pd-118, Ba-141, Kr-92, les chaînes Ca, Sn et Pb des tableaux S₂ₙ) ;
- les statistiques sur la formule semi-empirique et le résidu de la figure
  sont calculés sur les **seules masses mesurées**.

> **Correction.** La première version de ce document parlait de « 2825
> nuclides mesurés » et comptait C-21 et N-24 parmi les pires échecs de la
> goutte liquide. Ce sont deux masses estimées. Les chiffres concernés ont été
> recalculés plus bas ; la conclusion ne change pas.

### Les abondances naturelles

Source : **CIAAW** (Commission on Isotopic Abundances and Atomic Weights de
l'IUPAC), *Isotopic compositions of the elements 2021*, via le même paquet.
Pour l'uranium, le CIAAW recommande la composition d'un minerai de référence
de Namibie (Meija *et al.*, *Pure Appl. Chem.* 88 (2016) 293).

> **Une anomalie du paquet.** Le fichier de `periodictable` contient bien les
> abondances de l'uranium (0,0054 / 0,7204 / 99,2742 %), mais son analyseur
> les perd : `periodictable.U[235].abundance` renvoie 0. Je les réinjecte dans
> [`donnees_liaison.py`]({{ '/assets/code/donnees_liaison.py' | relative_url }}) avec les valeurs exactes du fichier. C'est le **seul**
> élément concerné : les autres à abondance nulle (Tc, Pm, Po, At, Rn, Fr, Ra,
> Ac) n'ont effectivement aucun isotope stable. Une version précédente
> corrigeait aussi Th-232 à 100 % — à tort, le paquet donnait la bonne valeur
> (99,98 %) ; c'est rectifié.

### Les constantes

| constante | valeur | source |
|---|---|---|
| unité de masse atomique | 931,494 102 42 MeV | CODATA 2018 |
| masse du neutron | 1,008 664 915 95 u | CODATA 2018 |
| masse atomique de ¹H | 1,007 825 032 07 u | AME2020 |

### Le calcul

```
B(A,Z) = [ Z·m(¹H) + N·m(n) − M_atomique(A,Z) ] × 931,494 MeV/u
```

On utilise la masse **atomique** de ¹H et celle du nuclide : les Z électrons
se compensent. On néglige ainsi l'énergie de liaison électronique (~0,8 MeV
pour l'uranium, soit 4·10⁻⁴ de B).

### Validation

Recalcul confronté aux valeurs de B/A publiées :

| nuclide | calculé | référence | écart |
|---|---|---|---|
| He-4 | 7,0739 | 7,0739 | 0,02 keV |
| C-12 | 7,6801 | 7,6801 | 0,04 keV |
| Fe-56 | 8,7904 | 8,7903 | 0,06 keV |
| Ni-62 | 8,7946 | 8,7946 | −0,05 keV |
| Pb-208 | 7,8675 | 7,8675 | −0,05 keV |
| U-235 | 7,5909 | 7,5909 | 0,02 keV |

Accord à mieux que 0,1 keV — les écarts résiduels viennent de l'arrondi à
quatre décimales des valeurs de référence.

---

## Vue d'ensemble : mesures brutes et modèle

![Énergie de liaison par nucléon : mesures brutes et modèle de la goutte liquide]({{ '/assets/img/nucleaire/energie_liaison_modeles.png' | relative_url }})

*Généré par [`trace_liaison.py`]({{ '/assets/code/trace_liaison.py' | relative_url }}).*

### En haut : données brutes et modèle

Points gris : les 2450 masses mesurées d'AME2020. Croix orange : les 375
masses estimées (`#`), concentrées aux deux bords de la carte. Points bleus :
les 289 isotopes naturels, tous mesurés. La source est rappelée en pied de
figure. Courbe rouge : le modèle de la goutte liquide, tracé le long de la
vallée de stabilité qu'il prédit lui-même,

```
Z₀(A) = A / (2 + (a_C / 2a_A) · A^(2/3))
```

c'est-à-dire le Z qui maximise B à A fixé (on annule dB/dZ).

L'accord est remarquable au-delà de A ≈ 30 :

| nuclide | mesuré | modèle |
|---|---|---|
| Ni-62 | 8,7946 | 8,7820 |
| U-235 | 7,5909 | 7,5996 |
| **He-4** | **7,0739** | **5,4863** |

Le modèle place son maximum à A = 63, juste à côté du vrai maximum mesuré
(Ni-62).

Les deux flèches montrent pourquoi fission et fusion libèrent de l'énergie :

- **Fission** — couper U-236 en deux Pd-118 fait passer B/A de 7,587 à
  8,405 MeV, soit un gain de 0,8 MeV par nucléon, **193 MeV** au total sur
  236 nucléons.
- **Fusion** — réunir deux deutérons en He-4 fait passer B/A de 1,11 à
  7,07 MeV, soit **6 MeV par nucléon** : sept fois plus que la fission, par
  unité de masse.

### En bas à gauche : le modèle construit terme par terme

On part du terme de volume seul (15,75 MeV, constant, ligne pointillée), puis
on retranche les autres termes un à un. Chaque zone colorée est ce que coûte
un terme.

Le maximum de la courbe naît d'un **duel entre deux termes** :

- le terme de **surface** (zone verte) pénalise lourdement les petits noyaux,
  qui ont une grande part de leurs nucléons au bord, donc moins de voisins ;
- le terme de **Coulomb** (zone orange) pénalise les gros noyaux, car la
  répulsion entre protons croît en Z² alors que la force nucléaire ne lie
  qu'aux voisins immédiats.

Le terme d'asymétrie (zone rose) ne fait qu'affiner la forme. Le sommet se
trouve là où la surface cesse de dominer et où Coulomb n'a pas encore pris le
dessus. C'est aussi pour cela que l'uranium est presque instable : il est
situé loin dans la zone dominée par Coulomb, et il suffit de le déformer un
peu pour que la répulsion l'emporte.

### En bas à droite : ce que le modèle ignore

L'écart entre mesure et modèle, en keV par nucléon, en fonction du nombre de
neutrons N. Si le modèle était complet, on verrait du bruit autour de zéro.

Au lieu de ça, on voit des **bosses nettes exactement sur N = 50, 82 et 126** :
à ces nombres, les noyaux sont plus liés que prévu parce qu'ils ont une couche
de neutrons complète. La goutte liquide traite le noyau comme un fluide
continu, sans niveaux d'énergie individuels, et ne peut donc pas voir cet
effet. La section « Les nombres magiques se lisent dans les données brutes »
plus bas le confirme sans passer par aucun modèle.

À gauche, la zone dispersée correspond aux noyaux légers, où le modèle est le
moins fiable — c'est là que se trouve l'écart de 1,6 MeV sur He-4.

---

## Valeurs de référence

### Le maximum n'est pas le fer 56

| rang | nuclide | B/A (MeV) | abondance |
|---|---|---|---|
| 1 | **Ni-62** | **8,79455** | 3,63 % |
| 2 | Fe-58 | 8,79225 | 0,28 % |
| 3 | Fe-56 | 8,79036 | 91,75 % |
| 4 | Ni-60 | 8,78078 | 26,22 % |
| 5 | Cr-54 | 8,77797 | 2,37 % |

Fe-56 arrive **troisième**, à 4,2 keV/nucléon derrière Ni-62.

La confusion vient de deux critères différents. Fe-56 est le noyau de plus
faible **masse par nucléon** (m/A minimal), parce qu'il contient
proportionnellement plus de protons, plus légers que les neutrons. Ni-62 est
le plus **lié** (B/A maximal). C'est le premier critère qu'on cite
d'ordinaire en astrophysique, d'où le raccourci.

Fe-56 domine quand même dans la nature (91,75 % du fer) pour une raison
cinétique et non thermodynamique : la nucléosynthèse stellaire passe par
Ni-56, qui décroît en Fe-56.

### Noyaux utiles au réacteur

| nuclide | B/A (MeV) | B (MeV) |
|---|---|---|
| H-2 | 1,11228 | 2,22 |
| He-4 | 7,07392 | 28,30 |
| C-12 | 7,68014 | 92,16 |
| O-16 | 7,97620 | 127,62 |
| Fe-56 | 8,79036 | 492,26 |
| Ni-62 | 8,79455 | 545,26 |
| Kr-92 | 8,51272 | 783,17 |
| Ba-141 | 8,32614 | 1173,99 |
| Pb-208 | 7,86746 | 1636,43 |
| U-233 | 7,60434 | 1771,81 |
| U-235 | 7,59092 | 1783,87 |
| U-236 | 7,58654 | 1790,42 |
| U-238 | 7,57013 | 1801,69 |
| Pu-239 | 7,56030 | 1806,92 |

### Énergies de séparation du neutron — la fissilité

| capture | composé | S_n mesurée | barrière | fissile ? |
|---|---|---|---|---|
| U-235 + n | U-236 | **6,546 MeV** | 6,2 | **OUI** |
| U-238 + n | U-239 | **4,806 MeV** | 6,6 | non |
| Pu-239 + n | Pu-240 | **6,534 MeV** | 6,0 | **OUI** |
| U-233 + n | U-234 | **6,845 MeV** | 5,9 | **OUI** |

Les mesures confirment la prédiction : seuls les noyaux à nombre de neutrons
**impair** (U-235 : 143, Pu-239 : 145, U-233 : 141) libèrent assez d'énergie
d'appariement en capturant un neutron pour franchir leur barrière.

### Q de fission sur masses mesurées

| partition de U-236 | Q mesuré | Q goutte liquide |
|---|---|---|
| Pd-118 + Pd-118 (symétrique) | **193,2 MeV** | 185,0 |
| Xe-137 + Sr-97 + 2n | 177,3 | 165,2 |
| Ba-141 + Kr-92 + 3n | 166,7 | 156,4 |

---

## Fission symétrique ou asymétrique : ce que disent les masses

**L'asymétrie de la fission n'est pas une erreur du modèle.** J'avais d'abord
écrit dans la leçon sur la [fission de l'uranium]({{ '/posts/la-fission-de-l-uranium/' | relative_url }}) que la goutte liquide « prédit un maximum pour la
fission symétrique, alors que la fission réelle est asymétrique », en
attribuant l'écart aux couches nucléaires.

Les masses mesurées montrent que la partition symétrique libère **effectivement
plus d'énergie** (193 MeV contre 167). Sur ce point le modèle avait raison.

L'asymétrie n'est donc **pas** un effet de bilan énergétique. Elle vient de la
**dynamique** du noyau au point de scission : la surface d'énergie potentielle
au point selle favorise des fragments proches des couches fermées Z=50 et
N=82. Le noyau ne choisit pas le partage le plus exothermique, il suit le
chemin le plus facile. Les deux leçons sont corrigées.

## Qualité réelle de la formule semi-empirique

Masses mesurées seulement.

| échantillon | écart RMS |
|---|---|
| 2363 nuclides mesurés, A ≥ 20 | 0,0719 MeV/nucléon |
| 269 isotopes **naturels**, A ≥ 20 | **0,0270 MeV/nucléon** |

(En incluant les 374 masses estimées, le RMS global monterait à 0,0842 : les
estimations portent sur les noyaux les plus exotiques, là où le modèle est le
plus mauvais.)

Près de trois fois meilleur près de la vallée de stabilité — normal, c'est là que les
paramètres ont été ajustés.

Les huit pires cas sont tous des noyaux **légers et très riches en neutrons** :

| nuclide | Z | N | mesuré | SEMF | écart |
|---|---|---|---|---|---|
| C-22 | 6 | 16 | 5,421 | 4,194 | −1,227 |
| C-20 | 6 | 14 | 5,962 | 5,054 | −0,908 |
| N-23 | 7 | 16 | 6,237 | 5,329 | −0,907 |
| N-22 | 7 | 15 | 6,379 | 5,590 | −0,788 |
| O-26 | 8 | 18 | 6,497 | 5,729 | −0,768 |

Exactement le domaine des effets de couche et des peaux de neutrons, que le
modèle ignore par construction. **Sans conséquence pour le réacteur** :
l'uranium et ses produits de fission sont lourds et proches de la stabilité,
où le modèle vaut 0,03 MeV/nucléon.

> **Pour aller plus loin** : la [leçon 3]({{ '/posts/au-dela-de-la-goutte-liquide/' | relative_url }}) construit des modèles qui
> corrigent ces défauts — goutte liquide étendue, termes de couches,
> correction statistique — et les juge sur leur capacité à prédire.

## Les nombres magiques se lisent dans les données brutes

Comment prouver que le noyau a une structure en couches, **sans aucun modèle** ?
En suivant l'énergie de séparation de **deux** neutrons,
`S_2n = B(A,Z) − B(A−2,Z)`, le long d'une chaîne isotopique. On en prend deux
à la fois pour éliminer l'oscillation pair-impair.

| chaîne | N | S_2n (MeV) | chute |
|---|---|---|---|
| **Sn** (Z=50) | 81 | 12,82 | 0,10 |
| | 82 | 12,56 | 0,26 |
| | **83** | **9,75** | **2,81** ← |
| **Pb** (Z=82) | 125 | 14,82 | −0,01 |
| | 126 | 14,11 | 0,72 |
| | **127** | **11,31** | **2,80** ← |
| **Ca** (Z=20) | 19 | 30,29 | 1,46 |
| | 20 | 28,93 | 1,36 |
| | **21** | **24,00** | **4,93** ← |

Dans les trois chaînes, la cassure tombe exactement **un cran après** le
nombre magique. Tant que la couche se remplit, chaque paire de neutrons coûte
à peu près pareil ; dès qu'elle est pleine, la paire suivante tombe dans une
couche supérieure et l'énergie chute de 3 à 5 MeV d'un coup. C'est le saut
d'énergie d'ionisation entre un gaz rare et l'alcalin qui le suit.

Les nombres magiques 2, 8, 20, 28, 50, 82, 126 ne sont donc pas une hypothèse
théorique ajoutée après coup : **ils se lisent directement dans les masses
mesurées.**

> **Note de méthode.** J'ai d'abord cherché ces couches en moyennant le résidu
> (mesure − SEMF) par valeur de N. N=50 et N=82 ressortaient, mais N=20 et
> N=28 restaient noyés : à bas N l'échantillon contient surtout des noyaux
> très exotiques, où le modèle échoue pour d'autres raisons. Le S_2n le long
> d'une chaîne est supérieur sur les deux plans — il ne dépend d'aucun modèle,
> et il compare des noyaux voisins au lieu de mélanger tout un isotone. Le
> graphique des résidus (`courbe_liaison.png`, bas-gauche) montre d'ailleurs
> très bien les bosses à N=50, 82 et 126 : c'est le moyennage qui était en
> cause, pas l'indicateur.

---

## À retenir

- `B(A,Z) = [Z·m(¹H) + N·m(n) − M(A,Z)]·c²` : l'énergie de liaison se calcule
  directement à partir des masses atomiques mesurées (AME2020, 2450 masses
  mesurées sur 2825 nucléides retenus).
- `B/A` culmine à **Ni-62** (8,7946 MeV), pas à Fe-56 : couper un noyau lourd
  (fission) ou réunir des noyaux légers (fusion) libère de l'énergie.
- La **goutte liquide** (volume, surface, Coulomb, asymétrie, appariement)
  reproduit `B/A` à ~0,03 MeV par nucléon au voisinage de la vallée de
  stabilité ; le maximum
  naît du duel surface contre Coulomb.
- Ce qu'elle ignore se lit dans les données : des bosses en N = 50, 82, 126
  et des chutes de `S₂ₙ` juste après les **nombres magiques**. Le noyau a une
  structure en couches.

---

## Utilisation

```python
from donnees_liaison import table_liaison, B_sur_A, Q_fission, energie_separation_neutron

t = table_liaison()                        # 2825 nuclides (2450 mesurés)
tm = table_liaison(mesurees_seulement=True) # 2450, sans les estimations '#'
B_sur_A(t, "U", 235)                       # 7.59092 MeV
energie_separation_neutron(t, "U", 236)    # 6.546 MeV
Q_fission(t, ("U", 236), ("Ba", 141), ("Kr", 92))   # 166.7 MeV
```

Deux scripts régénèrent tout :

| commande | produit |
|---|---|
| `python analyse_liaison.py` | [`energies_liaison.csv`]({{ '/assets/code/energies_liaison.csv' | relative_url }}), `courbe_liaison.png` et l'analyse détaillée |
| `python trace_liaison.py` | `energie_liaison_modeles.png` (figure de la vue d'ensemble) |

Ils ont besoin de [`donnees_liaison.py`]({{ '/assets/code/donnees_liaison.py' | relative_url }}) et [`fission.py`]({{ '/assets/code/fission.py' | relative_url }}) dans le même dossier.

## Pour aller vers la source primaire

Les données utilisées **sont** celles d'AME2020, mais lues à travers une copie
embarquée dans un paquet Python plutôt que téléchargées à la source : le réseau
de cet environnement n'atteint pas les serveurs de l'AMDC ou du NNDC. Pour un
travail publiable, citer AME2020 et prendre directement :

- **AME2020** — la table pleine précision ou `massround.mas20.txt`, sur le
  site de l'AMDC (AIEA) ou du NNDC (Brookhaven). L'article
  *Chinese Physics C* 45 (2021) 030003 décrit le format et la signification
  du `#`.
- **AME2020 partie I** — Huang *et al.*, *Chinese Physics C* 45 (2021) 030002 :
  la méthode d'évaluation, utile pour comprendre comment les mesures sont
  combinées.
- **NUBASE2020** — Kondev *et al.*, *Chinese Physics C* 45 (2021) 030001 :
  demi-vies, spins, modes de désintégration.
- **KAERI Table of Nuclides** ou **NNDC NuDat** — consultation en ligne d'un
  nuclide isolé.

---

## Bibliographie

{% include bibliographie.html %}

{% include cours-pied.html %}
