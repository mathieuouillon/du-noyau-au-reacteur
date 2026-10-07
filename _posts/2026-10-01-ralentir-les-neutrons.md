---
title: "Ralentir les neutrons : modération et quatre facteurs"
date: 2026-10-01 09:00:00 +0200
categories: ["III · Du noyau au réacteur", "Neutronique"]
tags: [modération, quatre facteurs, uranium, sous-modération]
description: "Pourquoi il faut ralentir les neutrons, comment choisir un modérateur, et pourquoi un REP ne peut pas fonctionner à l'uranium naturel."
image:
  path: /assets/img/nucleaire/fission_bilan.png
  alt: "Ralentir les neutrons"
lecon: 7
partie: "III"
objectifs:
  - "Comparer les modérateurs : vitesse de ralentissement et rapport de modération."
  - "Écrire le facteur de multiplication infini avec la formule des quatre facteurs."
  - "Comprendre pourquoi le combustible est en crayons (autoprotection des résonances)."
  - "Expliquer la sous-modération des REP et pourquoi l'uranium naturel ne diverge pas dans l'eau légère."
prerequis: [2]
code: [fission.py, etude_fission.py]
sources: >-
  Sections efficaces à 2200 m/s, ν, ξ et intégrales de résonance : valeurs classiques saisies dans `fission.py`, sans référence précise dans le code ; pour des valeurs évaluées, voir ENDF/B-VIII.0 ou JEFF-3.3. Calculs : `fission.py`, `etude_fission.py`.
bibliographie:
  - cle: reuss2003
    note: "La meilleure entrée en matière en français ; couvre exactement cette leçon."
  - cle: bussac1985
    note: "La référence complète."
  - cle: duderstadt1976
    note: "Chapitres 2-3 et 10 : quatre facteurs et autoprotection."
  - cle: lamarsh1966
    note: "Un classique, très pédagogique sur le ralentissement."
  - cle: fermi1952
    note: "La pile de Chicago, racontée par Fermi."
  - cle: brown2018
    note: "Données nucléaires évaluées américaines."
  - cle: plompen2020
    note: "Données nucléaires évaluées européennes."
  - cle: nudat
    note: "Pour vérifier une section efficace ou une énergie."
---

{% include cours-entete.html %}

La [leçon 2]({{ '/posts/la-fission-de-l-uranium/' | relative_url }}) a montré pourquoi l'U-235 fissionne et combien
d'énergie sort. Reste à **entretenir** la réaction en chaîne : chaque fission
doit en provoquer au moins une autre. Cette leçon reprend la suite de
[`fission.py`]({{ '/assets/code/fission.py' | relative_url }}) et [`etude_fission.py`]({{ '/assets/code/etude_fission.py' | relative_url }}) ; tous les chiffres sont produits par le code.
[`diffusion.py`]({{ '/assets/code/diffusion.py' | relative_url }}) résout ensuite le transport des neutrons
([annexe]({{ '/posts/solveur-de-diffusion/' | relative_url }})) et [`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }}) conçoit le cœur ([leçon 8]({{ '/posts/neutronique-du-coeur/' | relative_url }})).

![Courbe de liaison, ralentissement, courbe de modération et quatre facteurs](/assets/img/nucleaire/fission_bilan.png)
_Figure produite par `etude_fission.py` : efficacité de ralentissement (en haut à droite), courbe de modération (en bas à gauche) et quatre facteurs (en bas à droite)._

---

## 1. Ralentir les neutrons : le choix qui décide de tout

Un neutron naît à ~2 MeV. La section de fission de l'U-235 vaut 583 barns à
0,025 eV contre ~1 barn à 2 MeV : **500 fois plus grande une fois
thermalisé**. Il faut donc ralentir — et vite, car entre les deux se trouvent
les résonances de capture de l'U-238 (6 à 200 eV).

| noyau | ξ | chocs 2 MeV → 0,025 eV |
|---|---|---|
| H | 0,9998 | **18** |
| D | 0,7219 | 25 |
| C | 0,1576 | 115 |
| U-238 | 0,0084 | 2 171 |

L'hydrogène est imbattable, et la raison est purement mécanique : neutron et
proton ont la même masse, donc un choc frontal transfère **toute** l'énergie.
C'est le billard. Sur l'U-238, le neutron rebondit comme une bille sur un mur.

Mais ralentir vite ne suffit pas — il faut ralentir **sans absorber** :

| modérateur | ξΣs [1/cm] | rapport ξΣs/Σa |
|---|---|---|
| eau légère | **0,987** | **62** |
| eau lourde | 0,162 | **4 720** |
| graphite | 0,060 | 214 |

**Le classement s'inverse.** L'eau légère ralentit 6 fois plus vite, mais son
rapport de modération est 76 fois plus mauvais, parce que l'hydrogène capture
les neutrons (H + n → deutérium).

Toute l'architecture des filières sort de ce tableau :

- **Eau légère** — ralentit vite → réseau compact, forte densité de
  puissance. Mais absorbe → il **faut** enrichir. C'est le REP.
- **Eau lourde** — n'absorbe presque pas → marche à l'uranium **naturel**,
  donc pas d'usine d'enrichissement. Mais ralentit mal → réseau très ouvert,
  cœur énorme. C'est CANDU.
- **Graphite** — intermédiaire, cœur encore plus grand. UNGG, RBMK, HTR.

Ce n'est pas un détail technique : c'est ce qui décide si un pays a besoin ou
non d'une usine d'enrichissement.

## 2. Les quatre facteurs

$$
k_\infty = \eta \cdot \varepsilon \cdot p \cdot f
$$

Le cycle du neutron, lu à l'envers : $$f$$ il est absorbé dans le combustible
plutôt qu'ailleurs · $$p$$ il a survécu aux résonances · $$\varepsilon$$ quelques fissions
rapides l'ont précédé · $$\eta$$ son absorption produit η nouveaux neutrons.

$$\eta$$ ne dépend **que** de l'enrichissement, jamais de la géométrie. C'est le
plafond absolu :

| enrichissement | η |
|---|---|
| 0,72 % (naturel) | 1,343 |
| 3,00 % (REP) | 1,838 |
| 100 % | 2,072 |

### Pourquoi le combustible est en crayons

Si l'on dissolvait l'uranium dans l'eau, à 3 % d'enrichissement :

| Vm/Vf | p homogène | p hétérogène | k_inf hom. | k_inf hét. |
|---|---|---|---|---|
| 1,0 | 0,546 | 0,600 | 1,004 | 1,102 |
| 2,0 | 0,651 | **0,775** | 1,162 | **1,383** |
| 3,0 | 0,702 | 0,843 | 1,219 | 1,465 |

Le réseau hétérogène gagne systématiquement. Deux auto-protections se
cumulent :

- **Énergétique** — à l'énergie exacte d'une résonance, les premiers atomes
  d'U-238 absorbent tout et creusent le flux ; les suivants ne voient plus
  rien. L'intégrale de résonance effective tombe de ~280 barns à dilution
  infinie à ~23 barns.
- **Spatiale** (hétérogène seulement) — les neutrons ralentissent dans
  l'**eau**, loin de l'uranium. Ils ne traversent la zone dangereuse 6-200 eV
  qu'à l'extérieur du crayon, où il n'y a pas d'U-238 pour les capturer.

La pile de Fermi, en 1942, était un empilement hétérogène de blocs d'uranium
dans du graphite. Homogénéisée, elle n'aurait jamais divergé.

## 3. La sous-modération : le choix de sûreté fondamental

Ajouter de l'eau fait monter $$p$$ (on ralentit mieux) et baisser $$f$$ (l'eau
absorbe). Il existe donc un optimum.

| Vm/Vf | p | f | k_inf |
|---|---|---|---|
| 1,0 | 0,600 | 0,971 | 1,102 |
| **2,0** | 0,775 | 0,944 | **1,383** ← REP réel |
| 3,0 | 0,843 | 0,918 | 1,465 |
| **4,4** | 0,890 | 0,884 | **1,489** ← optimum |
| 8,0 | 0,938 | 0,807 | 1,433 |

**Optimum à Vm/Vf = 4,4. Un REP fonctionne à 2,0 — nettement à gauche.**

Ce n'est pas une contrainte subie, c'est un **choix de sûreté**. Le REP est
volontairement **sous-modéré** : si la température monte, l'eau se dilate,
Vm/Vf diminue, et comme on est à gauche du maximum, $$k_\infty$$ **diminue**. La
puissance baisse d'elle-même. Le coefficient de température modérateur est
négatif, le réacteur est intrinsèquement stable.

À droite de l'optimum ce serait l'inverse : dilatation → plus de réactivité →
plus de puissance → plus de dilatation. Emballement. Un cœur sur-modéré serait
interdit.

C'est aussi la vraie raison de la limite en bore vue dans
[Neutronique du cœur : comprendre et concevoir]({{ '/posts/neutronique-du-coeur/' | relative_url }}) : le bore étant dans l'eau, trop de bore rend ce
coefficient positif.

## 4. Pourquoi un REP ne peut pas marcher à l'uranium naturel

| cas | Vm/Vf | p | f | k_inf | |
|---|---|---|---|---|---|
| U 3,0 % + H₂O (REP) | 2,0 | 0,775 | 0,944 | **1,383** | diverge |
| U naturel + H₂O | 2,0 | 0,770 | 0,846 | **0,901** | impossible |
| U naturel + H₂O très dilué | 6,0 | 0,917 | 0,647 | **0,820** | impossible |
| U naturel + D₂O, réseau serré | 2,0 | 0,230 | 0,999 | **0,318** | impossible |
| U naturel + D₂O (CANDU) | 16,0 | 0,832 | 0,996 | **1,147** | diverge |

**Résultat central : l'uranium naturel avec de l'eau ordinaire donne
k_inf < 1 quelle que soit la quantité d'eau.** Le réacteur est impossible, et
aucun raffinement de conception n'y changera rien. Il faut soit enrichir, soit
changer de modérateur.

Noter l'avant-dernière ligne : la **même** eau lourde dans un réseau serré
donne 0,318. On ne change pas de modérateur sans redessiner tout le réseau —
c'est l'erreur que j'ai faite au premier essai en imposant Vm/Vf = 2 à l'eau
lourde.

## 5. Bouclage : d'où viennent vraiment les constantes de [`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }})

[`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }}) utilisait $$\Sigma_{a2} = 0{,}082$$ et $$\nu\Sigma_{f2} = 0{,}1375$$ sans les justifier.
Tentative naïve, $$\Sigma = N\sigma$$ avec σ à 2200 m/s :

| | naïf | corrigé Maxwell | [`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }}) |
|---|---|---|---|
| Σa2 | 0,178 | 0,112 | 0,082 |
| νΣf2 | 0,327 | 0,206 | 0,1375 |

Facteur 2 d'écart au départ. Deux étapes manquaient :

**(a) Moyenne sur le spectre.** La valeur à 2200 m/s est une *convention*,
pas une moyenne. Les neutrons thermiques suivent une maxwellienne à la
température du modérateur ; pour une section en 1/v, la moyenne vaut
$$\sigma(2200)\cdot\frac{\sqrt{\pi}}{2}\cdot\sqrt{293{,}6/T}$$, soit **0,63** à 580 K.

**(b) Facteur de désavantage.** L'homogénéisation suppose le flux uniforme.
Il ne l'est pas : le crayon absorbe, donc le flux thermique y est déprimé.
Les sections homogénéisées correctes sont pondérées par le **flux**, pas par
le volume. Ce rapport vaut typiquement 1,1 à 1,3 dans un REP.

Après correction spectrale on est dans un facteur 1,3 des valeurs réelles,
avec un calcul de coin de table. C'est l'ordre de grandeur attendu — et c'est
exactement pourquoi la génération des constantes de groupe est un métier à
part entière (APOLLO, CASMO, WIMS) et non une multiplication $$N\sigma$$.

**La chaîne complète :**

```
données nucléaires évaluées (JEFF, ENDF/B)
  → traitement des résonances et autoprotection
  → calcul de réseau 2D fin (assemblage, multigroupe)
  → condensation à 2 groupes + homogénéisation pondérée par le flux
  → constantes utilisées par coeur.py
```

Chaque flèche est un domaine de recherche.

## À retenir

- Un neutron naît à ~2 MeV, mais la fission de l'U-235 est **500 fois** plus
  probable une fois le neutron thermalisé : il faut le ralentir sans le perdre
  dans les résonances de capture de l'U-238.
- L'hydrogène ralentit le mieux (18 chocs), mais absorbe : l'**eau légère**
  impose d'enrichir, l'**eau lourde** permet l'uranium naturel au prix d'un
  cœur bien plus grand.
- $$k_\infty = \eta\, \varepsilon\, p\, f$$ ; $$\eta$$ ne dépend que de l'enrichissement et fixe le
  plafond.
- Le réseau **hétérogène** (crayons dans l'eau) gagne sur le mélange homogène
  grâce à l'autoprotection énergétique et spatiale des résonances.
- Un REP est volontairement **sous-modéré** (Vm/Vf ≈ 2 contre un optimum à
  4,4) : son coefficient de température modérateur est négatif.
- Uranium naturel + eau légère donne $$k_\infty < 1$$ quelle que soit la quantité
  d'eau : il faut enrichir ou changer de modérateur.

---

## Bibliographie

{% include bibliographie.html %}

{% include cours-pied.html %}
