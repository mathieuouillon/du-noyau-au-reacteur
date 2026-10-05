---
title: "Au-delà de la goutte liquide : modèles de masse"
date: 2026-09-28 09:00:00 +0200
categories: ["II · Couches et formes des noyaux", "Modèles de masse"]
tags: [goutte liquide, couches, apprentissage automatique, extrapolation, peau de neutrons]
description: "Cinq modèles de masse, de la goutte liquide aux couches, jugés sur ce qu'ils prédisent : la physique extrapole, la statistique non."
image:
  path: /assets/img/nucleaire/modeles_masse.png
  alt: "Au-delà de la goutte liquide : modèles de masse"
lecon: 3
partie: "II"
objectifs:
  - "Construire une hiérarchie de modèles de masse, de la goutte liquide aux termes de couches."
  - "Juger un modèle sur ce qu'il prédit : ajustement, interpolation, extrapolation."
  - "Interpréter la symétrie de surface (peau de neutrons) et le terme de couches."
  - "Comprendre pourquoi une correction statistique n'extrapole pas, et ce qui manque encore : la déformation."
prerequis: [1]
code: [modeles_masse.py, etude_modeles.py]
---

{% include cours-entete.html %}

La [leçon 1]({{ '/posts/energies-de-liaison/' | relative_url }}) a montré ce que la goutte liquide ignore : les couches.
Cette leçon ouvre la partie II en la complétant. Cinq modèles de complexité croissante, tous
ajustés sur les **2367 masses mesurées** d'AME2020 (Z ≥ 8, N ≥ 8, sans les
valeurs estimées `#`). Tous les chiffres sont produits par [`etude_modeles.py`]({{ '/assets/code/etude_modeles.py' | relative_url }}).

![Comparaison des modèles de masse]({{ '/assets/img/nucleaire/modeles_masse.png' | relative_url }})

---

## La hiérarchie

| | modèle | ce qu'il ajoute | paramètres |
|---|---|---|---|
| **M0** | goutte liquide, manuel | — | 5 (fixés) |
| **M1** | goutte liquide réajustée | coefficients refaits sur AME2020 | 5 |
| **M2** | goutte liquide étendue | symétrie de surface, échange coulombien, Wigner | 8 |
| **M3** | M2 + couches | comptage des nucléons de valence | 11 |
| **M4** | M3 + régression à noyau | correction statistique du résidu | ~2400 |

M0 à M3 sont **linéaires** dans leurs coefficients : l'ajustement est un
moindres carrés exact, sans hasard ni réglage. M4 est un apprentissage
automatique.

### M2 : trois termes de physique en plus

- **Symétrie de surface**, `+(N−Z)²/A^(4/3)`. L'excès de neutrons coûte
  **moins cher en surface** qu'au cœur : il s'y loge en formant une peau de
  neutrons. Le coefficient de symétrie effectif d'un noyau fini est donc
  plus petit que celui de la matière nucléaire infinie.
- **Échange coulombien**, `+Z^(4/3)/A^(1/3)`. Correction quantique : le
  principe de Pauli tient les protons de même spin à distance, ce qui réduit
  un peu leur répulsion.
- **Terme de Wigner**, `−|N−Z|/A`. Surliaison des noyaux N ≈ Z, où protons et
  neutrons occupent les mêmes orbitales et interagissent au maximum.

### M3 : les couches, avec trois nombres

Pour chaque type de nucléon, situé dans une couche `[M_bas, M_haut[` de
dégénérescence `D`, avec `ν` nucléons de valence :

```
x = ν (D − ν) / D
```

`x` vaut 0 sur une couche fermée et culmine à mi-couche. Il compte en fait le
plus petit du nombre de particules ou de trous de valence (idée de Casten,
forme de Dieperink et Van Isacker). Avec `S = x_p + x_n` :

```
B_couches = a₁·S + a₂·S² + a₃·x_p·x_n
```

Les nombres magiques utilisés sont 2, 8, 20, 28, 50, 82, 126 et 184 — ce
dernier est une fermeture **prédite**, jamais observée, nécessaire pour les
neutrons au-delà de 126 (actinides).

### M4 : apprendre ce que la physique a raté

Le résidu de M3 n'est pas du bruit : c'est une structure lisse sur la carte
(Z, N) — déformations, couches secondaires (carte en bas à gauche de la
figure). Une **régression à noyau gaussien** l'apprend : chaque prédiction
est une moyenne pondérée des erreurs de M3 sur les noyaux voisins, avec un
poids `exp(−d²/2l²)`. La longueur de corrélation `l` et la régularisation
sont choisies par validation croisée **sur les seules données
d'entraînement** ; on obtient `l = 3` nucléons.

C'est la méthode de correction par fonctions de base radiales utilisée dans
la littérature pour améliorer FRDM, HFB et les modèles WS.

---

## Trois épreuves

Écart RMS sur l'énergie de liaison **totale** B (la métrique standard du
domaine), en MeV :

| modèle | paramètres | ajustement | interpolation | extrapolation | rapport |
|---|---|---|---|---|---|
| M0 | 5 | 3,76 | 3,76 | 5,47 | 1,46× |
| M1 | 5 | 3,08 | 3,09 | 4,33 | 1,40× |
| M2 | 8 | 2,47 | 2,48 | 2,82 | 1,14× |
| M3 | 11 | 1,21 | 1,22 | 1,42 | 1,16× |
| M4 | ~2400 | **0,23** | **0,28** | **0,49** | **1,71×** |

- **Ajustement** — erreur sur les noyaux qui ont servi à ajuster.
- **Interpolation** — validation croisée à 5 plis : on cache 20 % des noyaux
  tirés au hasard, on les prédit, et on recommence.
- **Extrapolation** — on retire les 2 isotopes les plus riches et les 2 plus
  pauvres en neutrons de chaque élément (372 noyaux), puis on les prédit.

**Trois leçons :**

1. **Pour M0-M3, ajustement ≈ interpolation.** Avec 5 à 11 paramètres pour
   2367 noyaux, le modèle ne peut pas mémoriser les données. Il n'y a pas de
   sur-apprentissage possible.
2. **Chaque pièce de physique paie.** Du manuel aux couches, l'erreur est
   divisée par 3. Le seul terme de couches (M2 → M3) la divise par 2 avec
   trois paramètres de plus.
3. **M4 écrase tout en interpolation mais se dégrade le plus en
   extrapolation** (1,71× contre 1,16× pour M3). Ce rapport est le chiffre à
   regarder quand on choisit un modèle pour **prédire**.

---

## Ce que disent les coefficients

### La goutte liquide réajustée (M1)

| terme | manuel | réajusté |
|---|---|---|
| volume | 15,750 | 15,570 |
| surface | 17,800 | 17,360 |
| Coulomb | 0,711 | 0,705 |
| asymétrie | 23,700 | 22,859 |
| appariement | 11,180 | 12,121 |

Les coefficients bougent peu : ceux du manuel étaient déjà ajustés, sur des
tables plus anciennes et avec `Z²` au lieu de `Z(Z−1)`. Le gain M0 → M1 vient
de là, pas d'une physique nouvelle.

### Le modèle complet (M3)

| terme | coefficient (MeV) | lecture |
|---|---|---|
| volume | 15,54 | liaison d'un nucléon au cœur du noyau |
| surface | 18,42 | coût des nucléons de surface |
| Coulomb | 0,723 | répulsion entre protons |
| asymétrie (J) | 28,88 | coût de l'excès de neutrons au cœur |
| appariement | 11,29 | bonus des noyaux pair-pair |
| symétrie de surface (Q) | 40,18 | ce que la surface rend sur l'asymétrie |
| échange coulombien | 1,16 | correction quantique à Coulomb |
| Wigner | 24,45 | surliaison des noyaux N ≈ Z |
| couches S | −1,340 | pénalité d'éloignement des couches fermées |
| couches S² | 0,0212 | courbure |
| couches p-n | 0,0606 | interaction p-n de valence |

**La symétrie.** Le volume coûte J = 28,9 MeV, mais la surface en rend
Q = 40,2 MeV. Le coefficient effectif d'un noyau fini est donc

```
a_sym(A) = J − Q · A^(−1/3)
```

soit **17,1 MeV pour le calcium 40** et **22,1 MeV pour le plomb 208**. Les
neutrons en excès se logent en surface, où ils coûtent moins cher : c'est la
**peau de neutrons**, mesurée sur le plomb 208 par l'expérience PREX et
directement liée à la physique des étoiles à neutrons.

> **Correction.** La première version du code présentait ce terme comme
> « l'asymétrie coûte plus en surface », avec le signe opposé. Le moindres
> carrés compensait en rendant le coefficient négatif : l'ajustement était
> juste, mon interprétation fausse. Le terme est maintenant écrit avec le
> signe physique et le coefficient sort positif.

**Les couches**, évaluées sur quelques noyaux :

| noyau | Z | N | x_p | x_n | terme de couches |
|---|---|---|---|---|---|
| Pb-208 | 82 | 126 | 0 | 0 | 0 |
| Sn-132 | 50 | 82 | 0 | 0 | 0 |
| Ca-48 | 20 | 28 | 0 | 0 | 0 |
| Sm-152 | 62 | 90 | 7,50 | 6,55 | −11,7 MeV |
| Er-166 | 68 | 98 | 7,88 | 10,18 | −12,4 MeV |
| U-238 | 92 | 146 | 7,73 | 13,10 | −12,6 MeV |

Les noyaux doublement magiques sont à zéro par construction ; ceux de milieu
de couche perdent plus de 10 MeV par rapport à eux. Le terme p-n est
**positif** : l'interaction entre protons et neutrons de valence rend de la
liaison, et c'est elle qui pousse les noyaux de milieu de couche à se
**déformer**.

---

## L'erreur en fonction de la distance au connu

L'épreuve qui compte vraiment. Pour chaque élément ayant au moins 14 isotopes
mesurés, on retire les **8 plus riches en neutrons**, on ajuste sur le reste,
puis on prédit les noyaux retirés. `d` = nombre de neutrons au-delà du dernier
isotope connu de l'élément. 696 noyaux prédits, sur 87 éléments.

| d | M2 | M3 | M4 |
|---|---|---|---|
| 1 | 2,95 | 1,27 | **0,32** |
| 2 | 3,08 | 1,31 | 0,31 |
| 4 | 3,38 | 1,42 | 0,54 |
| 6 | 3,76 | 1,55 | 0,78 |
| 8 | 4,35 | 1,77 | **1,10** |

C'est exactement la situation réelle : prédire la masse de noyaux de plus en
plus exotiques, vers la limite d'existence des noyaux (*drip line*), là où se
joue le processus r de nucléosynthèse des éléments lourds.

M4 est **quatre fois** meilleur que M3 à un pas du connu, mais son erreur
croît bien plus vite : ×3,5 de d = 1 à d = 8, contre ×1,4 pour M3. L'avantage
fond de 4× à 1,6× sur la plage testée.

C'est mécanique. La correction de M4 est une moyenne pondérée des erreurs des
voisins ; plus on s'éloigne, moins il y a de voisins, plus elle s'efface, et
plus il ne reste que la physique de M3. Elle ne s'efface pas tout à fait à
d = 8, car les éléments voisins (Z−1, Z+1) fournissent encore un peu
d'information en diagonale sur la carte.

> **La leçon générale.** Un modèle statistique ne sait rien au-delà de ses
> données ; seule la physique extrapole. C'est pourquoi les tables de masse
> utilisées en astrophysique reposent sur des modèles physiques, et que les
> corrections par apprentissage sont réservées au voisinage des noyaux
> mesurés.

---

## Où se situent ces modèles

Écart RMS sur les masses d'AME2020 de modèles publiés :

| modèle | type | RMS |
|---|---|---|
| Gogny D1M | HFB, interaction de portée finie | 0,81 MeV |
| Duflo-Zuker DZ10 | formule à 10 paramètres issue du modèle en couches | 0,56 MeV |
| HFB-27 | HFB avec interaction de Skyrme | 0,52 MeV |
| WS4 | macroscopique-microscopique Weizsäcker-Skyrme | ~0,3 MeV |
| réseaux de neurones | correction du résidu | ~0,2 MeV |
| *nos modèles* | | |
| M3 | goutte étendue + couches | 1,22 MeV |
| M4 | M3 + régression à noyau | 0,28 MeV |

M3 — 11 paramètres, une trentaine de lignes — fait deux fois l'erreur de
Duflo-Zuker avec une physique bien plus fruste. Pour descendre encore avec
de la physique, il manque deux choses :

- la **déformation** : la plupart des noyaux de milieu de couche sont en
  ballon de rugby, et l'énergie de déformation se chiffre en MeV. C'est la
  structure en bandes visible sur la carte du résidu de M3, entre les
  nombres magiques ;
- un vrai calcul des **niveaux individuels** (méthode de Strutinsky, ou
  Hartree-Fock-Bogoliubov) au lieu d'un simple comptage de valence.

> **Suite** : la [leçon 4]({{ '/posts/vers-frdm/' | relative_url }}) construit ce modèle macroscopique-microscopique
> (goutte déformable, niveaux de Nilsson, correction de Strutinsky, BCS).

M4 rivalise avec WS4 en interpolation, **mais la comparaison n'est pas
équitable** : WS4 est un modèle physique à une dizaine de paramètres, qui
extrapole bien mieux qu'une correction statistique. Sa valeur se juge sur
l'épreuve de la distance, pas sur le tableau ci-dessus.

---

## À retenir

- Chaque pièce de physique paie : de la goutte liquide de manuel (M0) au
  modèle avec couches (M3), l'erreur sur `B` passe de 3,8 à 1,2 MeV avec
  seulement 11 paramètres.
- La **symétrie de surface** traduit la **peau de neutrons** : l'excès de
  neutrons coûte moins cher en surface qu'au cœur.
- Le terme de couches est nul pour les noyaux doublement magiques et coûte
  plus de 10 MeV en milieu de couche ; l'interaction proton-neutron de valence
  pousse ces noyaux à se **déformer**.
- Un modèle se juge sur l'**extrapolation**, pas sur l'ajustement : la
  correction statistique (M4) excelle près des données mais se dégrade vite en
  s'en éloignant. Seule la physique extrapole.
- Pour aller plus loin, il faut calculer la **déformation** et les **niveaux
  individuels** : c'est l'objet de la leçon suivante.

---

## Utilisation

```python
from donnees_liaison import table_liaison
from modeles_masse import ModeleLineaire, ModeleCorrige, colonnes_M3

t = table_liaison(mesurees_seulement=True)
t = t[(t["Z"] >= 8) & (t["N"] >= 8)]
Z, N, B = t["Z"].astype(float), t["N"].astype(float), t["B_MeV"]

m3 = ModeleLineaire("M3", colonnes_M3).ajuster(Z, N, B)
m3.predire(92, 146)            # U-238, en MeV
dict(zip(m3.termes, m3.coefs)) # les 11 coefficients

m4 = ModeleCorrige().ajuster(Z, N, B)   # ~15 s
```

`python etude_modeles.py` refait toute l'étude et la figure (environ 2 min,
l'essentiel pour la validation croisée de M4).

## Références

- **AME2020** — M. Wang *et al.*, *Chinese Physics C* **45** (2021) 030003.
- **Termes de couches** — A.E.L. Dieperink, P. Van Isacker, *Eur. Phys. J. A*
  **42** (2009) 269 ; idée du comptage de valence : R.F. Casten.
- **Correction RBF** — N. Wang, M. Liu, *Phys. Rev. C* **84** (2011) 051303.
- **WS4** — N. Wang, M. Liu, X. Wu, J. Meng, *Phys. Lett. B* **734** (2014) 215.
- **Duflo-Zuker** — J. Duflo, A.P. Zuker, *Phys. Rev. C* **52** (1995) R23.
- **HFB-27** — S. Goriely, N. Chamel, J.M. Pearson, *Phys. Rev. C* **88**
  (2013) 061302.
- **FRDM2012** — P. Möller *et al.*, *At. Data Nucl. Data Tables* **109-110**
  (2016) 1. Modèle macroscopique-microscopique de référence, qui ajoute la
  déformation et la correction de Strutinsky.
- **Valeurs RMS sur AME2020** citées ci-dessus : comparaisons publiées pour
  WS4, HFB-27, D1M et DZ10 (A&A 2025 ; arXiv 2011.07904 ; arXiv 2609.19578),
  réseaux de neurones (« Machine learning the nuclear mass », 2021).

{% include cours-pied.html %}
