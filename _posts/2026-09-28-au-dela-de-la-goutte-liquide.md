---
title: "Au-delà de la goutte liquide : modèles de masse"
date: 2026-09-28 09:00:00 +0200
categories: ["II · Couches et formes des noyaux", "Modèles de masse"]
tags: [goutte liquide, couches, apprentissage automatique, extrapolation, peau de neutrons]
description: "Cinq modèles de masse, de la goutte liquide aux couches : la physique de chaque terme, les mathématiques de l'ajustement, et un jugement sur ce qu'ils prédisent."
image:
  path: /assets/img/nucleaire/modeles_masse.png
  alt: "Au-delà de la goutte liquide : modèles de masse"
lecon: 3
partie: "II"
objectifs:
  - "Retrouver l'origine physique de chaque terme de la goutte liquide : saturation, tension de surface, sphère chargée, gaz de Fermi, appariement."
  - "Construire une hiérarchie de modèles de masse, jusqu'à la symétrie de surface et aux termes de couches."
  - "Écrire un ajustement par moindres carrés et lire, dans les corrélations, ce que les masses contraignent mal."
  - "Juger un modèle sur ce qu'il prédit : ajustement, interpolation, extrapolation."
  - "Comprendre la régression à noyau, et pourquoi une correction statistique n'extrapole pas."
prerequis: [1]
code: [modeles_masse.py, etude_modeles.py]
sources: >-
  Masses : les 2367 masses mesurées d'AME2020 avec Z, N ≥ 8 (sans les valeurs estimées `#`). Calculs : `modeles_masse.py`, `etude_modeles.py`. Écarts RMS des modèles publiés : comparaisons sur AME2020 publiées pour WS4, HFB-27, D1M et DZ10 (A&A 2025 ; arXiv 2011.07904 ; arXiv 2609.19578) et pour les réseaux de neurones (« Machine learning the nuclear mass », 2021), telles que relevées par l'auteur.
bibliographie:
  - cle: weizsacker1935
    note: "La formule à cinq termes."
  - cle: bethe1936
    note: "La goutte liquide et le terme d'asymétrie."
  - cle: wigner1937
    note: "La symétrie de spin-isospin, origine du terme de Wigner."
  - cle: myers1966
    note: "La goutte liquide de Myers et Swiatecki et ses corrections de couches."
  - cle: myers1969
    note: "Le modèle de la gouttelette : symétrie de surface et peau de neutrons."
  - cle: dieperink2009
    note: "Le terme de couches par comptage de valence utilisé ici (M3)."
  - cle: duflo1995
    note: "Une formule de masse à dix paramètres issue du modèle en couches."
  - cle: wang2014
    note: "WS4, modèle macroscopique-microscopique de référence."
  - cle: wang2011
    note: "La correction par fonctions de base radiales, modèle de M4."
  - cle: goriely2013
    note: "HFB-27 : une table de masses entièrement microscopique."
  - cle: goriely2009
    note: "Le premier modèle de masse Gogny-HFB (D1M)."
  - cle: moller2016
    note: "FRDM(2012), le modèle de la leçon suivante."
  - cle: adhikari2021
    note: "La mesure de la peau de neutrons du plomb 208."
  - cle: wang2021
    note: "Les masses utilisées."
  - cle: hastie2009
    note: "Moindres carrés, validation croisée, régression ridge et à noyau."
  - cle: rasmussen2006
    note: "La régression à noyau gaussien, vue comme un processus gaussien."
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
automatique. Les sections 1 à 3 détaillent la physique de chaque terme, les
sections 4 et 5 les mathématiques de l'ajustement.

En notant $$I = (N-Z)/A$$, et en comptant $$B > 0$$ pour un noyau lié, M2 s'écrit :

$$
\begin{aligned}
B(Z,N) = {} & a_V A - a_S A^{2/3} - a_C \frac{Z(Z-1)}{A^{1/3}} - J \frac{(N-Z)^2}{A} + \delta_{\text{pair}} \\
& + Q \frac{(N-Z)^2}{A^{4/3}} + c_{\text{ex}} \frac{Z^{4/3}}{A^{1/3}} - W \frac{\lvert N-Z \rvert}{A}
\end{aligned}
$$

---

## 1. La goutte liquide, terme par terme

**Volume : la saturation.** La force nucléaire est attractive mais de
**courte portée** (~1 fm) et répulsive à très courte distance. Chaque nucléon
n'interagit donc qu'avec un nombre fixe de voisins, et la densité au centre des
noyaux est la même pour tous : $$\rho_0 \approx 0{,}16~\text{nucléon/fm}^3$$. L'énergie de liaison
est alors proportionnelle au nombre de nucléons, $$a_V A$$, et le rayon croît
comme $$R = r_0 A^{1/3}$$. Si la force agissait entre toutes les paires, comme la
gravitation, l'énergie croîtrait comme $$A^2$$ et les noyaux s'effondreraient.

**Surface : la tension superficielle.** Un nucléon de surface a moins de
voisins : on retire une énergie proportionnelle à l'aire, $$4\pi R^2 \sigma = a_S A^{2/3}$$.
Avec $$a_S$$ = 17,36 MeV (M1) et le $$r_0$$ déduit ci-dessous, la **tension de
surface** nucléaire vaut $$\sigma = a_S / 4\pi r_0^2$$ = **0,92 MeV/fm²**. C'est la
grandeur qui fixe le coût d'une déformation et la barrière de fission
([leçon 6]({{ '/posts/barrieres-de-fission/' | relative_url }})).

**Coulomb : la sphère chargée.** Construisons une sphère de charge $$Q$$
uniforme couche par couche. Quand le rayon vaut $$r$$, la charge déjà en place
est $$q(r) = Q\,(r/R)^3$$, et apporter la couche suivante $$dq = 3Q\, r^2\, dr / R^3$$
coûte $$q(r)\, dq / r$$ :

$$
E_C = \int_0^R \frac{q(r)\, dq}{r} = \frac{3Q^2}{R^6} \int_0^R r^4\, dr = \frac{3}{5}\,\frac{Q^2}{R}
$$

Avec $$Q = Ze$$, on remplace $$Z^2$$ par $$Z(Z-1)$$ : un proton ne se repousse
pas lui-même. Donc $$a_C = \tfrac{3}{5}\, e^2/r_0$$, avec $$e^2 = 1{,}44~\text{MeV·fm}$$. Le
coefficient ajusté de M1, $$a_C$$ = 0,705 MeV, donne $$r_0$$ = **1,225 fm**. On
peut le comparer au rayon de charge mesuré du plomb 208, 5,50 fm en moyenne
quadratique (Angeli et Marinova) : la sphère uniforme équivalente a pour rayon
$$\sqrt{5/3} \times 5{,}50 = 7{,}10~\text{fm}$$, soit $$r_0 = 7{,}10 / 208^{1/3}$$ = 1,20 fm. Un
ajustement sur les seules masses retrouve donc la taille des noyaux à 2 % près.

**Asymétrie : le principe de Pauli.** Les protons et les neutrons remplissent
chacun leurs niveaux jusqu'à un niveau de Fermi. Dans le modèle du gaz de
Fermi, $$n$$ nucléons d'une espèce dans un volume $$V$$ ont une énergie cinétique
$$\tfrac{3}{5}\, n\, \varepsilon_F(n)$$, avec $$\varepsilon_F \propto (n/V)^{2/3}$$. Écrivons $$N = \tfrac{A}{2}(1+I)$$ et
$$Z = \tfrac{A}{2}(1-I)$$ :

$$
\begin{aligned}
E_{\text{cin}} &= \frac{3}{5}\,\frac{A}{2}\,\varepsilon_F \Bigl[ (1+I)^{5/3} + (1-I)^{5/3} \Bigr] \\
&= \frac{3}{5}\,\varepsilon_F A + \frac{\varepsilon_F}{3}\, A I^2 + O(I^4)
\end{aligned}
$$

Le second terme est le terme d'asymétrie, $$(\varepsilon_F/3)\,(N-Z)^2/A$$. À la densité
de saturation, $$k_F = (3\pi^2 \rho_0/2)^{1/3}$$ = 1,333 fm⁻¹ et
$$\varepsilon_F = \hbar^2 k_F^2 / 2m$$ = 36,8 MeV : la seule énergie cinétique donne **12,3 MeV**.
L'ajustement en demande 22,9 (M1), et même $$J$$ = 28,9 MeV une fois la surface
séparée (M3). Le reste vient de l'**interaction** : la force entre un proton
et un neutron est plus attractive qu'entre deux nucléons identiques (le
deutéron est lié, le dineutron ne l'est pas), et un noyau $$N = Z$$ en profite
au maximum.

**Appariement.** Les nucléons identiques gagnent de l'énergie à se grouper
par paires de moments opposés (leçon 4, BCS). D'où un terme
$$\delta_{\text{pair}} = +a_P / A^{1/2}$$ pour les noyaux pair-pair, 0 pour $$A$$ impair, et
$$-a_P / A^{1/2}$$ pour les noyaux impair-impair. La loi en $$A^{-1/2}$$ est
empirique : c'est celle du gap moyen $$\Delta \approx 12/\sqrt{A}$$ MeV de Bohr et Mottelson,
et l'ajustement redonne bien $$a_P$$ ≈ 12 MeV.

---

## 2. La goutte étendue (M2)

**Le développement leptoderme.** Un noyau lourd est une goutte « à peau
fine » : l'épaisseur de surface (~2 fm) est petite devant le rayon. Toute
grandeur extensive se développe alors en puissances de $$A^{-1/3}$$, un terme
de volume ($$\propto A$$), puis un terme de surface ($$\propto A^{2/3}$$), puis de courbure
($$\propto A^{1/3}$$). La goutte à cinq termes applique ce développement à l'énergie
de liaison, mais **pas** à l'énergie d'asymétrie, qui n'a qu'une partie de
volume.

**Symétrie de surface.** Appliquer le même développement à l'asymétrie
donne :

$$
\begin{aligned}
E_{\text{sym}} &= \bigl(J A - Q A^{2/3}\bigr)\, I^2 = J\,\frac{(N-Z)^2}{A} - Q\,\frac{(N-Z)^2}{A^{4/3}} \\
&= a_{\text{sym}}(A)\,\frac{(N-Z)^2}{A}, \qquad a_{\text{sym}}(A) = J - Q\,A^{-1/3}
\end{aligned}
$$

$$Q > 0$$ : l'asymétrie coûte **moins** cher en surface. L'excès de neutrons
s'y loge et forme une **peau de neutrons**. Le coefficient effectif d'un
noyau fini est donc plus petit que celui de la matière nucléaire infinie, $$J$$.

**Échange coulombien.** Le principe de Pauli tient les protons de même spin
à distance : ils se repoussent un peu moins que des charges classiques. Pour
une sphère uniforme, l'approximation de Slater donne :

$$
\begin{gathered}
E_{\text{ex}} = -\frac{3}{4}\left(\frac{3}{2\pi}\right)^{2/3} \frac{e^2 Z^{4/3}}{R} \\[4pt]
\Longrightarrow\quad c_{\text{ex}} = \frac{3}{4}\left(\frac{3}{2\pi}\right)^{2/3} \frac{e^2}{r_0} = 0{,}54~\text{MeV}
\end{gathered}
$$

L'ajustement donne 0,96 MeV (M2) et 1,16 MeV (M3), environ le double. Le
terme en $$Z^{4/3}/A^{1/3}$$ absorbe donc d'autres effets que la seule forme
fonctionnelle ne permet pas de distinguer, comme la diffusivité de la surface
de charge. C'est un premier exemple de ce que montre la section 4 : un bon
ajustement ne garantit pas que chaque coefficient ait le sens physique qu'on
lui prête.

**Terme de Wigner.** Pour $$N \approx Z$$, protons et neutrons occupent les mêmes
orbitales et forment des paires proton-neutron. Wigner (1937) a montré qu'une
interaction indépendante du spin et de l'isospin regroupe les états en
« supermultiplets », avec une énergie qui varie en $$\lvert N-Z \rvert$$ et non en
$$(N-Z)^2$$. D'où une **pointe** de liaison en $$N = Z$$, $$-W \lvert N-Z \rvert / A$$. Son
coefficient change beaucoup d'un modèle à l'autre (39 MeV dans M2, 24 MeV
dans M3) : il est porté par une poignée de noyaux légers, où il absorbe aussi
des effets de couches.

---

## 3. Les couches, avec trois nombres (M3)

**La variable de valence.** Pour chaque type de nucléon, situé dans une
couche $$[M_{\text{bas}},\, M_{\text{haut}}[$$ de dégénérescence $$D = M_{\text{haut}} - M_{\text{bas}}$$, avec
$$\nu = n - M_{\text{bas}}$$ nucléons de valence et $$D - \nu$$ trous :

$$
x = \frac{\nu\,(D-\nu)}{D}
$$

C'est le produit du nombre de particules par le nombre de trous, normalisé.
$$x$$ est **symétrique** particule-trou (ν ↔ D − ν), nul sur une couche fermée,
et maximal à mi-couche, où il vaut $$D/4$$. Près d'une fermeture, $$x \approx \nu$$ ou
$$x \approx D - \nu$$ : il compte le plus petit des deux nombres (idée de Casten, forme
de Dieperink et Van Isacker). Exemple : l'U-238 a $$N = 146$$ neutrons dans la
couche $$[126,\, 184[$$, donc $$D = 58$$, $$\nu = 20$$ et $$x_n = 20 \times 38 / 58$$ = 13,10.

**Trois termes.** Avec $$S = x_p + x_n$$ :

$$
B_{\text{couches}} = a_1 S + a_2 S^2 + a_3\, x_p x_n
$$

- $$a_1 S$$ ($$a_1$$ = −1,34 MeV) : un noyau perd de la liaison en s'éloignant des
  couches fermées, ce qui est la définition même d'une couche ;
- $$a_2 S^2$$ (0,021 MeV) : une courbure, qui sature la pénalité en milieu de
  couche ;
- $$a_3\, x_p x_n$$ (0,061 MeV) : l'interaction entre **protons et neutrons de
  valence**. C'est elle qui, en milieu de couche pour les deux types de
  nucléons, rend la forme sphérique instable et **déforme** le noyau.

Les nombres magiques utilisés sont 2, 8, 20, 28, 50, 82, 126 et 184. Le
dernier est une fermeture **prédite**, jamais observée, nécessaire pour les
neutrons au-delà de 126 (actinides).

---

## 4. Ajuster : les moindres carrés

**Le problème.** Chaque modèle M0 à M3 s'écrit $$B \approx X a$$, où $$X$$ est une
matrice $$n \times p$$ (une ligne par noyau, une colonne par terme : $$A$$,
$$-A^{2/3}$$, …) et $$a$$ le vecteur des $$p$$ coefficients. On minimise la somme
des carrés des écarts :

$$
\begin{gathered}
a^* = \underset{a}{\operatorname{argmin}}\ \lVert B - X a \rVert^2 \\[4pt]
\Longleftrightarrow\quad X^\mathsf{T} X\, a^* = X^\mathsf{T} B \qquad \text{(équations normales)}
\end{gathered}
$$

Le problème est **convexe** : la solution est unique, sans point de départ ni
hasard. On ne forme pas $$X^\mathsf{T} X$$, ce qui élèverait au carré les erreurs
d'arrondi : `np.linalg.lstsq` passe par la décomposition en valeurs
singulières de $$X$$. La qualité se mesure par l'écart quadratique moyen,
$$\mathrm{RMS} = \sqrt{\tfrac{1}{n} \sum_i \bigl(B_i - B_i^{\text{modèle}}\bigr)^2}$$.

**Ce que les masses ne savent pas séparer.** Les erreurs formelles des
coefficients s'obtiennent par $$\operatorname{cov}(a) = s^2\, (X^\mathsf{T} X)^{-1}$$, avec
$$s^2 = \sum \text{résidus}^2 / (n-p)$$. Elles supposent des écarts indépendants. Ici,
les écarts sont des erreurs de **modèle**, corrélées sur la carte : ces
erreurs sont donc des bornes inférieures. Mais les **corrélations** entre
coefficients sont instructives :

| | M1 | M2 | M3 |
|---|---|---|---|
| conditionnement de X (colonnes normées) | 135 | 313 | 366 |
| corrélation (surface, Coulomb) | +0,96 | +0,95 | +0,93 |
| corrélation (volume, surface) | +0,99 | +0,26 | +0,26 |
| corrélation (asymétrie, symétrie de surface) | — | +0,95 | +0,94 |
| $$a_S$$ (MeV) | 17,36 ± 0,08 | 18,23 ± 0,14 | 18,42 ± 0,07 |
| $$x$$ de l'U-236 | 0,721 ± 0,002 | 0,727 ± 0,012 | 0,717 ± 0,006 |

Surface et Coulomb sont corrélés à 0,95 : si l'un augmente, l'autre suit, et
les masses changent à peine. Les masses fixent leur **combinaison**, pas
chacun séparément. Or la barrière de fission dépend de leur **rapport**, le
paramètre de fissilité $$x = E_C / 2E_S$$. Nos trois modèles s'accordent sur
$$x \approx 0{,}72$$ à 0,01 près, mais les gouttes publiées étudiées à la
[leçon 6]({{ '/posts/barrieres-de-fission/' | relative_url }}), qui répartissent autrement l'énergie
entre surface, courbure et Coulomb, vont de 0,76 à 0,84. Cet écart, bien
au-delà des erreurs formelles, est une erreur de **modèle**, et il suffit à
doubler la barrière de fission.

---

## 5. M4 : la régression à noyau

Le résidu de M3, $$r_i = B_i - B_i^{\text{M3}}$$, n'est pas du bruit : c'est une
structure lisse sur la carte (Z, N), faite de déformations et de couches
secondaires (carte en bas à gauche de la figure). On l'apprend par
**régression à noyau** (*kernel ridge regression*).

**Le modèle.** Avec $$z_i = (Z_i, N_i)$$ et un noyau gaussien
$$k(z, z') = \exp\bigl(-\lvert z - z' \rvert^2 / 2l^2\bigr)$$, la correction en un point $$z$$ est :

$$
\delta(z) = \sum_i \alpha_i\, k(z, z_i), \qquad (K + \lambda I)\,\alpha = r, \qquad K_{ij} = k(z_i, z_j)
$$

C'est la fonction qui minimise $$\sum_i \bigl(r_i - \delta(z_i)\bigr)^2 + \lambda \lVert \delta \rVert^2$$, la norme étant
celle de l'espace de fonctions engendré par le noyau (théorème du
représentant). C'est aussi la moyenne a posteriori d'un **processus
gaussien** de covariance $$k$$, avec un bruit de variance $$\lambda$$. La matrice
$$K + \lambda I$$ est symétrique définie positive : on résout par une factorisation de
Cholesky. Il y a **un poids $$\alpha_i$$ par noyau** d'entraînement : environ 2 400
paramètres.

**Les hyperparamètres** $$l$$ (longueur de corrélation, en nucléons) et $$\lambda$$
(régularisation) sont choisis par validation croisée à trois plis **sur les
seules données d'entraînement**, parmi $$l$$ ∈ {1,5 ; 2 ; 3 ; 4 ; 6} et
$$\lambda$$ ∈ {0,03 ; 0,3}. On obtient $$l = 3$$ et $$\lambda = 0{,}03$$.

**Pourquoi M4 n'extrapole pas.** Le poids d'un noyau connu à la distance $$d$$
vaut $$\exp(-d^2/2l^2)$$ :

| d (nucléons) | 1 | 2 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| poids, $$l = 3$$ | 0,95 | 0,80 | 0,41 | 0,14 | 0,03 |

Loin des données, $$\delta(z) \to 0$$ et M4 redevient M3 : la correction **s'efface**
exactement là où l'on en aurait besoin. C'est ce que mesure l'épreuve de la
distance, plus bas. C'est la méthode de correction par fonctions de base
radiales utilisée dans la littérature pour améliorer FRDM, HFB et les modèles
WS (Wang et Liu, 2011).

---

## 6. Valider : trois épreuves

Un modèle se juge sur des noyaux qu'il n'a **pas** vus. Trois épreuves, de
difficulté croissante, toutes mesurées par l'écart RMS sur l'énergie de
liaison **totale** B (la métrique standard du domaine), en MeV :

- **Ajustement** : on ajuste sur tous les noyaux et on mesure l'erreur sur
  ces mêmes noyaux.
- **Interpolation** (validation croisée à $$K = 5$$ plis) : on partage les
  noyaux au hasard en cinq groupes $$G_k$$, on ajuste sur quatre et on prédit le
  cinquième, cinq fois. $$\mathrm{RMS}_{\text{CV}} = \sqrt{\tfrac{1}{n} \sum_k \sum_{i \in G_k} \bigl(B_i - B_i^{(-k)}\bigr)^2}$$,
  où $$B^{(-k)}$$ est le modèle ajusté sans le groupe $$k$$.
- **Extrapolation** : on retire les 2 isotopes les plus riches et les 2 plus
  pauvres en neutrons de chaque élément (372 noyaux), puis on les prédit.

| modèle | paramètres | ajustement | interpolation | extrapolation | rapport |
|---|---|---|---|---|---|
| M0 | 5 | 3,76 | 3,76 | 5,47 | 1,46× |
| M1 | 5 | 3,08 | 3,09 | 4,33 | 1,40× |
| M2 | 8 | 2,47 | 2,48 | 2,82 | 1,14× |
| M3 | 11 | 1,21 | 1,22 | 1,42 | 1,16× |
| M4 | ~2400 | **0,23** | **0,28** | **0,49** | **1,71×** |

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

## 7. Ce que disent les coefficients

### La goutte liquide réajustée (M1)

| terme | manuel | réajusté |
|---|---|---|
| volume | 15,750 | 15,570 |
| surface | 17,800 | 17,360 |
| Coulomb | 0,711 | 0,705 |
| asymétrie | 23,700 | 22,859 |
| appariement | 11,180 | 12,121 |

Les coefficients bougent peu : ceux du manuel étaient déjà ajustés, sur des
tables plus anciennes et avec $$Z^2$$ au lieu de $$Z(Z-1)$$. Le gain M0 → M1 vient
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

$$
a_{\text{sym}}(A) = J - Q\, A^{-1/3}
$$

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

## 8. L'erreur en fonction de la distance au connu

L'épreuve qui compte vraiment. Pour chaque élément ayant au moins 14 isotopes
mesurés, on retire les **8 plus riches en neutrons**, on ajuste sur le reste,
puis on prédit les noyaux retirés. $$d$$ = nombre de neutrons au-delà du dernier
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

## 9. Où se situent ces modèles

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

- Chaque terme a une origine physique qu'on peut chiffrer : saturation
  (volume), tension de surface (σ = 0,92 MeV/fm²), sphère chargée
  (r₀ = 1,225 fm), principe de Pauli (ε_F/3 = 12,3 MeV d'asymétrie
  cinétique), appariement.
- Chaque pièce de physique paie : de la goutte liquide de manuel (M0) au
  modèle avec couches (M3), l'erreur sur $$B$$ passe de 3,8 à 1,2 MeV avec
  seulement 11 paramètres.
- La **symétrie de surface** traduit la **peau de neutrons** : l'excès de
  neutrons coûte moins cher en surface qu'au cœur.
- Le terme de couches est nul pour les noyaux doublement magiques et coûte
  plus de 10 MeV en milieu de couche ; l'interaction proton-neutron de valence
  pousse ces noyaux à se **déformer**.
- Les moindres carrés fixent bien les **combinaisons** de coefficients que
  les masses voient, mal les autres : surface et Coulomb sont corrélés à
  0,95, et leur rapport, le paramètre de fissilité, reste incertain.
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

## Bibliographie

{% include bibliographie.html %}

{% include cours-pied.html %}
