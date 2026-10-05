---
title: "Les barrières de fission : de la goutte liquide à la double bosse"
date: 2026-09-30 21:00:00 +0200
categories: ["II · Couches et formes des noyaux", "Barrières de fission"]
tags: [fission, barrière, point selle, isomère de fission, Strutinsky, goutte liquide]
description: "Suivre la forme d'un noyau jusqu'à la scission : barrière de la goutte liquide, formules de Bohr et Wheeler, puis la double bosse et les isomères de fission créés par les couches."
image:
  path: /assets/img/nucleaire/barrieres_fission.png
  alt: "Les barrières de fission"
lecon: 6
partie: "II"
objectifs:
  - "Définir une barrière de fission comme le col d'un paysage d'énergie potentielle."
  - "Retrouver la stabilité de la sphère et la barrière de Bohr et Wheeler, (98/135)(1−x)³ E_s."
  - "Calculer exactement la surface et l'énergie coulombienne d'une forme à col, jusqu'à la scission."
  - "Comprendre pourquoi deux gouttes ajustées sur les mêmes masses donnent des barrières qui varient du simple au double."
  - "Expliquer la double bosse des actinides et les isomères de fission par la correction de couches."
prerequis: [2, 4]
code: [barrieres.py, etude_barrieres.py, nilsson.py, strutinsky.py, frdm.py]
sources: >-
  Tous les chiffres de cette leçon sont produits par `etude_barrieres.py` (environ une minute).
  Paramètres de la goutte : Pomorski et Dudek (2003), table I et éq. (3.5). Barrières de
  référence : RIPL-3, via la table 1 de Ryssens et al. (2023) ; énergies des isomères :
  Samyn et al. (2004), via la même table. Partie microscopique : modèle de la leçon 4.
bibliographie:
  - cle: bohr1939
    note: "L'article fondateur : paramètre de fissilité, stabilité de la sphère, barrière de la goutte."
  - cle: cohen1963
    note: "Les formes et énergies du point selle de la goutte chargée, calculées sans développement."
  - cle: brack1972
    note: "La revue « Funny Hills » : paramétrisation (c, h, α), méthode de Strutinsky, double bosse."
  - cle: strutinsky1967
    note: "La correction de couches, et l'explication de la double bosse."
  - cle: polikanov1962
    note: "La découverte du premier isomère de fission (Am-242, 14 ms)."
  - cle: bjornholm1980
    note: "La revue de référence sur la double bosse, les isomères et les résonances de fission."
  - cle: hill1953
    note: "La pénétrabilité d'une barrière parabolique."
  - cle: myers1966
    note: "Le jeu de paramètres MS-LD."
  - cle: pomorski2003
    note: "Le jeu LSD, et la comparaison des barrières macroscopiques aux barrières mesurées."
  - cle: moller2001
    note: "Le paysage de la fission dans un espace de formes à cinq dimensions."
  - cle: moller2009
    note: "Barrières de fission calculées pour les noyaux lourds."
  - cle: scamps2018
    note: "L'origine de la fission asymétrique : couches des fragments déformés en poire."
  - cle: capote2009
    note: "RIPL-3 : barrières de fission recommandées pour l'évaluation des données nucléaires."
  - cle: ryssens2023
    note: "Table de valeurs de référence (RIPL-3) et calculs HFB récents des barrières."
  - cle: samyn2004
    note: "Source des énergies d'isomères utilisées ici."
  - cle: krappe2012
    note: "Manuel moderne de théorie de la fission."
  - cle: schunck2016
    note: "Revue des approches microscopiques de la fission."
  - cle: vandenbosch1973
    note: "Le manuel classique sur la fission."
---

{% include cours-entete.html %}

La [leçon 2]({{ '/posts/la-fission-de-l-uranium/' | relative_url }}) a affirmé que l'U-235 fissionne avec des
neutrons lents parce que la capture apporte au noyau composé U-236 plus que sa
**barrière de fission**, d'environ 6 MeV. La [leçon 4]({{ '/posts/vers-frdm/' | relative_url }}) a montré
pourquoi les noyaux se déforment, mais avec des ellipsoïdes seulement : impossible
d'y former un col, donc d'aller jusqu'à la scission.

Cette leçon calcule enfin cette barrière. D'abord avec la goutte liquide seule,
sur une famille de formes qui va de la sphère jusqu'à deux fragments, puis avec
les couches, qui creusent dans la barrière un second puits : c'est la **double
bosse** des actinides.

![Barrières de fission : formes, surface d'énergie, barrière de la goutte, double bosse]({{ '/assets/img/nucleaire/barrieres_fission.png' | relative_url }})
_Figure produite par `etude_barrieres.py`. En haut : la goutte liquide seule (U-236). En bas : la loi de Bohr et Wheeler, puis l'effet des couches._

---

## 1. Le paysage d'énergie potentielle

On décrit la forme du noyau par quelques paramètres `q = (q₁, q₂, …)`
(élongation, col, asymétrie) et on trace son énergie `E(q)`. Ce **paysage
d'énergie potentielle** contient tout :

- l'**état fondamental** est le fond d'un puits ;
- la **scission** est la région où le col se rompt et où le noyau devient deux
  fragments ;
- entre les deux, le noyau doit passer par un **col** au sens géographique,
  un **point selle** : minimum dans toutes les directions sauf une, maximum le
  long du chemin de fission.

La **barrière de fission** est la hauteur de ce col au-dessus du fond :

```
B_f = E(point selle) − E(état fondamental)
```

Trouver un point selle est plus subtil que trouver un minimum. Sur une surface à
deux dimensions `E(c, h)`, on utilise la propriété suivante : tout chemin qui va
du fond à la scission traverse chaque droite `c = cte`, donc y passe au moins par
l'énergie `min_h E(c, h)`. D'où une borne inférieure :

```
B_f ≥ max_c [ min_h E(c, h) ]
```

Elle est **exacte** quand la vallée `h*(c)` est continue, ce qui est le cas des
actinides. Elle n'est plus qu'une borne quand deux vallées se croisent (noyaux
plus légers, dont le col est proche de la scission) : on s'en tiendra donc aux
actinides.

---

## 2. La goutte de Bohr et Wheeler

**La stabilité de la sphère.** Déformons une sphère de rayon `R₀` à volume
constant, `R(θ) = R₀ [1 + Σ αₙ Pₙ(cos θ)]` (à une constante près qui conserve
le volume). Au second ordre, Bohr et Wheeler (1939) obtiennent :

```
E_s = E_s⁰ [ 1 + Σₙ (n−1)(n+2) / (2(2n+1)) · αₙ² ]
E_c = E_c⁰ [ 1 − Σₙ 5(n−1) / (2n+1)²    · αₙ² ]
```

Pour la déformation quadrupolaire (`n = 2`), avec `E_c⁰ = 2x E_s⁰` :

```
ΔE = E_s⁰ [ (2/5) − (2x)(1/5) ] α₂²  =  (2/5) (1 − x) E_s⁰ α₂²
```

C'est la définition même du **paramètre de fissilité** `x = E_c⁰ / 2E_s⁰` : la
sphère est stable tant que `x < 1`. Pour `n = 4`, les coefficients valent 1 et
5/27 : la sphère est stable contre toutes les déformations si elle l'est contre
la déformation quadrupolaire.

**La barrière.** En poussant le développement au troisième ordre en α₂ et en
laissant α₄ s'ajuster, Bohr et Wheeler obtiennent la barrière pour `x` proche
de 1 :

```
B_f ≈ (98/135) (1 − x)³ E_s⁰  ≈  0,726 (1 − x)³ E_s⁰
```

La barrière s'effondre comme le **cube** de `1 − x`. Pour l'uranium, `1 − x`
vaut 0,15 à 0,25 : la barrière ne représente que 0,3 à 1 % de l'énergie de
surface, qui fait plusieurs centaines de MeV. D'où une barrière de quelques MeV, et une extrême
sensibilité à `x`. On vérifie ce résultat numériquement à la section 4.

---

## 3. Une famille de formes jusqu'à la scission

**Les formes.** On prend la famille à trois paramètres de Brack *et al.* (1972),
en coordonnées cylindriques, longueurs en unités de `R₀` :

```
ρ²(z) = (c² − z²) · (A + B z²/c² + α z/c),      −c ≤ z ≤ c
B = 2h + (c − 1)/2
```

- `c` est l'**élongation** (demi-longueur ; 1 pour la sphère) ;
- `h` contrôle le **col** : plus `h` est grand, plus la taille se creuse ;
- `α` introduit une **asymétrie** gauche-droite, c'est-à-dire des fragments
  de masses différentes ;
- `A` est fixé par la **conservation du volume**. Comme
  `∫(c² − z²) dz = 4c³/3` et `∫(c² − z²) z² dz = 4c⁵/15`, on a
  `V = (4π/3) c³ (A + B/5)` ; `V = 4π/3` impose `A = 1/c³ − B/5`. Le terme en
  `α`, impair en `z`, ne change pas le volume.

`c = 1, h = 0` donne la sphère ; `h = −(c − 1)/4` annule `B` et donne un
sphéroïde de demi-axes `c` et `c^(−1/2)` ; la forme se **rompt** quand
`ρ²` s'annule à l'intérieur de `]−c, c[`.

**Les trois énergies de forme.** On écrit `f(z) = ρ²(z)`. Tout s'exprime avec
`f`, `f'` et `f''`, ce qui évite les singularités aux pointes, où `ρ → 0` mais
`ρρ' = f'/2` reste fini :

```
surface   :  S = 2π ∫ √(f + f'²/4) dz                                       B_s = S / 4π
courbure  :  ∮ (k₁ + k₂) dS = 2π ∫ [ 1 − (2 f f'' − f'²) / (4f + f'²) ] dz   B_k = … / 8π
```

**L'énergie coulombienne** d'une charge uniforme de densité `ρ_q` est une
intégrale de volume à six dimensions,
`E_C = (ρ_q²/2) ∫∫ d³r d³r' / |r − r'|`. Elle se ramène à une intégrale de
**surface** grâce à l'identité `∇²|r − r'| = 2/|r − r'|` et à deux applications
du théorème de Gauss :

```
∫_V d³r / |r − r'| = ½ ∮ dS · (r − r')/|r − r'|
∫_V d³r' (r − r')/|r − r'| = − ∮ dS' |r − r'|

⇒   E_C = − (ρ_q²/4) ∮∮ (dS · dS') |r − r'|
```

Avec la symétrie axiale, `dS = ρ (cos φ, sin φ, −ρ') dz dφ`, et l'intégrale
sur l'un des deux angles donne 2π. Il reste trois intégrales (`z`, `z'` et
l'angle relatif `ψ`), dont l'intégrande est **continu** : plus de singularité
en `1/|r − r'|`. Le changement de variable `z = c sin t` absorbe le
comportement en racine carrée aux pointes, et une quadrature de Gauss converge
très vite.

**Vérifications.**

| test | calcul | attendu |
|---|---|---|
| sphère : `B_s`, `B_k`, `B_c` | 1,000000 ; 1,000000 ; 0,999998 | 1 |
| sphéroïde `c = 1,5` : `B_s` | 1,059829 | 1,059829 (formule exacte) |
| sphéroïde `c = 1,5` : `B_c` | 0,967453 | 0,967451 (formule exacte) |
| sphéroïde `c = 2,0` : `B_c` | 0,908721 | 0,908711 |
| `α₂ = 0,02` : `(B_s − 1)/α₂²` | 0,390 | 2/5 au second ordre |
| `α₂ = 0,02` : `(B_c − 1)/α₂²` | −0,200 | −1/5 |

La surface est exacte à la précision machine ; l'énergie coulombienne à
quelques 10⁻⁶ près, soit **un centième de MeV** sur les ~1 000 MeV de
l'uranium. C'est largement suffisant pour des barrières de quelques MeV.

---

## 4. La barrière de la goutte

**L'énergie macroscopique** d'une forme s'écrit, à partir de la sphère :

```
E_def = E_s (B_s − 1) + E_k (B_k − 1) + E_c (B_c − 1)
E_s = b_s (1 − κ_s I²) A^(2/3),   E_k = b_k (1 − κ_k I²) A^(1/3),
E_c = (3/5) e² Z² / (r₀ A^(1/3)),   I = (N − Z)/A
```

On compare deux jeux de paramètres **publiés et ajustés sur les masses** :

| jeu | `b_s` (MeV) | `κ_s` | `b_k` (MeV) | `κ_k` | `r₀` (fm) |
|---|---|---|---|---|---|
| MS-LD (Myers et Swiatecki, 1966) | 18,56 | 1,79 | — | — | 1,2049 |
| LSD (Pomorski et Dudek, 2003) | 16,9707 | 2,2938 | 3,8602 | −2,3764 | 1,21725 |

**La surface d'énergie de l'U-236** (panneau b) a la forme attendue : un
fond près de la sphère, une vallée de fission qui monte doucement le long de
`h ≈ 0`, un col en `c = 1,52`, `h = −0,02` (l'étoile), puis la descente vers
la scission. La forme du point selle (panneau a, au milieu) est allongée, à
peine creusée : un ballon de rugby aux flancs droits. Le col se creuse ensuite
rapidement.

**Bohr et Wheeler, vérifiés.** Pour une goutte sans courbure, la barrière
réduite `ξ(x) = B_f / E_s` ne dépend que de `x` :

| x | ξ calculé | (98/135)(1 − x)³ | rapport | `c` au point selle |
|---|---|---|---|---|
| 0,75 | 0,01243 | 0,01134 | 1,096 | 1,673 |
| 0,80 | 0,00605 | 0,00581 | 1,041 | 1,507 |
| 0,85 | 0,00247 | 0,00245 | 1,008 | 1,362 |
| 0,90 | 0,00072 | 0,00073 | 0,991 | 1,235 |
| 0,95 | 0,00009 | 0,00009 | 1,008 | 1,117 |

Le calcul rejoint le développement quand `x → 1`, où il est exact, et s'en
écarte de 10 % à `x = 0,75`, où le point selle est déjà loin de la sphère.
Deux méthodes indépendantes donnent le même résultat : c'est le meilleur test
du code.

**La même goutte, deux barrières.**

| noyau | x (MS-LD) | B_f (MS-LD) | x (LSD) | B_f (LSD) | c au selle (LSD) |
|---|---|---|---|---|---|
| U-234 | 0,761 | 6,94 | 0,845 | 3,68 | 1,519 |
| U-236 | 0,759 | 7,15 | 0,844 | 3,76 | 1,524 |
| U-238 | 0,757 | 7,35 | 0,843 | 3,82 | 1,528 |
| Pu-240 | 0,776 | 5,66 | 0,863 | 2,78 | 1,464 |
| Pu-242 | 0,774 | 5,83 | 0,862 | 2,84 | 1,468 |
| Cm-244 | 0,794 | 4,41 | 0,882 | 1,99 | 1,405 |
| Cf-252 | 0,808 | 3,61 | 0,899 | 1,44 | 1,359 |

(MeV)

Deux gouttes qui décrivent les masses avec une précision comparable donnent des
barrières qui varient **du simple au double**. La raison tient en une ligne :
`B_f ∝ (1 − x)³`, et les deux jeux diffèrent sur `x` de 0,08. Ils répartissent
différemment l'énergie entre volume, surface, courbure et asymétrie de surface.
Les masses contraignent mal cette répartition : Pomorski et Dudek montrent que
la qualité de l'ajustement varie à peine quand on échange de l'énergie entre
surface et courbure. La déformation, elle, la révèle.

C'est la leçon de Pomorski et Dudek : les masses seules fixent mal le terme de
surface. Une fois les énergies microscopiques de l'état fondamental ajoutées,
MS-LD surestime les barrières des noyaux lourds de 3 à 4 MeV, alors que LSD,
qui possède un terme de **courbure** en `A^(1/3)`, les reproduit à 0,88 MeV
près (écart RMS, noyaux Z > 70), sans avoir été ajustée sur elles. Les
barrières de fission sont un **test** des modèles de masse, pas une
conséquence automatique.

> **Retour sur la leçon 4.** Les coefficients ajustés à la leçon 4 donnent
> `x = 0,709` pour U-236, encore plus bas que MS-LD : avec des ellipsoïdes et
> un ajustement sur les seules masses, ce modèle surestimerait la barrière.
> C'est la même difficulté, et une raison supplémentaire de ne pas lui
> demander de barrières.
{: .prompt-warning }

---

## 5. Symétrique ou asymétrique ?

Au point selle de l'U-236, on fait varier l'asymétrie `α` :

| α | 0 | 0,05 | 0,10 | 0,15 | 0,20 |
|---|---|---|---|---|---|
| E_def (MeV) | 3,76 | 4,06 | 4,95 | 6,41 | 8,43 |

L'énergie **croît** avec l'asymétrie : la goutte est stable contre elle. La
goutte liquide prédit donc une fission **symétrique**, alors que l'uranium se
casse en fragments inégaux (`A ≈ 95` et `A ≈ 139`, [leçon 2]({{ '/posts/la-fission-de-l-uranium/' | relative_url }})).

Le désaccord n'est pas une affaire de bilan d'énergie : la leçon 2 a montré que
le partage symétrique libère **plus** d'énergie. Il vient de la forme du
chemin. Au-delà de la barrière interne, les **couches des fragments naissants**
abaissent les formes asymétriques. Les calculs dans des espaces de formes à
cinq dimensions (Möller *et al.*, 2001) reproduisent ainsi l'asymétrie des
actinides. Une analyse plus récente (Scamps et Simenel, 2018) l'attribue
surtout aux couches du fragment lourd **déformé en poire**, autour de
Z = 52 à 56. La goutte seule ne peut rien en dire : c'est un effet quantique.

---

## 6. Les couches : la barrière à double bosse

**L'idée de Strutinsky.** La correction de couches de la leçon 4 n'est pas
une propriété de la sphère seulement : elle dépend de la **forme**. Elle est
positive quand le niveau de Fermi tombe dans une région dense en niveaux,
négative quand il tombe dans un trou du spectre. Or, en se déformant, le
spectre se réorganise et de nouveaux trous apparaissent. Ajoutée à la goutte
qui monte doucement, une correction qui **oscille** avec la déformation
découpe la barrière en deux bosses séparées par un second puits.

**Le calcul.** Même modèle qu'à la leçon 4 : `E(δ) = E_goutte(δ) + δE_couches(δ)
+ δE_appariement(δ)`, sur des **sphéroïdes**, avec la goutte LSD. La correction
microscopique de l'U-236 (panneau e) :

| δ | 0 | 0,10 | 0,20 | 0,30 | 0,35 | 0,40 | 0,45 | 0,50 | 0,55 |
|---|---|---|---|---|---|---|---|---|---|
| couches + appariement (MeV) | +3,5 | +2,2 | −0,4 | +0,9 | +1,5 | +1,2 | −1,6 | −3,5 | −2,0 |

Trois régions se dessinent : la sphère est **pénalisée** (+3,5 MeV, l'U-236
est en milieu de couche), un **premier puits** à `δ ≈ 0,2`, l'état
fondamental déformé de la leçon 4, une **bosse** vers `δ ≈ 0,35–0,40`, puis un
**second puits** profond à `δ ≈ 0,5`.

**Pourquoi δ = 0,5 ?** Pour un oscillateur à symétrie axiale,
`ω⊥² = ω₀²(1 + 2δ/3)` et `ω_z² = ω₀²(1 − 4δ/3)`. À `δ = 0,5`,
`ω⊥/ω_z = √(4/3 ÷ 1/3) = 2` : les fréquences sont dans un **rapport
rationnel**, les niveaux se regroupent en couches très dégénérées, et de
grands trous s'ouvrent dans le spectre, comme pour la sphère. C'est la
**superdéformation**, avec un rapport d'axes de 2 pour 1. Le spin-orbite
déplace un peu ces couches, mais le mécanisme est celui-là.

**Convergence.** Aux grandes déformations, il faut beaucoup de couches
d'oscillateur dans la base :

| N_max | δ = 0,20 | 0,40 | 0,50 | 0,55 | 0,60 |
|---|---|---|---|---|---|
| 14 | −0,41 | 1,21 | −4,05 | −2,68 | 0,17 |
| 18 | −0,42 | 1,14 | −3,48 | −1,90 | 0,15 |
| 20 | −0,43 | 1,18 | −3,48 | −1,96 | 0,47 |

Avec `N_max = 20`, le résultat est stable à 0,1 MeV près jusqu'à `δ = 0,55`.
On s'arrête à 0,60.

**Résultats, face aux valeurs de référence** (panneau f) :

| noyau | δ fondamental | barrière interne | isomère : énergie | isomère : δ | réf. E_I | réf. E_II | réf. isomère |
|---|---|---|---|---|---|---|---|
| U-234 | 0,175 | 3,75 | 1,14 | 0,475 | 5,50 | 4,80 | — |
| U-236 | 0,200 | 4,22 | 1,46 | 0,475 | 5,67 | 5,00 | 2,30 |
| U-238 | 0,200 | 4,79 | 2,47 | 0,475 | 6,30 | 5,50 | 2,60 |
| Pu-240 | 0,200 | 5,07 | 2,48 | 0,475 | 6,05 | 5,45 | 2,25 |
| Pu-242 | 0,200 | 5,63 | 3,36 | 0,500 | 5,85 | 5,05 | — |
| Cm-244 | 0,225 | 6,08 | 3,42 | 0,500 | 6,18 | 5,10 | 1,04 |

(MeV, depuis l'état fondamental du modèle.) Les références sont les barrières
**primaire** (E_I, la plus haute) et **secondaire** (E_II, la plus basse)
recommandées par RIPL-3. Ce sont des quantités extraites des sections efficaces
de fission à l'aide d'un modèle, pas des mesures directes.

Ce qui est **juste** : la structure à deux bosses, la position des deux puits,
et l'ordre de grandeur de la barrière interne (4 à 6 MeV) comme de l'énergie
de l'isomère (1 à 3,5 MeV).

Ce qui ne l'est **pas** :

- les barrières internes de l'uranium sont plus basses que les deux barrières
  de référence, de 0,7 à 1,8 MeV selon le noyau ;
- la barrière du modèle croît régulièrement de U-234 à Cm-244, alors que les
  références restent entre 5,5 et 6,3 MeV ;
- l'isomère du Cm-244 est 2,4 MeV trop haut.

On ne pouvait pas espérer mieux avec un seul paramètre de forme, sans
hexadécapole ni triaxialité.

**Les isomères de fission.** Un noyau piégé dans le second puits est un état
excité à longue durée de vie, qui fissionne en traversant **une seule** bosse
au lieu de deux. Sa période de fission spontanée est donc énormément plus
courte que celle de l'état fondamental. C'est ce qu'ont observé Polikanov *et
al.* en 1962, avec un état de l'Am-242 qui fissionne en 14 ms. Il a fallu
attendre la méthode de Strutinsky (1967) pour comprendre qu'il s'agissait
d'un noyau superdéformé.

---

## 7. Où s'arrête ce modèle

La goutte le long des sphéroïdes et la goutte dans la vallée de fission ne
coûtent pas la même chose. À moment quadrupolaire `Q₂` égal (U-236, LSD) :

| δ | Q₂ (unités réduites) | E sphéroïde | E vallée (c, h) | écart |
|---|---|---|---|---|
| 0,30 | 0,317 | 1,85 | 1,82 | 0,03 |
| 0,40 | 0,492 | 3,43 | 2,86 | 0,57 |
| 0,50 | 0,756 | 5,88 | 3,69 | 2,19 |
| 0,60 | 1,255 | 10,51 | 2,95 | 7,56 |

(MeV)

Jusqu'à `δ ≈ 0,4`, l'ellipsoïde est une bonne approximation : la barrière
interne est fiable, à la physique du modèle près. Au-delà, un noyau qui peut
creuser un col coûte beaucoup moins cher. Le second puits est donc
probablement **trop peu profond** dans notre calcul, et la **barrière externe**
n'y est pas du tout : il faudrait un potentiel microscopique défini sur des
formes à col, avec asymétrie, ce que fait FRDM avec son potentiel de Yukawa
replié dans un espace de formes à cinq dimensions. On sait aussi que la
**triaxialité** abaisse la barrière interne des actinides, et que l'asymétrie
abaisse la barrière externe.

---

## 8. Traverser la barrière

Une barrière n'est pas un mur : on la traverse par effet tunnel. Pour une
barrière parabolique de hauteur `B_f` et de courbure `ħω`, Hill et Wheeler
(1953) donnent la probabilité de passage :

```
T(E) = 1 / (1 + exp(2π (B_f − E) / ħω))
```

| E − B_f (MeV) | −2 | −1 | −0,5 | 0 | +0,5 | +1 |
|---|---|---|---|---|---|---|
| ħω = 0,5 MeV | 1,2·10⁻¹¹ | 3,5·10⁻⁶ | 1,9·10⁻³ | 0,5 | ≈ 1 | ≈ 1 |
| ħω = 1,0 MeV | 3,5·10⁻⁶ | 1,9·10⁻³ | 0,041 | 0,5 | 0,96 | ≈ 1 |

Deux conséquences :

- le **seuil** de fission n'est pas une marche : la section efficace de
  fission de l'U-238 monte progressivement autour de son seuil, au lieu de
  sauter de zéro à sa valeur maximale ;
- la **fission spontanée** est un effet tunnel à travers toute la barrière,
  d'où des périodes qui s'étalent sur des dizaines d'ordres de grandeur selon la
  hauteur et la largeur de la barrière.

C'est cette formule, appliquée à une barrière à deux bosses, que les
évaluateurs ajustent pour reproduire les sections efficaces de fission. Les
« barrières recommandées » de RIPL-3 en sont le résultat.

---

## À retenir

- La barrière de fission est la hauteur du **point selle** du paysage
  d'énergie potentielle au-dessus de l'état fondamental.
- Pour une goutte, `ΔE = (2/5)(1 − x) E_s⁰ α₂²` : la sphère est stable tant que
  `x < 1`, et la barrière vaut `(98/135)(1 − x)³ E_s⁰` près de `x = 1`. Le
  calcul numérique le confirme.
- Une forme à col (c, h, α), avec l'énergie coulombienne réduite à une
  intégrale de surface, suit le noyau jusqu'à la scission. Pour l'U-236, la
  goutte LSD place le col en `c = 1,52` et donne une barrière de 3,8 MeV.
- La barrière dépend de `x` au cube : deux gouttes ajustées sur les masses
  donnent 3,8 ou 7,2 MeV. Les barrières **testent** les modèles de masse.
- La goutte seule préfère la fission **symétrique** ; l'asymétrie vient des
  couches des fragments.
- La correction de couches, qui oscille avec la déformation, crée la **double
  bosse** et un second puits superdéformé (`ω⊥/ω_z = 2`, `δ ≈ 0,5`) : les
  **isomères de fission**.

---

## Bibliographie

{% include bibliographie.html %}

{% include cours-pied.html %}
