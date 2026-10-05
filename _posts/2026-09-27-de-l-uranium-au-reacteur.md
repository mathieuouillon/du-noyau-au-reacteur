---
title: "De l'uranium au réacteur"
date: 2026-09-27 09:00:00 +0200
categories: [Physique nucléaire, Fission]
tags: [fission, uranium, modération, quatre facteurs]
description: "Pourquoi l'uranium 235 fissionne, combien d'énergie sort, et pourquoi un REP ne peut pas fonctionner à l'uranium naturel."
image:
  path: /assets/img/nucleaire/fission_bilan.png
  alt: "De l'uranium au réacteur"
---

Document accompagnant [`fission.py`]({{ '/assets/code/fission.py' | relative_url }}) et [`etude_fission.py`]({{ '/assets/code/etude_fission.py' | relative_url }}).
Tous les chiffres sont produits par le code.

C'est le point de départ logique de la série : [`fission.py`]({{ '/assets/code/fission.py' | relative_url }}) explique d'où
viennent les sections efficaces, [`diffusion.py`]({{ '/assets/code/diffusion.py' | relative_url }}) résout le transport,
[`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }}) conçoit le cœur.


![Courbe de liaison, ralentissement, courbe de modération et quatre facteurs]({{ '/assets/img/nucleaire/fission_bilan.png' | relative_url }})
_Courbe de liaison, ralentissement, courbe de modération et quatre facteurs._

---

## 1. Pourquoi un noyau lourd libère de l'énergie en se cassant

Une seule courbe explique la fission **et** la fusion : l'énergie de liaison
par nucléon `B/A`. Elle monte de l'hydrogène jusqu'à un maximum vers A ≈ 60
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

**Le terme coupable** est le terme coulombien, en `Z²/A^(1/3)`. Les protons
se repoussent **tous** mutuellement (en `Z²`), alors que la force nucléaire
ne lie qu'aux voisins immédiats (en `A`). Passé le fer, la répulsion
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
cinétique — seuil ~1 MeV. C'est le facteur `ε ≈ 1,03` : 3 % des fissions
d'un REP ont lieu sur l'U-238, par des neutrons encore rapides.

## 3. Combien d'énergie

`Q = B(fragment 1) + B(fragment 2) − B(U-236)`. Les neutrons libres ont une
énergie de liaison nulle et ne comptent pas.

| partition de U-236 | Q (goutte liquide) |
|---|---|
| 118/46 + 118/46 (symétrique) | **185,0 MeV** |
| 133/52 + 100/40 + 3n | 164,6 |
| 141/56 + 92/36 + 3n | 156,4 |

**Pourquoi la fission réelle est-elle asymétrique ?** Le calcul donne un
maximum pour la partition **symétrique**, alors que la fission réelle fait
deux bosses vers A=95 et A=139.

Ce n'est **pas** une erreur du modèle — les masses mesurées donnent elles
aussi plus d'énergie au partage symétrique (193 MeV contre 167). L'asymétrie
n'est donc pas un effet de bilan énergétique : elle vient de la **dynamique**
au point de scission, où la surface d'énergie potentielle favorise des
fragments proches des couches fermées Z=50 et N=82. Le noyau ne choisit pas le
partage le plus exothermique, il suit le chemin le plus facile.

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

## 4. Ralentir les neutrons : le choix qui décide de tout

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

## 5. Les quatre facteurs

```
k_inf = η · ε · p · f
```

Le cycle du neutron, lu à l'envers : `f` il est absorbé dans le combustible
plutôt qu'ailleurs · `p` il a survécu aux résonances · `ε` quelques fissions
rapides l'ont précédé · `η` son absorption produit η nouveaux neutrons.

`η` ne dépend **que** de l'enrichissement, jamais de la géométrie. C'est le
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

## 6. La sous-modération : le choix de sûreté fondamental

Ajouter de l'eau fait monter `p` (on ralentit mieux) et baisser `f` (l'eau
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
Vm/Vf diminue, et comme on est à gauche du maximum, `k_inf` **diminue**. La
puissance baisse d'elle-même. Le coefficient de température modérateur est
négatif, le réacteur est intrinsèquement stable.

À droite de l'optimum ce serait l'inverse : dilatation → plus de réactivité →
plus de puissance → plus de dilatation. Emballement. Un cœur sur-modéré serait
interdit.

C'est aussi la vraie raison de la limite en bore vue dans
[Neutronique du cœur : comprendre et concevoir]({{ '/posts/neutronique-du-coeur/' | relative_url }}) : le bore étant dans l'eau, trop de bore rend ce
coefficient positif.

## 7. Pourquoi un REP ne peut pas marcher à l'uranium naturel

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

## 8. Bouclage : d'où viennent vraiment les constantes de [`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }})

[`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }}) utilisait `Σa2 = 0,082` et `νΣf2 = 0,1375` sans les justifier.
Tentative naïve, `Σ = N·σ` avec σ à 2200 m/s :

| | naïf | corrigé Maxwell | [`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }}) |
|---|---|---|---|
| Σa2 | 0,178 | 0,112 | 0,082 |
| νΣf2 | 0,327 | 0,206 | 0,1375 |

Facteur 2 d'écart au départ. Deux étapes manquaient :

**(a) Moyenne sur le spectre.** La valeur à 2200 m/s est une *convention*,
pas une moyenne. Les neutrons thermiques suivent une maxwellienne à la
température du modérateur ; pour une section en 1/v, la moyenne vaut
`σ(2200)·(√π/2)·√(293,6/T)`, soit **0,63** à 580 K.

**(b) Facteur de désavantage.** L'homogénéisation suppose le flux uniforme.
Il ne l'est pas : le crayon absorbe, donc le flux thermique y est déprimé.
Les sections homogénéisées correctes sont pondérées par le **flux**, pas par
le volume. Ce rapport vaut typiquement 1,1 à 1,3 dans un REP.

Après correction spectrale on est dans un facteur 1,3 des valeurs réelles,
avec un calcul de coin de table. C'est l'ordre de grandeur attendu — et c'est
exactement pourquoi la génération des constantes de groupe est un métier à
part entière (APOLLO, CASMO, WIMS) et non une multiplication `N·σ`.

**La chaîne complète :**

```
données nucléaires évaluées (JEFF, ENDF/B)
  → traitement des résonances et autoprotection
  → calcul de réseau 2D fin (assemblage, multigroupe)
  → condensation à 2 groupes + homogénéisation pondérée par le flux
  → constantes utilisées par coeur.py
```

Chaque flèche est un domaine de recherche.

---

## Références

- **P. Reuss, _Précis de neutronique_** (EDP Sciences) — la meilleure entrée
  en matière en français ; couvre exactement ce document.
- **J. Bussac & P. Reuss, _Traité de neutronique_** — la référence complète.
- **Duderstadt & Hamilton, _Nuclear Reactor Analysis_**, ch. 2-3 et 10 — pour
  les quatre facteurs et l'auto-protection.
- **Table of Nuclides** (KAERI, NNDC) — pour vérifier les énergies de liaison
  et sections efficaces citées ici.
