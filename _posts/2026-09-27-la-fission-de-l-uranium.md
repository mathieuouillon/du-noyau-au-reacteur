---
title: "La fission de l'uranium : énergie, fissilité, fragments"
date: 2026-09-27 09:00:00 +0200
categories: ["I · Énergie du noyau et fission", "Fission"]
tags: [fission, uranium, appariement, Q de fission]
description: "Pourquoi un noyau lourd libère de l'énergie en se cassant, pourquoi l'U-235 fissionne avec des neutrons lents et pas l'U-238, et où passent les 200 MeV."
image:
  path: /assets/img/nucleaire/fission_bilan.png
  alt: "La fission de l'uranium"
redirect_from:
  - /posts/de-l-uranium-au-reacteur/
lecon: 2
partie: "I"
objectifs:
  - "Expliquer, à partir de la courbe B/A, pourquoi un noyau lourd libère de l'énergie en se cassant."
  - "Distinguer noyau fissile et fissionnable, et relier la fissilité à l'appariement (S_n face à la barrière)."
  - "Calculer le Q d'une fission et détailler le bilan des ~200 MeV."
  - "Comprendre pourquoi la fission réelle est asymétrique alors que le bilan d'énergie favorise la symétrie."
prerequis: [1]
code: [fission.py, etude_fission.py]
---

{% include cours-entete.html %}

Les sections 1 à 3 de [`fission.py`]({{ '/assets/code/fission.py' | relative_url }}) et de
[`etude_fission.py`]({{ '/assets/code/etude_fission.py' | relative_url }}). Tous les chiffres sont produits par le code.
La suite de ces deux scripts, qui ralentit les neutrons et assemble un réseau de
réacteur, fait l'objet de la [leçon 6]({{ '/posts/ralentir-les-neutrons/' | relative_url }}).

![Courbe de liaison, ralentissement, courbe de modération et quatre facteurs]({{ '/assets/img/nucleaire/fission_bilan.png' | relative_url }})
_Figure produite par `etude_fission.py`. Le panneau en haut à gauche (courbe de liaison) concerne cette leçon ; les trois autres, la [leçon 6]({{ '/posts/ralentir-les-neutrons/' | relative_url }})._

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
cinétique — seuil ~1 MeV. C'est le facteur `ε ≈ 1,03` de la [leçon 6]({{ '/posts/ralentir-les-neutrons/' | relative_url }}) : 3 % des fissions
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

## À retenir

- La courbe `B/A` monte jusqu'à Ni-62 puis redescend : couper un noyau lourd
  en deux fragments de masse moyenne libère environ 0,8 MeV par nucléon.
- Le responsable est le terme **coulombien** en `Z²/A^(1/3)` : la répulsion
  entre protons croît plus vite que la liaison nucléaire.
- **Fissile** (U-235, Pu-239, U-233 : N impair) : la capture d'un neutron
  lent apporte `S_n` > barrière, grâce à l'énergie d'appariement.
  **Fissionnable** seulement (U-238 : N pair) : il faut un neutron de plus de ~1 MeV.
- Une fission libère ~207 MeV, dont ~200 MeV récupérables : l'essentiel en
  énergie cinétique des fragments, 12 MeV perdus en antineutrinos.
- Le partage **symétrique** libère le plus d'énergie, mais la fission de
  l'uranium est **asymétrique** : elle suit la surface d'énergie potentielle,
  façonnée par les couches Z = 50 et N = 82 des fragments. Comprendre ces
  effets de couches et de forme est l'objet de la partie II.

---

## Références

- **N. Bohr, J.A. Wheeler**, *The mechanism of nuclear fission*,
  **Phys. Rev. 56 (1939) 426** — l'article fondateur de la théorie de la
  fission par la goutte liquide.
- **P. Reuss, _Précis de neutronique_** (EDP Sciences), chapitres
  d'introduction — la fission vue du neutronicien.
- **Table of Nuclides** (KAERI, NNDC) — pour vérifier les énergies de liaison
  citées ici.

{% include cours-pied.html %}
