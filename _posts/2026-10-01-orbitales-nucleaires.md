---
title: "Les orbitales nucléaires, en 2D et en 3D"
date: 2026-10-01 09:00:00 +0200
categories: [Physique nucléaire, Structure du noyau]
tags: [Nilsson, orbitales, visualisation, déformation]
description: "Tracer les orbitales des nucléons comme celles de l'électron, puis la forme du noyau entier."
image:
  path: /assets/img/nucleaire/orbitales_nucleaires.png
  alt: "Les orbitales nucléaires, en 2D et en 3D"
---

Suite de [Vers FRDM : le modèle macroscopique-microscopique]({{ '/posts/vers-frdm/' | relative_url }}). Le modèle de Nilsson ne donne pas seulement des
**énergies** : chaque niveau a une **fonction d'onde**, qu'on peut tracer
comme une orbitale électronique. En sommant les orbitales occupées, on
obtient la **forme du noyau**.

- `orbitales_nucleaires.png` — la galerie ci-dessous
- `orbitales_3d.html` — la version 3D interactive (publiée en page)
- [`orbitales.py`]({{ '/assets/code/orbitales.py' | relative_url }}) — le calcul ; [`galerie_orbitales.py`]({{ '/assets/code/galerie_orbitales.py' | relative_url }}) — la figure et les données 3D


<iframe src="{{ '/assets/visualisations/orbitales_3d.html' | relative_url }}" title="Orbitales nucléaires en 3D" loading="lazy" style="width:100%;height:640px;border:1px solid var(--main-border-color, #ddd);border-radius:6px;"></iframe>
_Version 3D interactive : faites tourner le nuage. [Ouvrir en plein écran]({{ '/assets/visualisations/orbitales_3d.html' | relative_url }})._

![Orbitales nucléaires]({{ '/assets/img/nucleaire/orbitales_nucleaires.png' | relative_url }})

## De l'atome au noyau

Pour l'électron de l'hydrogène, on trace |ψ|² des orbitales 1s, 2p, 3d. Pour un
nucléon, c'est le même geste, avec trois différences :

1. **Le puits n'est pas coulombien.** Le potentiel moyen ressemble à un puits
   à fond plat. Les niveaux suivent donc ceux d'un oscillateur : 1s, 1p, 1d,
   2s, 1f… Le nombre devant la lettre compte les nœuds radiaux plus un, et
   non la couche comme dans l'atome.
2. **Le spin-orbite est énorme**, de l'ordre du MeV, contre des
   millionièmes d'eV dans l'atome. On couple donc l et s en j dès le départ :
   1d5/2, 1g9/2, 1i13/2. C'est ce couplage qui crée les nombres magiques
   (M. Goeppert Mayer ; O. Haxel, J.H.D. Jensen, H.E. Suess, *Phys. Rev.* 75,
   1949).
3. **Le puits peut être déformé.** l n'est alors plus conservé : seule reste
   Ω, la projection de j sur l'axe du noyau. Chaque orbitale est un mélange
   de plusieurs (l, j) ; on la nomme par son « parent » sphérique dominant
   (S.G. Nilsson, *Mat. Fys. Medd. Dan. Vid. Selsk.* 29, n° 16, 1955).

## Le calcul

Dans la base de l'oscillateur |N l j Ω⟩, avec m = Ω − σ :

```
ψ_σ(r, θ, φ) = Σ_i c_i R_nl(r/b) ⟨l m ½ σ | j Ω⟩ Y_lm(θ, φ) / b^(3/2)
|ψ|² = |ψ_↑|² + |ψ_↓|²
```

Les coefficients c_i sont les vecteurs propres du hamiltonien de Nilsson.
`b = √(ħ/mω₀)` est la longueur de l'oscillateur, environ 1 A^(1/6) fm. La
densité ne dépend pas de l'angle φ (symétrie axiale) : une coupe dans un plan
contenant l'axe suffit, et la 3D s'obtient par révolution.

Contrôles : chaque orbitale est normée à 1 ; la somme des orbitales occupées
du plomb 208 contient exactement 82 protons et 126 neutrons.

## Ce que montrent les images

**Ligne 1, noyau sphérique.** Les analogues directs des orbitales atomiques.
Ω petit concentre le nuage vers les pôles ; Ω maximal le range dans un tore
équatorial. Tant que le noyau est sphérique, toutes ces orientations ont la
même énergie.

**Ligne 2, noyau allongé (δ = 0,25), famille 1i13/2.** Les mêmes quatorze
états ne sont plus dégénérés :

| Ω | 1/2 | 5/2 | 9/2 | 13/2 |
|---|---|---|---|---|
| énergie relative (MeV) | 0 | +1,4 | +4,1 | +7,6 |

Pour Ω petit, le moment angulaire est presque perpendiculaire à l'axe :
l'orbite passe par les pôles, là où le noyau allongé contient le plus de
matière, et l'énergie est basse. Pour Ω grand, le moment angulaire est aligné
sur l'axe et l'orbite tourne à l'équateur, là où le noyau est mince : énergie
haute. **C'est l'éventail du diagramme de Nilsson, rendu visible.** C'est
aussi le mécanisme de la déformation : en milieu de couche, les nucléons de
valence peuvent descendre sur les orbitales de Ω petit en allongeant le noyau.

**Ligne 3, noyaux entiers.**

| noyau | δ | rapport des rayons z / x | rayon rms des protons |
|---|---|---|---|
| Pb-208 | 0 | 1,00 | 5,56 fm |
| Er-166 | 0,24 | 1,32 | 5,27 fm |
| U-238 | 0,19 | 1,24 | 5,87 fm |

Les déformations δ sont celles trouvées par le modèle de [Vers FRDM : le modèle macroscopique-microscopique]({{ '/posts/vers-frdm/' | relative_url }}).

## Limites, à garder en tête

- **Rayon.** Pour le plomb, 5,56 fm, contre un rayon de charge mesuré de
  5,50 fm : bon accord, à 1 % près.
- **Densité centrale trop forte.** Environ 0,21 nucléon/fm³, contre ~0,16
  dans la matière nucléaire réelle (densité de saturation). Un oscillateur
  n'a pas de fond plat et concentre trop la matière au centre. Les
  ondulations visibles dans l'erbium et l'uranium sont en partie des effets
  de couches réels, en partie amplifiées par ce défaut.
- **Pas de peau de neutrons.** Le modèle donne le même rayon aux protons et
  aux neutrons, et c'est voulu : la correction d'isospin
  `ħω₀(1 ± (N−Z)/3A)` a été conçue pour égaliser les rayons. Or la peau de
  neutrons du plomb est mesurée positive (expérience PREX). Il faudrait un
  calcul auto-cohérent, de type Hartree-Fock, pour la prédire.
- **Déformation de l'uranium sous-estimée** par le modèle, voir
  [Vers FRDM : le modèle macroscopique-microscopique]({{ '/posts/vers-frdm/' | relative_url }}).

## Pour aller plus loin

- A. Bohr, B.R. Mottelson, *Nuclear Structure*, vol. II (Benjamin, 1975) :
  le chapitre sur les noyaux déformés, et les célèbres diagrammes de Nilsson.
- P. Ring, P. Schuck, *The Nuclear Many-Body Problem* (Springer, 1980),
  chap. 2.
- Rayons de charge : I. Angeli, K.P. Marinova, *At. Data Nucl. Data Tables*
  99 (2013) 69.
