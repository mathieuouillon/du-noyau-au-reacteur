---
title: Code
icon: fas fa-code
order: 2
---

Tous les chiffres du cours sortent de ces scripts Python, téléchargeables un
par un. Placez-les dans un même dossier : ils s'importent les uns les autres.

**Dépendances :** `numpy`, `scipy`, `matplotlib`, `sympy` et `periodictable`.

```bash
pip install numpy scipy matplotlib sympy periodictable endf remotezip
```

{% include cours-lecons.html %}

## Le code de chaque leçon

{% for l in cours_lecons -%}
- [{% unless l.annexe %}Leçon {{ l.lecon }} — {% endunless %}{{ l.title }}]({{ l.url | relative_url }}) :
  {% for f in l.code %}[`{{ f }}`]({{ '/assets/code/' | append: f | relative_url }}){% unless forloop.last %}, {% endunless %}{% endfor %}
{% endfor %}
## Relancer les calculs

- `python analyse_liaison.py` : `energies_liaison.csv`, `courbe_liaison.png`
  et l'analyse des masses (leçon 1).
- `python trace_liaison.py` : `energie_liaison_modeles.png` (leçon 1).
- `python etude_fission.py` : `fission_bilan.png` et les tableaux des leçons 2
  et 7.
- `python sections_efficaces.py` : noyau composé, taille quantique du neutron,
  résonances, loi en 1/v, passage de la barrière et partage des fragments
  (leçon 2), avec `sections_taille.png`, `fission_seuil.png` et
  `fission_partage.png`.
- `python sections_endf.py` : télécharge les fichiers U-235 et U-238 de
  ENDF/B-VIII.0 (~16 Mo, une seule fois, dans `donnees_endf/`), reconstruit
  les résonances avec `reconstruction.py`, valide la reconstruction et trace
  `endf_vue_ensemble.png`, `endf_resonances.png` et `endf_doppler.png`
  (leçon 2) ; environ 40 s. Demande en plus les paquets `endf` et
  `remotezip`.
- `python etude_modeles.py` : `modeles_masse.png` et les trois épreuves
  (leçon 3) ; environ 2 min.
- `python trace_termes.py` : `termes_goutte.png`, `termes_m2.png`,
  `termes_couches.png`, `termes_gain.png` et `termes_carte.png`, la
  contribution de chaque terme des modèles de masse (leçon 3).
- `python precalcul_micro.py` : `micro.npz`, la partie microscopique du modèle
  (leçon 4) ; environ 10 min, une seule fois.
- `python etude_frdm.py` : `modele_frdm.png` et l'étude complète (leçon 4) ;
  environ 2 min.
- `python etude_deformation.py` : les sections de détail de la leçon 4
  (niveaux sphériques, pentes de Nilsson, gap déformé N = 152, correction de
  couches, BCS face aux masses, samarium, moments quadrupolaires) et les dix
  figures `deformation_*.png` ; quelques secondes. Lit `micro.npz` et
  `be2_adopte.csv`, les B(E2) adoptés du NNDC (reconstruit depuis le site du
  NNDC s'il manque).
- `python galerie_orbitales.py` : `orbitales_nucleaires.png` et
  `orbitales_3d.json` (leçon 5).
- `python etude_barrieres.py` : `barrieres_fission.png`, les barrières de la
  goutte et la double bosse (leçon 6) ; environ 1 min.
- `python etude_coeur.py` : `cartes_puissance.png` et l'étude du cœur
  (leçon 8).
- `python run_tests.py` : la vérification de `diffusion.py` et `flux_maps.png`
  (annexe).

> [`micro.npz`]({{ '/assets/code/micro.npz' | relative_url }}) est le cache de la partie microscopique du modèle de la
> leçon 4 : le télécharger évite une dizaine de minutes de calcul.
{: .prompt-tip }

## Visualisation

- [Les orbitales nucléaires en 3D]({{ '/assets/visualisations/orbitales_3d.html' | relative_url }}),
  produite par `galerie_orbitales.py` (leçon 5).
