    ---
    title: "Du noyau au réacteur : présentation de la série"
    date: 2026-09-26 09:00:00 +0200
    categories: [Présentation]
    tags: [série, code]
    description: "Une série pour comprendre la neutronique et la physique nucléaire en construisant le code soi-même."
    image:
      path: /assets/img/nucleaire/energie_liaison_modeles.png
      alt: "Du noyau au réacteur : présentation de la série"
    pin: true
---


Cette série construit, pas à pas et en Python, de quoi comprendre un cœur de
réacteur : de la structure du noyau jusqu'au calcul neutronique d'un REP.
Chaque article s'appuie sur du code qui tourne, vérifié contre des solutions
exactes ou des données mesurées, et chaque erreur trouvée en chemin est
signalée plutôt qu'effacée.

## Les articles, dans l'ordre de lecture

1. [De l'uranium au réacteur]({{ '/posts/de-l-uranium-au-reacteur/' | relative_url }}) — Pourquoi l'uranium 235 fissionne, combien d'énergie sort, et pourquoi un REP ne peut pas fonctionner à l'uranium naturel.
2. [Énergies de liaison : les données mesurées]({{ '/posts/energies-de-liaison/' | relative_url }}) — 2825 nucléides de l'évaluation AME2020, leur provenance, et ce qu'ils révèlent : le maximum en Ni-62 et les couches nucléaires.
3. [Au-delà de la goutte liquide : modèles de masse]({{ '/posts/au-dela-de-la-goutte-liquide/' | relative_url }}) — Cinq modèles de masse jugés sur ce qu'ils prédisent : la physique extrapole, la statistique non.
4. [Vers FRDM : le modèle macroscopique-microscopique]({{ '/posts/vers-frdm/' | relative_url }}) — Goutte déformable, niveaux de Nilsson, correction de Strutinsky et appariement BCS : pourquoi les noyaux se déforment.
5. [Les orbitales nucléaires, en 2D et en 3D]({{ '/posts/orbitales-nucleaires/' | relative_url }}) — Tracer les orbitales des nucléons comme celles de l'électron, puis la forme du noyau entier.
6. [Le solveur de diffusion multigroupe (en anglais)]({{ '/posts/solveur-de-diffusion/' | relative_url }}) — Physique, discrétisation et itération d'un solveur de diffusion 2D vérifié contre une solution analytique.
7. [Neutronique du cœur : comprendre et concevoir]({{ '/posts/neutronique-du-coeur/' | relative_url }}) — Plan de chargement, bore critique, efficacité des grappes et anti-ombrage, sur un cœur REP modélisé.

Les cinq premiers remontent à la source : pourquoi l'uranium fissionne, ce que
disent les masses mesurées, et comment des modèles de plus en plus physiques les
expliquent. Les deux derniers redescendent vers le réacteur : le solveur de
diffusion, puis la conception d'un cœur.

## Le code

Tout le code est téléchargeable. Dépendances : `numpy`, `scipy`, `matplotlib`,
`sympy` et `periodictable`.

- [`analyse_liaison.py`]({{ '/assets/code/analyse_liaison.py' | relative_url }})
- [`coeur.py`]({{ '/assets/code/coeur.py' | relative_url }})
- [`diffusion.py`]({{ '/assets/code/diffusion.py' | relative_url }})
- [`donnees_liaison.py`]({{ '/assets/code/donnees_liaison.py' | relative_url }})
- [`energies_liaison.csv`]({{ '/assets/code/energies_liaison.csv' | relative_url }})
- [`etude_coeur.py`]({{ '/assets/code/etude_coeur.py' | relative_url }})
- [`etude_fission.py`]({{ '/assets/code/etude_fission.py' | relative_url }})
- [`etude_frdm.py`]({{ '/assets/code/etude_frdm.py' | relative_url }})
- [`etude_modeles.py`]({{ '/assets/code/etude_modeles.py' | relative_url }})
- [`fission.py`]({{ '/assets/code/fission.py' | relative_url }})
- [`frdm.py`]({{ '/assets/code/frdm.py' | relative_url }})
- [`galerie_orbitales.py`]({{ '/assets/code/galerie_orbitales.py' | relative_url }})
- [`micro.npz`]({{ '/assets/code/micro.npz' | relative_url }})
- [`modeles_masse.py`]({{ '/assets/code/modeles_masse.py' | relative_url }})
- [`nilsson.py`]({{ '/assets/code/nilsson.py' | relative_url }})
- [`orbitales.py`]({{ '/assets/code/orbitales.py' | relative_url }})
- [`precalcul_micro.py`]({{ '/assets/code/precalcul_micro.py' | relative_url }})
- [`run_tests.py`]({{ '/assets/code/run_tests.py' | relative_url }})
- [`strutinsky.py`]({{ '/assets/code/strutinsky.py' | relative_url }})
- [`trace_liaison.py`]({{ '/assets/code/trace_liaison.py' | relative_url }})

> `micro.npz` est le cache de la partie microscopique du modèle de type FRDM :
> le garder évite une dizaine de minutes de calcul.
{: .prompt-tip }
