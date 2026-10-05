# Du noyau au réacteur

Site Jekyll (thème [Chirpy](https://github.com/cotes2020/jekyll-theme-chirpy) 7.6),
construit à partir du modèle officiel `chirpy-starter`.

## Publier avec `gh`

Prérequis : [git](https://git-scm.com) et la [CLI GitHub](https://cli.github.com).

```bash
./publier.sh                 # site utilisateur : https://<vous>.github.io/
./publier.sh noyau-reacteur  # site de projet : https://<vous>.github.io/noyau-reacteur/
```

Le script vous connecte à GitHub si besoin (`gh auth login`), remplace les
champs `VOTRE_UTILISATEUR` et `VOTRE_NOM`, crée le dépôt public, active GitHub
Pages en mode « GitHub Actions », pousse le site et suit le déploiement.

## Tester en local

Prérequis : Ruby 3 et Bundler.

```bash
bundle install
bundle exec jekyll serve      # http://127.0.0.1:4000
```

## Ajouter une leçon

Créer un article dans `_posts/` avec, dans son en-tête, `lecon` (numéro),
`partie` (`I`, `II` ou `III`), `objectifs`, `prerequis` (numéros de leçons) et
`code` (fichiers de `assets/code/`). Placer `{% include cours-entete.html %}`
au début du texte et `{% include cours-pied.html %}` à la fin. La date de
l'article doit suivre l'ordre des leçons. Le plan du cours, la navigation et
la page Code se mettent à jour d'eux-mêmes.

## Organisation

| chemin | contenu |
|---|---|
| `_posts/` | les leçons ; leur en-tête porte `lecon`, `partie`, `objectifs`, `prerequis` et `code` |
| `_data/cours.yml` | les parties du cours |
| `index.md`, `_layouts/home.html` | la page d'accueil « Plan du cours » |
| `_includes/cours-*.html`, `_includes/plan-cours.html` | en-tête et navigation des leçons, plan du cours |
| `_tabs/glossaire.md`, `_tabs/code.md` | le glossaire et la page du code |
| `_tabs/about.md` | la page « À propos » |
| `assets/img/nucleaire/` | les figures |
| `assets/code/` | le code Python téléchargeable |
| `assets/visualisations/orbitales_3d.html` | la visualisation 3D des orbitales |
| `.github/workflows/pages-deploy.yml` | construction et déploiement (modèle Chirpy) |
