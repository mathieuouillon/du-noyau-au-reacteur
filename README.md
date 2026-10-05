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
`partie` (`I`, `II` ou `III`), `objectifs`, `prerequis` (numéros de leçons),
`code` (fichiers de `assets/code/`), `sources` (d'où viennent les chiffres) et
`bibliographie` (liste de `cle` et `note`). Placer `{% include cours-entete.html %}`
au début du texte, puis, à la fin, une section `## Bibliographie` suivie de
`{% include bibliographie.html %}` et de `{% include cours-pied.html %}`. La date de
l'article doit suivre l'ordre des leçons. Le plan du cours, la navigation et
la page Code se mettent à jour d'eux-mêmes.

## Écrire une formule

Les formules sont rendues par MathJax (activé pour toutes les leçons et tous
les onglets dans `_config.yml`). On écrit du LaTeX entre `$$ … $$` :

- **dans le texte** : `l'énergie $$B/A$$ culmine…` ;
- **en bloc** : `$$` seul sur une ligne, la formule, puis `$$`, avec une ligne
  vide avant et après.

Quelques règles : virgule décimale `0{,}726` ; valeur absolue `\lvert x \rvert`
(jamais `|`, qui casse les tableaux) ; ne jamais écrire deux accolades
ouvrantes collées, réservées à Liquid ; `\begin{aligned}` ou
`\begin{gathered}` pour couper une formule trop large. Les noms de fichiers,
de variables du code et les commandes restent entre accents graves.

## Ajouter une référence

Ajouter une entrée à `_data/bibliographie.yml` (`type` : `livre`, `article`,
`revue`, `donnees` ou `code` ; `auteurs`, `titre`, `revue`, `volume`, `pages`,
`annee`, `doi` ou `url`), puis citer sa clé dans l'en-tête des leçons
concernées. Une clé inconnue s'affiche en rouge sur la page.

## Organisation

| chemin | contenu |
|---|---|
| `_posts/` | les leçons ; leur en-tête porte `lecon`, `partie`, `objectifs`, `prerequis` et `code` |
| `_data/cours.yml` | les parties du cours |
| `index.md`, `_layouts/home.html` | la page d'accueil « Plan du cours » |
| `_includes/cours-*.html`, `_includes/plan-cours.html` | en-tête et navigation des leçons, plan du cours |
| `_tabs/glossaire.md`, `_tabs/code.md` | le glossaire et la page du code |
| `_data/bibliographie.yml` | toutes les références (DOI vérifiés), citées par leur clé |
| `_includes/bibliographie.html`, `_tabs/bibliographie.md` | bibliographie de chaque leçon et bibliographie complète |
| `_tabs/about.md` | la page « À propos » |
| `assets/img/nucleaire/` | les figures |
| `assets/code/` | le code Python téléchargeable |
| `assets/visualisations/orbitales_3d.html` | la visualisation 3D des orbitales |
| `.github/workflows/pages-deploy.yml` | construction et déploiement (modèle Chirpy) |
