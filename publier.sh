#!/usr/bin/env bash
# publier.sh -- cree le depot GitHub et publie le site sur GitHub Pages avec gh.
#
#   ./publier.sh               -> depot du-noyau-au-reacteur, https://<vous>.github.io/du-noyau-au-reacteur/
#   ./publier.sh <nom-depot>   -> autre nom de depot, https://<vous>.github.io/<nom-depot>/
set -euo pipefail
cd "$(dirname "$0")"

command -v git >/dev/null || { echo "git est requis : https://git-scm.com"; exit 1; }
command -v gh  >/dev/null || { echo "gh est requis : https://cli.github.com"; exit 1; }

if ! gh auth status >/dev/null 2>&1; then
  echo "Connexion a GitHub (gh auth login)..."
  gh auth login
fi

export UTIL NOM
UTIL=$(gh api user --jq .login)
NOM=$(gh api user --jq '.name // .login')
DEPOT="${1:-du-noyau-au-reacteur}"
if [ "$DEPOT" = "$UTIL.github.io" ]; then export BASE=""; else export BASE="/$DEPOT"; fi
echo "Compte : $UTIL   Depot : $DEPOT   Adresse : https://$UTIL.github.io$BASE/"

# 1. Personnaliser la configuration
perl -pi -e 's/VOTRE_UTILISATEUR/$ENV{UTIL}/g; s/VOTRE_NOM/$ENV{NOM}/g' _config.yml _tabs/about.md
perl -pi -e 's|^baseurl: .*|baseurl: "$ENV{BASE}"|' _config.yml

# 2. Depot local
[ -d .git ] || git init -q -b main
git add -A
git commit -qm "Site : du noyau au reacteur" || echo "(rien de nouveau a valider)"

# 3. Depot GitHub (sans pousser tout de suite)
if gh repo view "$UTIL/$DEPOT" >/dev/null 2>&1; then
  echo "Le depot $UTIL/$DEPOT existe deja : on le reutilise."
  git remote get-url origin >/dev/null 2>&1 || git remote add origin "https://github.com/$UTIL/$DEPOT.git"
else
  gh repo create "$DEPOT" --public --source=. --remote=origin \
    --description "Physique nucleaire et neutronique, calculees pas a pas"
fi

# 4. Activer GitHub Pages en mode "GitHub Actions" AVANT le premier push
if ! gh api -X POST "repos/$UTIL/$DEPOT/pages" -f build_type=workflow >/dev/null 2>&1; then
  gh api -X PUT "repos/$UTIL/$DEPOT/pages" -f build_type=workflow >/dev/null 2>&1 || {
    echo "Activation automatique de Pages impossible."
    echo "Faites-la a la main : Settings > Pages > Source : GitHub Actions, puis relancez."
    exit 1; }
fi

# 5. Pousser : le workflow de construction et de deploiement demarre
git push -u origin main

# 6. Suivre le deploiement
echo "Attente du demarrage du workflow..."
sleep 8
RUN=$(gh run list -R "$UTIL/$DEPOT" --workflow pages-deploy.yml -L 1 --json databaseId --jq '.[0].databaseId')
gh run watch -R "$UTIL/$DEPOT" "$RUN" --exit-status
echo
echo "Site publie : https://$UTIL.github.io$BASE/"
echo "(la premiere mise en ligne peut demander une ou deux minutes de plus)"
