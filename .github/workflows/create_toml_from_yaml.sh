#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit

github_workspace=$1

cp ${github_workspace}/.github/environment/pixi.toml ${github_workspace}/pixi.toml
yq .requirements.host < ${github_workspace}/.github/recipe/recipe.yaml | \
    sed 's/- //' | sed 's/^/"/' | sed 's/$/" = "*"/'  >> ${github_workspace}/pixi.toml
