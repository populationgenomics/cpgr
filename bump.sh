#!/usr/bin/env bash

set -euo pipefail
v=$1

commit_message="Bump version: ${v}"

sed -i "" "s/\(Version:\).*/\1 ${v}/" DESCRIPTION
sed -i "" "s/CPGR_VERSION/${v}/" deploy/conda-recipe/meta.yaml
git add DESCRIPTION deploy/conda-recipe/meta.yaml
git commit -m "${commit_message}"
