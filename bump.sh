#!/usr/bin/env bash

set -euo pipefail
v=$1

commit_message="Bump version: ${v}"

sed -i "" "s/\(Version:\).*/\1 ${v}/" DESCRIPTION
git add DESCRIPTION
git commit -m "${commit_message}"
git tag -a "${v}" -m "bump v${v}"
git push origin "${v}"
