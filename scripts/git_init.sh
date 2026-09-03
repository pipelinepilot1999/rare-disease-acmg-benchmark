#!/usr/bin/env bash
set -euo pipefail
cd ~/rdacmg/repo
if [ ! -d .git ]; then
  git init -q -b main
  git config user.email "sathvikcisco@gmail.com"
  git config user.name "sathviksai"
fi
git add -A
git status --short | head -40
echo "---"
git status --short | wc -l
