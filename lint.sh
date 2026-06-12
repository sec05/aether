#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

if [ "$#" -gt 0 ]; then
  files=("$@")
else
  mapfile -t files < <(find include src -type f \
    \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \))
fi

if [ "${#files[@]}" -eq 0 ]; then
  exit 0
fi

python3 cpplint.py --linelength=100 "${files[@]}"