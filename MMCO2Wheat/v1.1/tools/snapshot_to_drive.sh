#!/bin/bash
# Copy the committed code of this branch (no .git, no build files) to the
# Google Drive snapshot folder MMCO2Wheat/v1.1/code_snapshot/, replacing the
# previous snapshot. The live git copy stays in ~/code (GitHub-backed).
set -e
G="$HOME/Library/CloudStorage/GoogleDrive-David.Helman@mail.huji.ac.il/My Drive/0_Tools/3_GitHub/M-M-CO2Crop-Model/MMCO2Wheat/v1.1/code_snapshot"
cd "$(git rev-parse --show-toplevel)"
[ -z "$(git status --porcelain MMCO2Wheat/v1.1)" ] || { echo "commit first"; exit 1; }
rm -rf "$G.tmp" && mkdir -p "$G.tmp"
git archive HEAD MMCO2Wheat/v1.1 | tar -x -C "$G.tmp" --strip-components=2
echo "branch $(git rev-parse --abbrev-ref HEAD), commit $(git rev-parse --short HEAD), $(date '+%Y-%m-%d %H:%M')" > "$G.tmp/SNAPSHOT.txt"
rm -rf "$G" && mv "$G.tmp" "$G"
echo "snapshot -> $G"
