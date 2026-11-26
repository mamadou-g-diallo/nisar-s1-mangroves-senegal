#!/usr/bin/env bash
# Step 1: extract NISAR GCOV/GUNW and Sentinel-1 RTC on a common 20 m grid with nisar-s1 (needs ~/.netrc
# for NASA Earthdata), then per-class statistics with mangroves (WorldCover 95) as the target class.
# Outputs go to $SALOUM_WORK (default ./work): saloum_full/ and casamance/.
set -euo pipefail
cd "$(dirname "$0")/.."
WORK="${SALOUM_WORK:-work}"; mkdir -p "$WORK"
SITES=config/sites_senegal.yaml
echo "[1/2] Saloum, June-November 2026"
nisar-s1 run --sites-file "$SITES" --site saloum --start 2026-06-01 --end 2026-11-20 --out "$WORK/saloum_full" --target-class 95
cp "$WORK/saloum_full/separability.csv" "$WORK/saloum_full/separability_t95.csv"
echo "[2/2] Casamance (replication)"
nisar-s1 run --sites-file "$SITES" --site casamance --start 2026-06-01 --end 2026-11-23 --out "$WORK/casamance" --target-class 95
