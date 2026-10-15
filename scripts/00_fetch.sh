#!/usr/bin/env bash
# Fetch the public inputs. Re-running may produce a newer metadata.json than the one analysed;
# data/SHA256SUMS records the exact files used for the published results.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data
curl -sf https://vesuvius-challenge-open-data.s3.amazonaws.com/metadata.json | gunzip -c > data/metadata.json   # served gzip-encoded
curl -sf https://raw.githubusercontent.com/ScrollPrize/villa/main/scrollprize.org/src/data/prizeEligibility.json -o data/prizeEligibility.json
curl -sf https://scrollprize.org/prizes -o data/prizes_snapshot_$(date +%Y%m%d).html
# zarrsweep integrity report (https://github.com/sfc0324/zarrsweep, reports/full_sweep.json)
curl -sf https://raw.githubusercontent.com/sfc0324/zarrsweep/main/reports/full_sweep.json -o data/zarrsweep_full_sweep.json
