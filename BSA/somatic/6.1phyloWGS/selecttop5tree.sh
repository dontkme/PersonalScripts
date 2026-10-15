#!/bin/bash
# Select the top 5 phyloWGS trees by log-likelihood.
#
# Usage: bash selecttop5tree.sh <SAMPLE>
#
# Reads  <SAMPLE>.summ.json.gz and <SAMPLE>.mutass.zip
# Writes <SAMPLE>.top5_tree_ids.tsv, <SAMPLE>.top5.summ.json.gz, <SAMPLE>.top5.mutass.zip
#
# Uses the base python3 explicitly so it still works when the py27 env
# (which has no python3) is active.

set -euo pipefail

SAMPLE="${1:?usage: selecttop5tree.sh <SAMPLE>}"
PY3="${PYTHON3:-/usr/bin/python3}"

"$PY3" - "${SAMPLE}" <<'PY'
import sys
import gzip
import json
import zipfile

sample = sys.argv[1]

summary_file = f"{sample}.summ.json.gz"
mutass_file = f"{sample}.mutass.zip"

with gzip.open(summary_file, "rt") as f:
    summary = json.load(f)

trees = summary["trees"]

top5 = sorted(
    trees.items(),
    key=lambda item: float(item[1]["llh"]),
    reverse=True
)[:5]

print("rank\ttree_id\tllh")
for rank, (tree_id, tree) in enumerate(top5, start=1):
    print(f"{rank}\t{tree_id}\t{tree['llh']}")

with open(f"{sample}.top5_tree_ids.tsv", "w") as f:
    f.write("rank\ttree_id\tllh\n")
    for rank, (tree_id, tree) in enumerate(top5, start=1):
        f.write(f"{rank}\t{tree_id}\t{tree['llh']}\n")

top5_summary = dict(summary)
top5_summary["trees"] = {
    tree_id: trees[tree_id]
    for tree_id, _ in top5
}

with gzip.open(f"{sample}.top5.summ.json.gz", "wt") as f:
    json.dump(top5_summary, f)

with zipfile.ZipFile(mutass_file) as source, \
     zipfile.ZipFile(
         f"{sample}.top5.mutass.zip",
         "w",
         compression=zipfile.ZIP_DEFLATED
     ) as target:

    available = set(source.namelist())

    for tree_id, _ in top5:
        member = f"{tree_id}.json"
        if member not in available:
            raise RuntimeError(
                f"{member} not found in {mutass_file}"
            )
        target.writestr(member, source.read(member))
PY
