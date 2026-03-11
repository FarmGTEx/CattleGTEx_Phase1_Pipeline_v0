#!/bin/bash

set -euo pipefail
BFILE="$1"
K="$2"

/faststorage/project/cattle_gtexs/software/admixture_linux-1.3.0/admixture -j12 --cv "$BFILE" "$K" |tee log${K}.out
