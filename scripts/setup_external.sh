#!/usr/bin/env bash
# Re-download upstream repos into /tmp/hubble_work_full/ for any heavy data
# that was deliberately excluded from external/.
set -euo pipefail
DEST=${1:-/tmp/hubble_work_full}
mkdir -p "$DEST"
cd "$DEST"
echo "Cloning into $DEST ..."
git clone --depth 1 https://github.com/Wang-weiYu/A-No-Go-guide-for-the-Hubble-tension.git
git clone --depth 1 https://github.com/dscolnic/Pantheon.git
git clone --depth 1 https://github.com/ytcosmo/TomoPost-BAO.git
git clone --depth 1 https://github.com/PantheonPlusSH0ES/DataRelease.git
git clone --depth 1 https://github.com/CobayaSampler/bao_data.git
git clone --depth 1 https://github.com/CobayaSampler/sn_data.git
echo "Done. Heavy items:"
du -sh */
