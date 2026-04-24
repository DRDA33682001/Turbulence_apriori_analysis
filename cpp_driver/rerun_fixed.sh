#!/usr/bin/env bash
set -euo pipefail

BUILD=/Users/dron/Documents/notes/Turbulence/Assignment/final_project/les_apriori/build
DATA=/Users/dron/Documents/notes/Turbulence/Assignment/final_project/les_apriori/data
BIN=$BUILD/les_apriori

cd "$BUILD"

for snap in "iso iso_t2 isotropic_256_t2.h5" \
            "iso iso_t3 isotropic_256_t3.h5" \
            "iso iso_t4 isotropic_256_t4.h5" \
            "channel channel_t1 channel_256_t1.h5" \
            "channel channel_t2 channel_256_t2.h5" \
            "channel channel_t3 channel_256_t3.h5" \
            "channel channel_t4 channel_256_t4.h5"; do
    read -r ds tag file <<< "$snap"
    echo "=== $tag ==="
    "$BIN" "$ds" all --input "$DATA/$file" --tag "$tag" 2>&1 | \
        grep -E "smagorinsky|wale|dynamic" | grep -v "^\["
    echo ""
done

echo "ALL DONE"
