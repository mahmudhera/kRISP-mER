#!/usr/bin/env bash
# Fetch Azimuth 2.0 for Rule Set 2 on-target scoring.
#
# Azimuth 2.0 is distributed inside CRISPOR. kRISP-meR does not bundle it,
# so this script clones CRISPOR to the location kRISP-meR looks in first:
#
#     ~/crisporWebsite/bin/Azimuth-2.0
#
# Clone somewhere else by passing a directory as the first argument, then
# export AZIMUTH_DIR as printed at the end.
#
# Usage:  bash setup_azimuth.sh [target-directory]

set -euo pipefail

TARGET="${1:-$HOME/crisporWebsite}"
AZIMUTH="$TARGET/bin/Azimuth-2.0"

if [ -d "$AZIMUTH/azimuth" ]; then
    echo "Azimuth 2.0 is already present at:"
    echo "    $AZIMUTH"
elif [ -d "$TARGET" ]; then
    echo "ERROR: $TARGET exists but does not contain bin/Azimuth-2.0/azimuth." >&2
    echo "Remove it, or pass a different directory." >&2
    exit 1
else
    command -v git >/dev/null 2>&1 || { echo "ERROR: git is not installed." >&2; exit 1; }
    echo "Cloning CRISPOR into $TARGET ..."
    git clone --depth 1 https://github.com/maximilianh/crisporWebsite "$TARGET"
    if [ ! -d "$AZIMUTH/azimuth" ]; then
        echo "ERROR: cloned, but $AZIMUTH/azimuth is missing." >&2
        echo "CRISPOR's layout may have changed. Set AZIMUTH_DIR by hand." >&2
        exit 1
    fi
    echo "Done."
fi

echo
if [ "$TARGET" = "$HOME/crisporWebsite" ]; then
    echo "This is the default location, so no environment variable is needed."
else
    echo "Add this to your shell profile:"
    echo
    echo "    export AZIMUTH_DIR=\"$AZIMUTH\""
fi

echo
echo "Verify the scorer with:"
echo
echo "    python -m krispmer.rs2_azimuth"
echo
echo "It should print a score of 55 for the reference 30-mer."
