#!/usr/bin/env bash
# Stitch the phase-3 frame PNGs (out/phase3/frames/<run>/frame_%05d.png)
# into one mp4 per run directory.
#
#   scripts/make_videos.sh [fps]      (default 24)
set -euo pipefail
cd "$(dirname "$0")/.."

fps="${1:-24}"
mkdir -p out/videos

shopt -s nullglob
found=0
for dir in out/phase3/frames/*/ out/phase3/N*/frames/*/; do
    run=$(basename "$dir")
    # prefix resolution-scoped runs (out/phase3/N144/frames/k2_ball -> N144_k2_ball)
    parent=$(basename "$(dirname "$(dirname "$dir")")")
    [[ $parent == N* ]] && run="${parent}_${run}"
    frames=("$dir"frame_*.png)
    ((${#frames[@]})) || continue
    found=1
    out="out/videos/phase3_${run}.mp4"
    echo ">> $run: ${#frames[@]} frames -> $out"
    ffmpeg -y -loglevel error -framerate "$fps" \
        -i "${dir}frame_%05d.png" \
        -c:v libx264 -crf 20 -pix_fmt yuv420p "$out"
done

if ((found)); then
    echo "Done:"
    ls -lh out/videos/
else
    echo "No frame directories with PNGs under out/phase3/frames/" >&2
    exit 1
fi
