#!/bin/bash
set -euo pipefail

# Default values for parameters
mode="both"
trimming=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            mode="$2"
            shift 2
            ;;
        --trimming)
            trimming=true
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

conda activate iderare

# Pass parameters to Python script
## If trimming=True, then the script will run the trimming process --trimming, else, only pass mode
if [ "$trimming" = true ]; then
    python backbone/iderare_prep.py --mode "$mode" --trimming
else
    python backbone/iderare_prep.py --mode "$mode"
fi

source pipeline.sh