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
python backbone/iderare_prep.py --mode "$mode" --trimming "$trimming"

source pipeline.sh