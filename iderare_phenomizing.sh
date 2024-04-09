#!/bin/bash
set -euo pipefail

# Default values for parameters
threshold=0.4
differential=10
recommendation=20
add_yml=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threshold)
            threshold="$2"
            shift 2
            ;;
        --differential)
            differential="$2"
            shift 2
            ;;
        --recommendation)
            recommendation="$2"
            shift 2
            ;;
        --add-yml)
            add_yml=true
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
if "$add_yml"; then
    python backbone/iderare_phenomizing.py --threshold "$threshold" --differential "$differential" --recommendation "$recommendation" --add-yml
else
    python backbone/iderare_phenomizing.py --threshold "$threshold" --differential "$differential" --recommendation "$recommendation"
fi

conda deactivate