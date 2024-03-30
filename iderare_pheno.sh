#!/bin/bash
set -euo pipefail

conda activate iderare
python backbone/iderare_pheno.py

## Uncomment this line below if you would like to directly add the transformed HPO terms directly to iderare.yml
# python backbone/iderare_pheno_to_yml.py

conda deactivate