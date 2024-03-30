#!/bin/bash
set -euo pipefail

conda activate iderare

## This will use the default settings of threshold 0.4, differential 10, recommendation 20
python backbone/iderare_pheno.py

## Uncomment this line below if you would like to directly add the transformed HPO terms directly to iderare.yml and customize the parameters
# python backbone/iderare_pheno.py --threshold=0.5 --differential=10 --recommendation=25 --add-yml

conda deactivate