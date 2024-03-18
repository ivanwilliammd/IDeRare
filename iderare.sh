#!/bin/bash
set -euo pipefail

conda activate iderare
python backbone/iderare_prep.py

source pipeline.sh