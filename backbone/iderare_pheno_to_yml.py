#!/usr/bin/python

import os
import yaml

with open("iderare.yml") as i:
    y = yaml.safe_load(i)

    # Load data for trio analysis
    data_dir = y['analysis']['data_dir']
    proband = y['analysis']['proband']
    mother = y['analysis']['mother']
    father = y['analysis']['father']

    prob