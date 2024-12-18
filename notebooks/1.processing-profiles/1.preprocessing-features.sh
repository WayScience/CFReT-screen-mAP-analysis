#!/bin/bash

# This script updates the plate maps with MOA (Mechanism of Action) information and generates
# aggregate profiles at both the MOA level and the replicate level.

# activate conda env
conda activate cfret-map

# convert notebooks into python scripts
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# run the scripts
python nbconverted/1.preprocessing-features.py
