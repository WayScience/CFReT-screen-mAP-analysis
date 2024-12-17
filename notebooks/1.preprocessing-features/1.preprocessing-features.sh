#!/bin/bash

# activate conda env
conda activate cfret-map

# convert notebooks into python scripts
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# run the scripts
python nbconverted/1.preprocessing-features.py
