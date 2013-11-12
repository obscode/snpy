#!/bin/sh

rm *.eps
rm *.snpy

./test_templates.py
./test_plots.py
./test_fits.py
./test_kcorr.py

