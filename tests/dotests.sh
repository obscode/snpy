#!/bin/sh

rm *.eps
rm *.snpy

echo "Testing load"
~/snpy2/bin/python test_load.py
echo "Testing templates"
~/snpy2/bin/python test_templates.py
echo "Testing plots"
~/snpy2/bin/python test_plots.py
echo "Testing fitting"
~/snpy2/bin/python test_fits.py
echo "Testing k-corrections"
~/snpy2/bin/python test_kcorr.py

