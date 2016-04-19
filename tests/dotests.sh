#!/bin/sh

#export PY=/Users/cburns/snpy2/bin/python
export PY=python

rm *.eps
rm *.snpy
rm *.dat

echo "Testing load"
$PY test_load.py
echo "Testing templates"
$PY test_templates.py
echo "Testing plots"
$PY test_plots.py
echo "Testing fitting"
$PY test_fits.py
echo "Testing k-corrections"
$PY test_kcorr.py
echo "Testing bolometric calculator"
$PY test_bolometric.py
