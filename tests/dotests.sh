#!/bin/sh

rm *.eps
rm *.snpy

~/snpy2/bin/python test_templates.py
~/snpy2/bin/python test_plots.py
~/snpy2/bin/python test_fits.py
~/snpy2/bin/python test_kcorr.py

