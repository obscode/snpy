#!/usr/bin/env python

from snpy import get_sn

# Test out the plotting capabilities. This doesn't test out the interactive
# stuff, just wheter or not we can plot things up.

s = get_sn('SN2006ax_kcorr.snpy')
s.plot(outfile="SN2006ax.eps")
s.plot_kcorrs(outfile="SN2006ax_kcorrs.eps")
s.plot_filters(outfile="SN2006ax_filters.eps")
s.plot_color("B", "V", outfile="SN2006ax_BV.eps")
