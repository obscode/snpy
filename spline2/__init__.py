''' Spline2.py: wrapper for B. Thijsse et al.'s hyper-spline routines.

Yet another spline interpolation routine.  The problem:  given a set of
experimental data with noise, find the spline with the optimal number of
knots.

Solution :  They use the usual kind of routines to determine least-squares
            splines from a given set of knot points.  The problem REALLY
            boils down to:  how many knots do you use?  There are two 
            extremes:  put a knot point on each data point to get an
            interpolating spline (which sucks for experimental data with
            noise).  The other extreme is to have the minimal set of knots
            to define a polynomial of order k (e.g., a cubic).  This also
            sucks.  Somewhere between the two extremes is a number of
            knots that optimally recovers the information in the data and
            smooths out the noise.

            spline2 starts with a large number of knots (interpolating
            spline) and iteratively removes knots until a figure of merit
            reaches some prescribed value.  In this case, this figure of
            merit is the Durbin-Watson statistic, which measures the auto-
            correlation between the residuals of the spline fit.

For more details, see:
 *  Barend J. Thijsse et al., "A Practical Algorithm for Least-Squares 
    spline Approximation of Data Containing Noise", Computers in Physics, 
    vol 12 no. 4 July 1998
 *  http://structureandchange.3me.tudelft.nl/
'''
from spline2 import *
