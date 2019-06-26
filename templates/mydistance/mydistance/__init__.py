'''
This is a custom distance metric that I use for pymc's GP package. They
include a isotropic euclidean distance metric, but I wanted one where you
could scale the dimensions separately.
'''

from aniso_euclidean_s import aniso_euclidean_s
aniso_euclidean_s.extra_parameters = {'sr': 'Scale ratio for each dimension'}
aniso_euclidean_s.__name__ = "aniso_euclidean_s"
aniso_euclidean_s.__doc__ = """
   aniso_euclidean_s(D, x, y, cmin=0, cmax=y.shape[0], scalerat, symmm=False)

   :Arguments:

      - `x and y` are arrays of points in Euclidean coordinates
        formatted as follows:
        
        [[x_{0,0} ... x_{0,ndim}],
         [x_{1,0} ... x_{1,ndim}],
         ...
         [x_{N,0} ... x_{N,ndim}]]
       - `symm` indicates whether x and y are references to
          the same array.
       - `cmin' and `cmax' indicate which columns to compute.
          These are used for multithreaded evaluation.
       - `scalerat` should have shape ndim and represent a scale
          ratio for the 1 .. ndim-1 dimensions

    Return value is a matrix D, where D[i,j] gives the Euclidean
    distance between the point x[i,:] and y[j,:]."""

