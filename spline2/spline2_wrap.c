/* Implementation : PYTHON */

#include <string.h>
#include <stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "Python.h"
#ifdef __cplusplus
}
#endif

/* Definitions for Windows/Unix exporting */
#if defined(__WIN32__)
#   if defined(_MSC_VER)
#   define SWIGEXPORT(a,b) __declspec(dllexport) a b
#   else
#   if defined(__BORLANDC__)
#       define SWIGEXPORT(a,b) a _export b
#   else
#       define SWIGEXPORT(a,b) a b
#   endif
#   endif
#else
#   define SWIGEXPORT(a,b) a b
#endif

/* #include <numpy/oldnumeric.h> */
#include <numpy/arrayobject.h>
#include "spline2.h"

static PyObject *_wrap_spline2(PyObject *self, PyObject *args) {
    PyArrayObject *x_a, *y_a, *dy_a, *tt_a, *cc_a;
    PyObject *ret=NULL;
    int npoints, order, acfsearch, acf_ind, xflag, Xflag, logtrans;
    int nset, ksiset, n_min, n_max, rel, fixed_sigma, lind, lopt;
    int allownonopt, interactive, verbose;
    int n_knots, ret_lfin, kk;
    double *tt, *cc, ret_rms, ret_dws, ret_ksi, ret_acffit;
    double *dptr_x, *dptr_y, *dptr_dy;
    int  result, i;
    npy_intp dims[1];

    double xbegin, xend, Xbegin, Xend,ksibegin, ksiend, rejlev, fixval;

    self = self;

    if(!PyArg_ParseTuple(args,"OOOiiiiddiddiiiiidddiidiiiii",
             &x_a, &y_a, &dy_a, &order, &acfsearch, &acf_ind, &xflag,
             &xbegin, &xend, &Xflag, &Xbegin, &Xend, &logtrans,
             &nset, &ksiset, &n_min, &n_max, &ksibegin, &ksiend,
             &rejlev, &rel, &fixed_sigma, &fixval, &lind, &lopt, 
             &allownonopt, &interactive, &verbose))
        return NULL;
    /* Test the inputs */
    if (! PyArray_Check(x_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 1 must be an array");
       return NULL;
    }
    if (! PyArray_Check(y_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 2 must be an array");
       return NULL;
    }
    if (! PyArray_Check(dy_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 3 must be an array");
       return NULL;
    }

    if (PyArray_DIMS(x_a)[0] != PyArray_DIMS(y_a)[0] ||
        PyArray_DIMS(dy_a)[0] != PyArray_DIMS(x_a)[0]) {
       PyErr_SetString(PyExc_TypeError, "All arrays must have the same length");
       return NULL;
    }

    dptr_x = (double *) PyArray_DATA(x_a);
    dptr_y = (double *) PyArray_DATA(y_a);
    dptr_dy = (double *) PyArray_DATA(dy_a);
    npoints = PyArray_DIMS(x_a)[0];
    result = spline2(dptr_x, dptr_y, dptr_dy, npoints, order, acfsearch,
          acf_ind, xflag, xbegin, xend, Xflag, Xbegin, Xend, logtrans,
          nset, ksiset, n_min, n_max, ksibegin, ksiend, rejlev, rel,
          fixed_sigma, fixval, lind, lopt, allownonopt, interactive,
          verbose, &tt, &n_knots, &cc, &kk, &ret_rms, &ret_dws, &ret_lfin,
          &ret_ksi, &ret_acffit);
    if (result < 0) {
       PyErr_SetString(PyExc_RuntimeError, "spline2 failed");
       return NULL;
    }
    
    /* Create the output arrays */
    dims[0] = (npy_intp) n_knots;
    tt_a = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    for (i = 0 ; i < n_knots ; ++i) {
       dptr_x = (double *) PyArray_GETPTR1(tt_a, i);
       *dptr_x = tt[i];
    }
    dims[0] = (npy_intp) (n_knots + kk - 2);
    cc_a = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    for (i = 0 ; i < n_knots + kk -2 ; ++i) {
       dptr_x = (double *) PyArray_GETPTR1(cc_a, i);
       *dptr_x = cc[i];
    }
    /* these were malloc'ed by the code, so free them */
    free(tt);
    free(cc);
    
    ret = Py_BuildValue("NNiddidd", tt_a, cc_a, kk, ret_rms, ret_dws,
          ret_lfin, ret_ksi, ret_acffit);
    return(ret);
}


static PyObject *_wrap_evalsp(PyObject *self, PyObject *args) {
    PyArrayObject *x_a;
    PyArrayObject *in_t_a, *in_c_a, *out_y_a;
    int deriv;
    int numpoints, in_k, in_l;
    double *in_t, *in_c, *out_y, *in_x;
    int result, i;
    PyObject *ret;
    npy_intp dims[1];

    self = self;

    if(!PyArg_ParseTuple(args,"OiiOOi",
             &x_a, &in_k, &in_l, &in_t_a, &in_c_a, &deriv))
        return NULL;
    /* Test the inputs */
    if (! PyArray_Check(x_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 1 must be an array");
       return NULL;
    }
    if (! PyArray_Check(in_t_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 3 must be an array");
       return NULL;
    }
    if (! PyArray_Check(in_c_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 4 must be an array");
       return NULL;
    }

    numpoints = PyArray_DIMS(x_a)[0];

    in_x = (double *) PyArray_DATA(x_a);
    in_t = (double *) PyArray_DATA(in_t_a);
    in_c = (double *) PyArray_DATA(in_c_a);

    result = evalsp(in_x, numpoints, in_k, in_l, in_t, in_c, deriv, &out_y);

    /* Make an array */
    dims[0] = (npy_intp) numpoints;
    out_y_a = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    for (i = 0 ; i < numpoints ; ++i) {
       in_x = (double *) PyArray_GETPTR1(out_y_a, i);
       *in_x = out_y[i];
    }
    free(out_y);

    /* No need to DECREF, since "N" steals a reference */
    ret = Py_BuildValue("N", out_y_a);

    return(ret);
}


static PyObject *_wrap_eval_extrema(PyObject *self, PyObject *args) {
    PyArrayObject *in_t_a, *in_c_a, *out_xextr, *out_yextr, *out_signs;
    int in_k, in_l, n_extr, *signs;
    double *in_t, *in_c, *xextr, *yextr;
    double *dptr1, *dptr2;
    int *iptr;
    int result, i;
    PyObject *ret;
    npy_intp dims[1];

    self = self;

    if(!PyArg_ParseTuple(args,"iiOO",
             &in_k, &in_l, &in_t_a, &in_c_a))
        return NULL;
    /* Test the inputs */
    if (! PyArray_Check(in_t_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 3 must be an array");
       return NULL;
    }
    if (! PyArray_Check(in_c_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 4 must be an array");
       return NULL;
    }

    in_t = (double *) PyArray_DATA(in_t_a);
    in_c = (double *) PyArray_DATA(in_c_a);

    result = eval_extrema(in_k, in_l, in_t, in_c, &xextr, &yextr, &signs,&n_extr);

    /* Make the arrays */
    dims[0] = (npy_intp) n_extr;
    out_xextr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    out_yextr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    out_signs = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
    for (i = 0 ; i < n_extr ; ++i) {
       dptr1 = PyArray_GETPTR1(out_xextr, i);
       dptr2 = PyArray_GETPTR1(out_yextr, i);
       iptr = PyArray_GETPTR1(out_signs, i);
       *dptr1 = xextr[i];
       *dptr2 = yextr[i];
       *iptr = signs[i];
    }
    free(xextr);
    free(yextr);
    free(signs);

    /* No need to DECREF, since "N" steals a reference */
    ret = Py_BuildValue("NNN", out_xextr, out_yextr, out_signs);

    return(ret);
}


static PyObject *_wrap_eval_inflect(PyObject *self, PyObject *args) {
    PyArrayObject *in_t_a, *in_c_a, *out_x_a, *out_y_a, *out_dy_a;
    int in_k, in_l, n_inflect;
    double *in_t, *in_c, *out_y, *out_dy, *out_x;
    double *dptr1, *dptr2, *dptr3;
    int result, i;
    PyObject *ret;
    npy_intp dims[1];

    self = self;

    if(!PyArg_ParseTuple(args,"iiOO",
             &in_k, &in_l, &in_t_a, &in_c_a))
        return NULL;
    /* Test the inputs */
    if (! PyArray_Check(in_t_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 3 must be an array");
       return NULL;
    }
    if (! PyArray_Check(in_c_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 4 must be an array");
       return NULL;
    }

    in_t = (double *) PyArray_DATA(in_t_a);
    in_c = (double *) PyArray_DATA(in_c_a);

    result = eval_inflect(in_k, in_l, in_t, in_c, &out_x, &out_y, &out_dy,&n_inflect);

    /* Make the arrays */
    dims[0] = (npy_intp) n_inflect;
    out_x_a = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
    out_y_a = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
    out_dy_a = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
    for (i = 0 ; i < n_inflect ; ++i) {
       dptr3 = (double *) PyArray_GETPTR1(out_x_a, i);
       dptr2 = (double *) PyArray_GETPTR1(out_y_a, i);
       dptr1 = (double *) PyArray_GETPTR1(out_dy_a, i);
       *dptr3 = out_x[i];
       *dptr1 = out_y[i];
       *dptr2 = out_dy[i];
    }
    free(out_x);
    free(out_y);
    free(out_dy);

    /* No need to DECREF, since "N" steals a reference */
    ret = Py_BuildValue("NNN", out_x_a, out_y_a, out_dy_a);

    return(ret);
}

static PyObject *_wrap_eval_integ(PyObject *self, PyObject *args) {
    PyArrayObject *in_t_a, *in_c_a;
    double x1, x2;
    int in_k, in_l;
    double *in_t, *in_c;
    double result;
    PyObject *ret;

    self = self;

    if(!PyArg_ParseTuple(args,"ddiiOO",
             &x1, &x2, &in_k, &in_l, &in_t_a, &in_c_a))
        return NULL;
    /* Test the inputs */
    if (! PyArray_Check(in_t_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 5 must be an array");
       return NULL;
    }
    if (! PyArray_Check(in_c_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 6 must be an array");
       return NULL;
    }

    in_t = (double *) PyArray_DATA(in_t_a);
    in_c = (double *) PyArray_DATA(in_c_a);

    result = eval_integ(x1, x2, in_k, in_l, in_t, in_c);

    ret = Py_BuildValue("d", result);

    return(ret);
}

static PyObject *_wrap_eval_x(PyObject *self, PyObject *args) {
    PyArrayObject *in_t_a, *in_c_a, *out_roots;
    double value;
    int in_k, in_l, n_roots;
    double *in_t, *in_c, *roots;
    double *dptr;
    int result, i;
    PyObject *ret;
    npy_intp dims[1];

    self = self;

    if(!PyArg_ParseTuple(args,"diiOO",
             &value, &in_k, &in_l, &in_t_a, &in_c_a))
        return NULL;
    /* Test the inputs */
    if (! PyArray_Check(in_t_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 4 must be an array");
       return NULL;
    }
    if (! PyArray_Check(in_c_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 5 must be an array");
       return NULL;
    }

    in_t = (double *) PyArray_DATA(in_t_a);
    in_c = (double *) PyArray_DATA(in_c_a);

    result = eval_x(value, in_k, in_l, in_t, in_c, &roots, &n_roots);

    /* Make the arrays */
    dims[0] = (npy_intp) n_roots;
    out_roots = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    dptr = (double *) out_roots->data;
    for (i = 0 ; i < n_roots ; ++i) {
       dptr = (double *) PyArray_GETPTR1(out_roots, i);
       *dptr = roots[i];
    }
    free(roots);

    /* No need to DECREF, since "N" steals a reference */
    ret = Py_BuildValue("N", out_roots);

    return(ret);
}

static PyMethodDef spline2cMethods[] = {
    { "spline2", _wrap_spline2, 1 },
    { "evalsp", _wrap_evalsp, 1 },
    { "eval_extrema", _wrap_eval_extrema, 1 },
    { "eval_inflect", _wrap_eval_inflect, 1 },
    { "eval_integ", _wrap_eval_integ, 1 },
    { "eval_x", _wrap_eval_x, 1 },
    { NULL, NULL }
};
#ifdef __cplusplus
extern "C" 
#endif
SWIGEXPORT(void,initspline2c)(void) {
    PyObject *m, *d;
    m = Py_InitModule("spline2c", spline2cMethods);
    d = PyModule_GetDict(m);

   import_array();
}
