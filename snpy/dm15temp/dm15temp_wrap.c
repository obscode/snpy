/* Implementation : PYTHON */

#include <string.h>
#include <stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "Python.h"
#include <float.h>
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

#ifdef SWIG_GLOBAL
#ifdef __cplusplus
#define SWIGSTATIC extern "C"
#else
#define SWIGSTATIC
#endif
#endif

#ifndef SWIGSTATIC
#define SWIGSTATIC static
#endif

#define SWIG_init    initdm15tempc
#define SWIG_name    "dm15temp"
/* #include <numpy/oldnumeric.h> */
#include <numpy/arrayobject.h>
#include "dm15temp.h"

static PyObject *_wrap_dm15temp(PyObject *self, PyObject *args) {
    PyObject * _resultobj;
    int  _result;
    int size;
    double dm15;
    double shiftv, shiftr, shifti;
    char *path;
    int npts;
    int method;
    PyObject *t, *tsig, *b, *bsig, *v, *vsig, *r, *rsig, *i, *isig;
    PyArrayObject *t_a, *tsig_a, *b_a, *bsig_a, *v_a, *vsig_a, *r_a, *rsig_a,
                  *i_a, *isig_a;

    if(!PyArg_ParseTuple(args,"diOOOOOOOOOOis", &dm15, &size, &t, &tsig, &b, 
             &bsig, &v, &vsig, &r, &rsig, &i, &isig, &method, &path))
        return NULL;

    if (! (t_a = (PyArrayObject *) PyArray_FROM_OTF(t, NPY_DOUBLE, 
                NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 3 must be a 1D array or sequence");
       return NULL;
    }
    if (! ( tsig_a = (PyArrayObject *) PyArray_FROM_OTF(tsig, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 4 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( b_a = (PyArrayObject *) PyArray_FROM_OTF(b, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 5 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( bsig_a = (PyArrayObject *) PyArray_FROM_OTF(bsig, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 6 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( v_a = (PyArrayObject *) PyArray_FROM_OTF(v, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 7 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( vsig_a = (PyArrayObject *) PyArray_FROM_OTF(vsig, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 8 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( r_a = (PyArrayObject *) PyArray_FROM_OTF(r, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 9 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( rsig_a = (PyArrayObject *) PyArray_FROM_OTF(rsig, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 10 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( i_a = (PyArrayObject *) PyArray_FROM_OTF(i, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 11 must be a 1D array or sequence");
       goto fail;
    }
    if (! ( isig_a = (PyArrayObject *) PyArray_FROM_OTF(isig, 
                NPY_DOUBLE, NPY_INOUT_ARRAY))) {
       PyErr_SetString(PyExc_TypeError, 
             "parameter 12 must be a 1D array or sequence");
       goto fail;
    }

    _result = dm15temp(dm15, size, (double *)PyArray_DATA(t_a), 
          (double *)PyArray_DATA(tsig_a), (double *) PyArray_DATA(b_a), 
          (double *) PyArray_DATA(bsig_a), (double *) PyArray_DATA(v_a),
          (double *) PyArray_DATA(vsig_a), (double *) PyArray_DATA(r_a), 
          (double *) PyArray_DATA(rsig_a), (double *) PyArray_DATA(i_a), 
          (double *) PyArray_DATA(isig_a), 
          &shiftv, &shiftr, &shifti, &npts, method, path);
    if (_result < 0) {
       PyErr_SetString(PyExc_TypeError, "Error generated in dm15 code.");
       goto fail;   /* check to see if there was an error */
    }
    _resultobj = Py_BuildValue("lddd", npts, shiftv, shiftr, shifti);
    Py_DECREF(t_a);
    Py_DECREF(tsig_a);
    Py_DECREF(b_a);
    Py_DECREF(bsig_a);
    Py_DECREF(v_a);
    Py_DECREF(vsig_a);
    Py_DECREF(r_a);
    Py_DECREF(rsig_a);
    Py_DECREF(i_a);
    Py_DECREF(isig_a);
   return _resultobj;
fail:
    Py_XDECREF(t_a);
    Py_XDECREF(tsig_a);
    Py_XDECREF(b_a);
    Py_XDECREF(bsig_a);
    Py_XDECREF(v_a);
    Py_XDECREF(vsig_a);
    Py_XDECREF(r_a);
    Py_XDECREF(rsig_a);
    Py_XDECREF(i_a);
    Py_XDECREF(isig_a);
    return NULL;

}

static PyMethodDef dm15tempMethods[] = {
    { "dm15temp", _wrap_dm15temp, 1 }, 
    { NULL, NULL }
};
#ifdef __cplusplus
extern "C" 
#endif
SWIGEXPORT(void,initdm15tempc)(void) {
    PyObject *m, *d;
    m = Py_InitModule("dm15tempc", dm15tempMethods);
    d = PyModule_GetDict(m);

   import_array();
}
