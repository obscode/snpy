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

/* #include <numpy/oldnumeric.h> */
#include <numpy/arrayobject.h>
#include "dm15temp.h"

static PyObject *_wrap_load_data(PyObject *self, PyObject *args) {
    int res;
    char *path;

    PyObject *result=NULL;

    self = self;

    if(!PyArg_ParseTuple(args,"s", &path)) return NULL;

    res = load_data(path);

    result = PyInt_FromLong((long) res);
    return(result);
}

static PyObject *_set_data(PyObject *self, PyObject *args) {
   char *filter;
   int i, j;
   int nchars;
   PyObject *result;
   PyObject *data, *list;
   PyObject *st;
   double *dp;
   int status=0;
   npy_intp *dims;
   if(!PyArg_ParseTuple(args,"sOO", &filter, &data, &list))
      return NULL;
   if (strcmp(filter, "B") == 0) {
      i = 0;
   } else if (strcmp(filter, "V") == 0) {
      i = 1;
   } else if (strcmp(filter, "u") == 0) {
      i = 2;
   } else if (strcmp(filter, "g") == 0) {
      i = 3;
   } else if (strcmp(filter, "r") == 0) {
      i = 4;
   } else if (strcmp(filter, "i") == 0) {
      i = 5;
   } else if (strcmp(filter, "Y") == 0) {
      i = 6;
   } else if (strcmp(filter, "J") == 0) {
      i = 7;
   } else if (strcmp(filter, "H") == 0) {
      i = 8;
   } else if (strcmp(filter, "K") == 0) {
      i = 9;
   } else {
      PyErr_SetString(PyExc_ValueError, "Unrecognized filter name");
      return NULL;
   }
   if (! PyArray_Check(data)){
      PyErr_SetString(PyExc_TypeError, "second argument should be a 2D array");
      return NULL;
   }
   if (PyArray_NDIM(data) != 2) {
      PyErr_SetString(PyExc_TypeError, "second argument should be a 2D array");
      return NULL;
   }
   dims = PyArray_DIMS(data);
   if (dims[0] != 4){
      PyErr_SetString(PyExc_ValueError, "data array must by 4XN");
      return NULL;
   }
   if (dims[1] > max_points) {
      PyErr_SetString(PyExc_ValueError, "data array must by 4XN and N < max_points");
      return NULL;
   }

   /* Check on the list */
   if ( ! PyList_Check(list)) {
      PyErr_SetString(PyExc_ValueError, "second argument must be a list type");
      return NULL;
   }
   if ( PyList_Size(list) != (Py_ssize_t) dims[1] ) {
      PyErr_SetString(PyExc_ValueError, "length of list be equal second dimension of data");
      return NULL;
   }

   np[i] = dims[1];
   for ( j = 0 ; j < dims[1] ; ++j) {
      dp = (double *) PyArray_GETPTR2(data, 0, j);
      x[i][j] = *dp;
      dp = (double *) PyArray_GETPTR2(data, 1, j);
      y[i][j] = *dp;
      dp = (double *) PyArray_GETPTR2(data, 2, j);
      z[i][j] = *dp;
      dp = (double *) PyArray_GETPTR2(data, 3, j);
      wz[i][j] = *dp;
      st = PyList_GetItem(list, (Py_ssize_t) j);
      if (! PyString_Check(st) ) {
         PyErr_SetString(PyExc_ValueError, "list elements must be strings");
         return NULL;
      }
      nchars = snprintf(name[i][j], 80, "%s", PyString_AsString(st));
   }
   result = Py_BuildValue("i", status);
   return result;
}

static PyObject *_get_data(PyObject *self, PyObject *args) {
   char *filter;
   int i, j;
   PyObject *result;
   PyObject *data, *list;
   npy_intp dims[2];
   double *dp;

   if(!PyArg_ParseTuple(args,"s", &filter))
      return NULL;
   if (strcmp(filter, "B") == 0) {
      i = 0;
   } else if (strcmp(filter, "V") == 0) {
      i = 1;
   } else if (strcmp(filter, "u") == 0) {
      i = 2;
   } else if (strcmp(filter, "g") == 0) {
      i = 3;
   } else if (strcmp(filter, "r") == 0) {
      i = 4;
   } else if (strcmp(filter, "i") == 0) {
      i = 5;
   } else if (strcmp(filter, "Y") == 0) {
      i = 6;
   } else if (strcmp(filter, "J") == 0) {
      i = 7;
   } else if (strcmp(filter, "H") == 0) {
      i = 8;
   } else if (strcmp(filter, "K") == 0) {
      i = 9;
   } else {
      PyErr_SetString(PyExc_TypeError, "Unrecognized filter name");
      return NULL;
   }
   dims[0] = (npy_intp) 4;
   dims[1] = (npy_intp) np[i];

   data = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
   list = PyList_New((Py_ssize_t) 0);
   for (j = 0 ; j < np[i] ; ++j) {
      dp = (double *) PyArray_GETPTR2(data, 0, j);
      *dp = x[i][j];
      dp = (double *) PyArray_GETPTR2(data, 1, j);
      *dp = y[i][j];
      dp = (double *) PyArray_GETPTR2(data, 2, j);
      *dp = z[i][j];
      dp = (double *) PyArray_GETPTR2(data, 3, j);
      *dp = wz[i][j];
      PyList_Append(list, PyString_FromString(name[i][j]));
   }
   result = Py_BuildValue("(OO)", data, list);
   Py_DECREF(data);
   Py_DECREF(list);
   return(result);
}


static PyObject *_wrap_dm15temp(PyObject *self, PyObject *args) {
    int size, filter, normalize;
    double dm15;
    double *t;
    double *mag;
    double *emag;
    int res;
    double sigx0 = 3.0;
    double scalex = 0.1;
    double maxsigx = 10.0;
    double sigy0 = 0.3;

    PyArrayObject *t_a, *mag_a, *emag_a;
    PyObject *result=NULL;

    self = self;

    if(!PyArg_ParseTuple(args,"idOOOi|dddd",
            &filter, &dm15, &t_a, &mag_a, &emag_a, &normalize,
            &sigx0, &scalex, &maxsigx, &sigy0))
        return NULL;

    /* Test the inputs */
    if (! PyArray_Check(t_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 3 must be an array");
       return NULL;
    }
    if (! PyArray_Check(mag_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 4 must be an array");
       return NULL;
    }
    if (! PyArray_Check(emag_a)) {
       PyErr_SetString(PyExc_TypeError, "parameter 5 must be an array");
       return NULL;
    }

    size = (int) PyArray_DIMS(t_a)[0];
    if (PyArray_DIMS(mag_a)[0] != size || PyArray_DIMS(emag_a)[0] != size) {
       PyErr_SetString(PyExc_TypeError, "t, mag, and emag must have equal size");
       return NULL;
    }

    t = (double *) PyArray_DATA(t_a);
    mag = (double *) PyArray_DATA(mag_a);
    emag = (double *) PyArray_DATA(emag_a);

    res = dm15temp(filter, dm15, t, size, mag, emag, normalize, sigx0, scalex,
          maxsigx, sigy0);

    result = PyInt_FromLong((long) res);
    return(result);
}

static PyMethodDef dm15temp2cMethods[] = {
   { "load_data", _wrap_load_data, 1},
   { "dm15temp", _wrap_dm15temp, 1},
   { "get_data", _get_data, 1},
   { "set_data", _set_data, 1},
   { NULL,NULL}
};

#ifdef __cplusplus
extern "C" 
#endif
SWIGEXPORT(void,initdm15temp2c)(void) {
    PyObject *m, *d;
    m = Py_InitModule("dm15temp2c", dm15temp2cMethods);
    d = PyModule_GetDict(m);

   import_array();
}
