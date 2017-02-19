%module rjmcmc

%{

#include "rjmcmc_helper.h"

%}

void 
rjmcmc_seed(int s);
  
typedef struct {
  dataset1d_t *d;
} dataset1d;

%extend dataset1d {

  dataset1d(const char *filename)
  {
    dataset1d *d = malloc(sizeof(dataset1d));
    d->d = dataset1d_load_known(filename);

    return d;
  }

  dataset1d(const char *filename, double n)
  {
    dataset1d *d = malloc(sizeof(dataset1d));
    d->d = dataset1d_load_fixed(filename, n);

    return d;
  }

  dataset1d(PyObject *x,
	    PyObject *y,
	    PyObject *n)
  {
    dataset1d *d = malloc(sizeof(dataset1d));
    int i;
    PyObject *po;

    d->d = NULL;

    if (!PyList_Check(x) ||
	!PyList_Check(y) ||
	!PyList_Check(n)) {
      PyErr_SetString(PyExc_ValueError,
		      "Parameter must be a list of floating point values");
      goto create_from_array_cleanup;
    }

    if (PyList_Size(x) != PyList_Size(y) ||
	PyList_Size(x) != PyList_Size(n)) {
      PyErr_SetString(PyExc_ValueError,
		      "Parameters must have the same length.");
      goto create_from_array_cleanup;
    }

    d->d = dataset1d_create((int)PyList_Size(x));

    d->d->xmin = FLT_MAX;
    d->d->xmax = -FLT_MAX;
    d->d->ymin = FLT_MAX;
    d->d->ymax = -FLT_MAX;

    for (i = 0; i < PyList_Size(x); i ++) {
      po = PyList_GetItem(x, i);
      if (PyNumber_Check(po)) {
	d->d->points[i].x = PyFloat_AsDouble(po);
      } else {
	PyErr_SetString(PyExc_ValueError,
			"Invalid element in list, must be a number");
	goto create_from_array_cleanup;
      }

      po = PyList_GetItem(y, i);
      if (PyNumber_Check(po)) { 
	d->d->points[i].y = PyFloat_AsDouble(po);
      } else {
	PyErr_SetString(PyExc_ValueError,
			"Invalid element in y list, must be a number");
	goto create_from_array_cleanup;
      }

      po = PyList_GetItem(n, i);
      if (PyNumber_Check(po)) {
	d->d->points[i].n = PyFloat_AsDouble(po);

	if (d->d->points[i].n <= 0.0) {
	  PyErr_SetString(PyExc_ValueError,
			  "All values in the n array must be greater than zero\n");
	  goto create_from_array_cleanup;
	}

      } else {
	PyErr_SetString(PyExc_ValueError,
			"Invalid element in n list, must be a number");
	goto create_from_array_cleanup;
      }

      if (d->d->points[i].x < d->d->xmin) {
	d->d->xmin = d->d->points[i].x;
      }
      if (d->d->points[i].x > d->d->xmax) {
	d->d->xmax = d->d->points[i].x;
      }

      if (d->d->points[i].y < d->d->ymin) {
	d->d->ymin = d->d->points[i].y;
      }
      if (d->d->points[i].y > d->d->ymax) {
	d->d->ymax = d->d->points[i].y;
      }

    }

    return d;

  create_from_array_cleanup:
    dataset1d_destroy(d->d);
    free(d);

    return NULL;
  }

  void set_xrange(double xmin, double xmax)
  {
    
    $self->d->xmin = xmin;
    $self->d->xmax = xmax;
  }

  double get_xmin(void)
  {
    return $self->d->xmin;
  }

  double get_xmax(void)
  {
    return $self->d->xmax;
  }

  void set_yrange(double ymin, double ymax)
  {
    $self->d->ymin = ymin;
    $self->d->ymax = ymax;
  }

  double get_ymin(void)
  {
    return $self->d->ymin;
  }

  double get_ymax(void)
  {
    return $self->d->ymax;
  }

  void set_lambda_std(double std)
  {
    $self->d->lambdastd = std;
  }

  double get_lambda_std(void)
  {
    return $self->d->lambdastd;
  }

  void set_lambda_range(double lambdamin, double lambdamax)
  {
    $self->d->lambdamin = lambdamin;
    $self->d->lambdamax = lambdamax;
  }

  double get_lambda_min(void)
  {
    return $self->d->lambdamin;
  }

  double get_lambda_max(void)
  {
    return $self->d->lambdamax;
  }

  ~dataset1d()
  {
    dataset1d_destroy($self->d);
    free($self);
  }
};

typedef struct {
  resultset1d_t *r;
} resultset1d;

%extend resultset1d {

  ~resultset1d()
  {
    resultset1d_destroy($self->r);
    free($self);
  }

  PyObject *proposed(void)
  {
    int na;
    const int *a;
    PyObject *r;

    a = resultset1d_get_propose($self->r, &na);
    
    r = pyrjmcmc_make_int_list(a, na);

    return r;
  }

  PyObject *acceptance(void)
  {
    int na;
    const int *a;
    int i;
    PyObject *r;

    a = resultset1d_get_accept($self->r, &na);
    r = PyList_New(na);
    
    for (i = 0; i < na; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    return r;
  }

  PyObject *partitions(void)
  {
    const int *a;
    int n;
    int i;
    PyObject *r;

    a = resultset1d_get_partitions($self->r);

    if (resultset1d_get_max_partitions($self->r) == 0 ||
	a == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_total($self->r);
  
    r = PyList_New(n);

    for (i = 0; i < n; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    return r;
  }

  PyObject *order_histogram(void)
  {
    const int *a;
    int i;
    PyObject *r;

    int mp;

    mp = resultset1d_get_max_order($self->r);
    a = resultset1d_get_order($self->r);

    if (mp == 0 ||
	a == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    r = PyList_New(mp + 1);
    for (i = 0; i <= mp; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    Py_INCREF(r);
    return r;
  }    
    
  PyObject *partition_histogram(void)
  {
    const int *a;
    int n;
    int i;
    PyObject *r;

    int mp;
    int *h;

    mp = resultset1d_get_max_partitions($self->r);
    a = resultset1d_get_partitions($self->r);

    if (mp == 0 ||
	a == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    h = malloc(sizeof(int) * (mp + 1));
    for (i = 0; i <= mp; i ++) {
      h[i] = 0;
    }

    for (i = 0; i < n; i ++) {
      h[a[i]] ++;
    }

    r = PyList_New(mp + 1);
    for (i = 0; i <= mp; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    free(h);

    return r;
  }

  PyObject *partition_location_histogram(void)
  {
    const int *a;
    int n;
    int i;
    PyObject *r;
    

    n = resultset1d_get_xsamples($self->r);
    a = resultset1d_get_partition_x_histogram($self->r);

    if (n == 0 || a == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    r = PyList_New(n);
    for (i = 0; i < n; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    return r;
  }

  PyObject *x(void)
  {
    PyObject *r;
    double *xc;
    int n;

    n = resultset1d_get_xsamples($self->r);
    if (n <= 0) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    xc = malloc(sizeof(double) * n);
    resultset1d_fill_xcoord_vector($self->r, xc);

    r = pyrjmcmc_make_float_list(xc, n);

    free(xc);

    Py_INCREF(r);
    return r;
  }

  PyObject *y(void)
  {
    PyObject *r;
    double *yc;
    int n;

    n = resultset1d_get_ysamples($self->r);
    if (n <= 0) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    yc = malloc(sizeof(double) * n);
    resultset1d_fill_ycoord_vector($self->r, yc);

    r = pyrjmcmc_make_float_list(yc, n);

    free(yc);

    Py_INCREF(r);
    return r;
  }

  PyObject *mean(void)
  {
    const double *m = resultset1d_get_mean($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    Py_INCREF(r);
    return r;
  }

  PyObject *median(void)
  {
    const double *m = resultset1d_get_median($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *mode(void)
  {
    const double *m = resultset1d_get_mode($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *credible_min(void)
  {
    const double *m = resultset1d_get_credible_min($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *credible_max(void)
  {
    const double *m = resultset1d_get_credible_max($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *misfit(void) 
  {
    const double *m = resultset1d_get_misfit($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_total($self->r);

    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *lambda_history(void)
  {
    const double *m = resultset1d_get_lambda($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1d_get_total($self->r);

    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *histogram(void)
  {
    const int **h = resultset1d_get_histogram($self->r);
    PyObject *r;
    int xs;
    int ys;
    
    if (h == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    xs = resultset1d_get_xsamples($self->r);
    ys = resultset1d_get_ysamples($self->r);

    r = pyrjmcmc_make_int_list_2d(h, xs, ys);

    return r;
  }

};

typedef struct {
  resultset1dfm_t *r;
} resultset1dfm;

%extend resultset1dfm {

  ~resultset1dfm()
  {
    resultset1dfm_destroy($self->r);
    free($self);
  }

  PyObject *proposed(void)
  {
    int na;
    const int *a;
    PyObject *r;

    a = resultset1dfm_get_propose($self->r, &na);
    
    r = pyrjmcmc_make_int_list(a, na);

    return r;
  }

  PyObject *acceptance(void)
  {
    int na;
    const int *a;
    int i;
    PyObject *r;

    a = resultset1dfm_get_accept($self->r, &na);
    r = PyList_New(na);
    
    for (i = 0; i < na; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    return r;
  }

  PyObject *partitions(void)
  {
    const int *a;
    int n;
    int i;
    PyObject *r;

    a = resultset1dfm_get_partitions($self->r);

    if (resultset1dfm_get_max_partitions($self->r) == 0 ||
	a == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_total($self->r);
  
    r = PyList_New(n);

    for (i = 0; i < n; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    return r;
  }

  PyObject *partition_histogram(void)
  {
    const int *a;
    int n;
    int i;
    PyObject *r;

    int mp;
    int *h;

    mp = resultset1dfm_get_max_partitions($self->r);
    a = resultset1dfm_get_partitions($self->r);

    if (mp == 0 ||
	a == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    h = malloc(sizeof(int) * (mp + 1));
    for (i = 0; i <= mp; i ++) {
      h[i] = 0;
    }

    for (i = 0; i < n; i ++) {
      h[a[i]] ++;
    }

    r = PyList_New(mp + 1);
    for (i = 0; i <= mp; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    free(h);

    return r;
  }

  PyObject *partition_location_histogram(void)
  {
    const int *a;
    int n;
    int i;
    PyObject *r;
    

    n = resultset1dfm_get_xsamples($self->r);
    a = resultset1dfm_get_partition_x_histogram($self->r);

    if (n == 0 || a == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    r = PyList_New(n);
    for (i = 0; i < n; i ++) {
      PyList_SetItem(r, i, PyInt_FromLong(a[i]));
    }

    return r;
  }

  PyObject *x(void)
  {
    PyObject *r;
    double *xc;
    int n;
    int t;

    n = resultset1dfm_get_xsamples($self->r);
    if (n <= 0) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    xc = malloc(sizeof(double) * n);
    t = n;
    resultset1dfm_fill_xcoord_vector($self->r, xc, &t);

    r = pyrjmcmc_make_float_list(xc, n);

    free(xc);

    Py_INCREF(r);
    return r;
  }

  PyObject *mean(int li = 0)
  {
    const double *m = resultset1dfm_get_local_parameter_mean($self->r, li);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    Py_INCREF(r);
    return r;
  }

  PyObject *median(int li = 0)
  {
    const double *m = resultset1dfm_get_local_parameter_median($self->r, li);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *mode(int li = 0)
  {
    const double *m = resultset1dfm_get_local_parameter_mode($self->r, li);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *credible_min(int li = 0)
  {
    const double *m = 
      resultset1dfm_get_local_parameter_credible_min($self->r, li);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *credible_max(int li = 0)
  {
    const double *m = 
      resultset1dfm_get_local_parameter_credible_max($self->r, li);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_xsamples($self->r);
    
    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *global_parameter(int gi) 
  {
    const double *m = resultset1dfm_get_global_parameter($self->r, gi);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_total($self->r);

    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

  PyObject *misfit(void) 
  {
    const double *m = resultset1dfm_get_misfit($self->r);
    PyObject *r;
    int n;

    if (m == NULL) {
      Py_INCREF(Py_None);
      return Py_None;
    }

    n = resultset1dfm_get_total($self->r);

    r = pyrjmcmc_make_float_list(m, n);

    return r;
  }

};

%feature("autodoc", "1");

/*
 * Single Partition 
 */

resultset1d *
regression_single1d(dataset1d *dataset,
		    int burnin = 10000,
		    int total = 50000,
		    int max_order = 5,
		    int xsamples = 100,
		    int ysamples = 100,
		    double credible_interval = 0.95);

resultset1d *
regression_single1d_sampled(dataset1d *dataset,
			    PyObject *callback,
			    int burnin = 10000,
			    int total = 50000,
			    int max_order = 5,
			    int xsamples = 100,
			    int ysamples = 100,
			    double credible_interval = 0.95);


/*
 * 1D Partition Models 
 */

resultset1d *
regression_part1d_zero(dataset1d *dataset,
		       double pd,
		       int burnin = 10000,
		       int total = 50000,
		       int max_partitions = 20,
		       int xsamples = 100,
		       int ysamples = 100,
		       double credible_interval = 0.95);

resultset1d *
regression_part1d_natural(dataset1d *dataset,
			  double pv,
			  double pd,
			  int burnin = 10000,
			  int total = 50000,
			  int max_partitions = 20,
			  int xsamples = 100,
			  int ysamples = 100,
			  double credible_interval = 0.95);

resultset1d *
regression_part1d(dataset1d *dataset,
		  double pd,
		  int burnin = 10000,
		  int total = 50000,
		  int max_partitions = 20,
		  int max_order = 5,
		  int xsamples = 100,
		  int ysamples = 100,
		  double credible_interval = 0.95);

resultset1d *
regression_part1d_sampled(dataset1d *dataset,
			  PyObject *callback,
			  double pd,
			  int burnin = 10000,
			  int total = 50000,
			  int max_partitions = 20,
			  int max_order = 5,
			  int xsamples = 100,
			  int ysamples = 100,
			  double credible_interval = 0.95);

resultset1dfm *
forwardmodel_part1d(PyObject *local_parameters,
		    PyObject *global_parameters,
		    PyObject *loglikelihood_cb,
		    double minx,
		    double maxx,
		    double pd,
		    int burnin = 10000,
		    int total = 50000,
		    int max_partitions = 20,
		    int xsamples = 100,
		    int ysamples = 100,
		    double credible_interval = 0.95);
		    

		    
		    
