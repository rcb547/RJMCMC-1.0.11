
#include <float.h>

#include "rjmcmc_helper.h"

#include <rjmcmc/rjmcmc_util.h>

static const int DEFAULT_RESULTS = 
  RESULTSET1D_MEAN |
  RESULTSET1D_MODE |
  RESULTSET1D_MEDIAN |
  RESULTSET1D_CREDIBLE;

static const int DEFAULT_FM_RESULTS = 
  RESULTSET1DFM_MEAN | 
  RESULTSET1DFM_MODE | 
  RESULTSET1DFM_MEDIAN | 
  RESULTSET1DFM_CREDIBLE;

static double pyrjmcmc_fm_part1d_likelihood(void *user_arg,
					    int npartitions,
					    const double *partition_boundaries,
					    int nglobalparameters,
					    const double *global_parameters,
					    part1d_fm_likelihood_state_t *state,
					    part1d_fm_value_at_t value_at,
					    part1d_fm_value_at_t gradient_at);


struct pyrjmcmc_fm_part1d_cbdata {
  PyObject *cb;
  int nlocalparameters;
};

resultset1d *
regression_single1d(dataset1d *dataset,
		    int burnin,
		    int total,
		    int max_order,
		    int xsamples,
		    int ysamples,
		    double credible_interval)
{
  resultset1d *R;
  resultset1d_t *r = single1d_regression(dataset->d,
					 burnin,
					 total,
					 max_order,
					 xsamples,
					 ysamples,
					 credible_interval,
					 rjmcmc_uniform,
					 rjmcmc_normal,
					 DEFAULT_RESULTS,
					 NULL,
					 NULL);

  if (r == NULL) {
    return NULL;
  }

  R = malloc(sizeof(resultset1d));
  R->r = r;

  return R;
}

struct sampled_cb_data {
  PyObject *cb;

  int xsamples;
  double *x;
  double *y;
  
  int error_flag;
};

static void
sampled_cb(void *state,
	   double *boundaries,
	   int nboundaries,
	   regression1d_value_at_t value_at,
	   double lambda,
	   void *user_arg)
{
  struct sampled_cb_data *d = 
    (struct sampled_cb_data *)user_arg;

  PyObject *xlist;
  PyObject *ylist;
  
  PyObject *arglist;
  PyObject *result;

  int i;

  if (d->cb != NULL && PyCallable_Check(d->cb)) {

    /*
     * Evaluate the current model
     */
    xlist = PyList_New(d->xsamples);
    ylist = PyList_New(d->xsamples);
    
    for (i = 0; i < d->xsamples; i ++) {
      d->y[i] = value_at(state, d->x[i]);

      PyList_SetItem(xlist, i, PyFloat_FromDouble(d->x[i]));
      PyList_SetItem(ylist, i, PyFloat_FromDouble(d->y[i]));
    }
  
    arglist = Py_BuildValue("(OO)", xlist, ylist);
    result = PyObject_CallObject(d->cb, arglist);
    Py_DECREF(arglist);

    if (result == NULL) {
      /* Error */
      d->cb = NULL;
      d->error_flag = -1;
    }
  }
}


resultset1d *
regression_single1d_sampled(dataset1d *dataset,
			    PyObject *callback,
			    int burnin,
			    int total,
			    int max_order,
			    int xsamples,
			    int ysamples,
			    double credible_interval)
{
  resultset1d *R;
  struct sampled_cb_data d;
  resultset1d_t *r;

  d.cb = callback;

  d.xsamples = xsamples;

  d.x = malloc(sizeof(double) * xsamples);
  rjmcmc_fill_coord_vector(dataset->d->xmin,
			   dataset->d->xmax,
			   xsamples,
			   d.x);
  
  d.y = malloc(sizeof(double) * ysamples);

  d.error_flag = 0;

  r = single1d_regression(dataset->d,
			  burnin,
			  total,
			  max_order,
			  xsamples,
			  ysamples,
			  credible_interval,
			  rjmcmc_uniform,
			  rjmcmc_normal,
			  DEFAULT_RESULTS,
			  sampled_cb,
			  &d);

  free(d.x);
  free(d.y);

  if (r == NULL) {
    return NULL;
  }

  if (d.error_flag) {
    return NULL;
  }

  R = malloc(sizeof(resultset1d));
  R->r = r;

  return R;
}

resultset1d *
regression_part1d_zero(dataset1d *dataset,
		       double pd,
		       int burnin,
		       int total,
		       int max_partitions,
		       int xsamples,
		       int ysamples,
		       double credible_interval)
{
  resultset1d *R;

  resultset1d_t *r = part1d_zero_regression(dataset->d,
					    burnin,
					    total,
					    2,
					    max_partitions,
					    xsamples,
					    ysamples,
					    credible_interval,
					    pd,
					    rjmcmc_uniform,
					    rjmcmc_normal,
					    DEFAULT_RESULTS,
					    NULL,
					    NULL);
  if (r == NULL) {
    return NULL;
  }

  R = malloc(sizeof(resultset1d));
  R->r = r;

  return R;
}

resultset1d *
regression_part1d_natural(dataset1d *dataset,
			  double pv,
			  double pd,
			  int burnin,
			  int total,
			  int max_partitions,
			  int xsamples,
			  int ysamples,
			  double credible_interval)
{
  resultset1d *R;

  resultset1d_t *r = part1d_natural_regression(dataset->d,
					       burnin,
					       total,
					       2,
					       max_partitions,
					       xsamples,
					       ysamples,
					       credible_interval,
					       pv,
					       pd,
					       rjmcmc_uniform,
					       rjmcmc_normal,
					       DEFAULT_RESULTS,
					       NULL,
					       NULL);
  if (r == NULL) {
    return NULL;
  }

  R = malloc(sizeof(resultset1d));
  R->r = r;

  return R;
}

resultset1d *
regression_part1d(dataset1d *dataset,
		  double pd,
		  int burnin,
		  int total,
		  int max_partitions,
		  int max_order,
		  int xsamples,
		  int ysamples,
		  double credible_interval)
{
  resultset1d *R;
  resultset1d_t *r;

  r = part1d_regression(dataset->d,
			burnin,
			total,
			2, /* min_partitions */
			max_partitions,
			max_order,
			xsamples,
			ysamples,
			credible_interval,
			pd,
			rjmcmc_uniform,
			rjmcmc_normal,
			DEFAULT_RESULTS,
			NULL, /* callback */
			NULL  /* user_arg */);
  if (r == NULL) {
    return NULL;
  }

  R = malloc(sizeof(resultset1d));
  R->r = r;

  return R;
}

resultset1d *
regression_part1d_sampled(dataset1d *dataset,
			  PyObject *callback,
			  double pd,
			  int burnin,
			  int total,
			  int max_partitions,
			  int max_order,
			  int xsamples,
			  int ysamples,
			  double credible_interval)
{
  resultset1d *R;
  struct sampled_cb_data d;
  resultset1d_t *r;

  d.cb = callback;

  d.xsamples = xsamples;

  d.x = malloc(sizeof(double) * xsamples);
  rjmcmc_fill_coord_vector(dataset->d->xmin,
			   dataset->d->xmax,
			   xsamples,
			   d.x);
  
  d.y = malloc(sizeof(double) * ysamples);

  d.error_flag = 0;

  r = part1d_regression(dataset->d,
			burnin,
			total,
			2, /* min_partitions */
			max_partitions,
			max_order,
			xsamples,
			ysamples,
			credible_interval,
			pd,
			rjmcmc_uniform,
			rjmcmc_normal,
			DEFAULT_RESULTS,
			sampled_cb,
			&d);

  free(d.x);
  free(d.y);

  if (r == NULL) {
    return NULL;
  }

  if (d.error_flag) {
    return NULL;
  }

  R = malloc(sizeof(resultset1d));
  R->r = r;

  return R;
}

resultset1dfm *
forwardmodel_part1d(PyObject *local_parameters,
		    PyObject *global_parameters,
		    PyObject *loglikelihood_cb,
		    double minx,
		    double maxx,
		    double pd,
		    int burnin,
		    int total,
		    int max_partitions,
		    int xsamples,
		    int ysamples,
		    double credible_interval)
{
  resultset1dfm *R;
  resultset1dfm_t *r;

  int nglobal_parameters;
  forwardmodelparameter_t *fglobal_parameters;

  int nlocal_parameters;
  forwardmodelparameter_t *flocal_parameters;
  
  struct pyrjmcmc_fm_part1d_cbdata cbdata;

  if (!PyCallable_Check(loglikelihood_cb)) {
    return NULL;
  }

  if (pyrjmcmc_list_to_parameter_list(global_parameters,
				      &fglobal_parameters,
				      &nglobal_parameters) < 0) {
    return NULL;
  }

  if (pyrjmcmc_list_to_parameter_list(local_parameters,
				      &flocal_parameters,
				      &nlocal_parameters) < 0) {
    return NULL;
  }

  cbdata.cb = loglikelihood_cb;
  Py_INCREF(cbdata.cb);
  cbdata.nlocalparameters = nlocal_parameters;

  r = part1d_forwardmodel(burnin,
			  total,
			  2,
			  max_partitions,
			  minx,
			  maxx,
			  xsamples,
			  ysamples,
			  credible_interval,
			  pd,
			  rjmcmc_uniform,
			  rjmcmc_normal,
			  nglobal_parameters,
			  fglobal_parameters,
			  nlocal_parameters,
			  flocal_parameters,
			  pyrjmcmc_fm_part1d_likelihood,
			  &cbdata,
			  DEFAULT_FM_RESULTS);

  Py_DECREF(cbdata.cb);

  if (r == NULL) {
    return NULL;
  }

  R = malloc(sizeof(resultset1dfm));
  R->r = r;

  return R;
}


PyObject *pyrjmcmc_make_float_list(const double *data, int size)
{
  PyObject *r;
  int i;

  r = PyList_New(size);

  for (i = 0; i < size; i ++) {
    PyList_SetItem(r, i, PyFloat_FromDouble(data[i]));
  }

  return r;
}

PyObject *pyrjmcmc_make_int_list(const int *data, int size)
{
  PyObject *r;
  int i;

  r = PyList_New(size);

  for (i = 0; i < size; i ++) {
    PyList_SetItem(r, i, PyInt_FromLong(data[i]));
  }

  return r;
}

PyObject *pyrjmcmc_make_int_list_2d(const int **data, int xsize, int ysize)
{
  PyObject *r;
  int i;

  r = PyList_New(xsize);
  for (i = 0; i < xsize; i ++) {
    PyList_SetItem(r, i, pyrjmcmc_make_int_list(data[i], ysize));
  }

  return r;
}

double pyrjmcmc_value(PyObject *o)
{
  if (o == NULL) {
    return 0.0;
  }

  if (PyFloat_Check(o)) {
    return PyFloat_AsDouble(o);
  }

  if (PyInt_Check(o)) {
    return (double)PyInt_AsLong(o);
  }

  if (PyLong_Check(o)) {
    return (double)PyLong_AsLong(o);
  }

  return 0.0;
}

int
pyrjmcmc_list_to_parameter_list(PyObject *l,
				forwardmodelparameter_t **parameters,
				int *nparameters)
{
  Py_ssize_t size;
  Py_ssize_t i;

  PyObject *el;
  PyObject *tl;

  forwardmodelparameter_t *f;

  if (l == Py_None) {
    *parameters = NULL;
    *nparameters = 0;
    return 0;
  }

  if (!PyList_Check(l)) {
    return -1;
  }

  size = PyList_Size(l);
  if (size == 0) {
    *parameters = NULL;
    *nparameters = 0;
    return 0;
  }    
  
  f = forwardmodelparameter_create(size);
  if (f == NULL) {
    return -1;
  }

  for (i = 0; i < size; i ++) {
    
    el = PyList_GetItem(l, i);
    if (el == NULL) {
      return -1;
    }
    
    if (!PyTuple_Check(el)) {
      return -1;
    }

    if (PyTuple_Size(el) != 3) {
      return -1;
    }

    tl = PyTuple_GetItem(el, 0);
    f[i].fmin = pyrjmcmc_value(tl);

    tl = PyTuple_GetItem(el, 1);
    f[i].fmax = pyrjmcmc_value(tl);

    tl = PyTuple_GetItem(el, 2);
    f[i].fstd_value = pyrjmcmc_value(tl);
    f[i].fstd_bd = pyrjmcmc_value(tl);

  }
  
  *nparameters = size;
  *parameters = f;

  return 0;
}

static double pyrjmcmc_fm_part1d_likelihood(void *user_arg,
					    int npartitions,
					    const double *partition_boundaries,
					    int nglobalparameters,
					    const double *global_parameters,
					    part1d_fm_likelihood_state_t *state,
					    part1d_fm_value_at_t value_at,
					    part1d_fm_value_at_t gradient_at)
{
  struct pyrjmcmc_fm_part1d_cbdata *cbdata = 
    (struct pyrjmcmc_fm_part1d_cbdata*)user_arg;

  PyObject *globalvalues;
  PyObject *boundaries;
  PyObject *localvalues;

  Py_ssize_t i;

  double midx;
  const double *values;

  PyObject *arguments;
  PyObject *result;

  double r;

  if (cbdata->cb == NULL) {
    fprintf(stderr, "null cb\n");
    return 0.0;
  }

  if (nglobalparameters == 0) {
    globalvalues = Py_None;
  } else {
    globalvalues = 
      pyrjmcmc_make_float_list(global_parameters, nglobalparameters);
  }
  Py_INCREF(globalvalues);

  boundaries = pyrjmcmc_make_float_list(partition_boundaries, npartitions);
  Py_INCREF(boundaries);

  localvalues = PyList_New(npartitions - 1);

  for (i = 1; i < npartitions;  i++) {
    midx = (partition_boundaries[i] + partition_boundaries[i - 1])/2.0;

    values = value_at(state, midx);

    PyList_SetItem(localvalues, 
		   i - 1, 
		   pyrjmcmc_make_float_list(values, cbdata->nlocalparameters));
  }

  Py_INCREF(localvalues);

  arguments = Py_BuildValue("(OOO)", globalvalues, boundaries, localvalues);
  Py_INCREF(arguments);

  result = PyObject_CallObject(cbdata->cb, arguments);

  Py_DECREF(arguments);
  Py_DECREF(localvalues);
  Py_DECREF(boundaries);
  Py_DECREF(globalvalues);
  
  if (result == NULL) {
    PySys_WriteStderr("rjmcmc: error in forward model callback\n");
    if (PyErr_Occurred() != NULL) {
      PyErr_Print();
    }
    return 0.0;
  }

  if (PyFloat_Check(result)) {
    r = PyFloat_AsDouble(result);
    Py_DECREF(result);
    return r;
  }

  PySys_WriteStderr("rjmcmc: error, callback did not return a float\n");
  Py_DECREF(result);
  return 0.0;
}

