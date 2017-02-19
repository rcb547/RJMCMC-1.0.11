
#ifndef rjmcmc_helper_h
#define rjmcmc_helper_h

#include <Python.h>

#include <stdlib.h>

#include <rjmcmc/regression.h>
#include <rjmcmc/resultset1d.h>
#include <rjmcmc/dataset1d.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/resultset1dfm.h>

#include <rjmcmc/rjmcmc_random.h>

#include <float.h>

typedef struct {
  dataset1d_t *d;
} dataset1d;

typedef struct {
  resultset1d_t *r;
} resultset1d;

typedef struct {
  resultset1dfm_t *r;
} resultset1dfm;

resultset1d *
regression_single1d(dataset1d *dataset,
		    int burnin,
		    int total,
		    int max_order,
		    int xsamples,
		    int ysamples,
		    double credible_interval);

resultset1d *
regression_single1d_sampled(dataset1d *dataset,
			    PyObject *callback,
			    int burnin,
			    int total,
			    int max_order,
			    int xsamples,
			    int ysamples,
			    double credible_interval);

resultset1d *
regression_part1d_zero(dataset1d *dataset,
		       double pd,
		       int burnin,
		       int total,
		       int max_partitions,
		       int xsamples,
		       int ysamples,
		       double credible_interval);

resultset1d *
regression_part1d_natural(dataset1d *dataset,
			  double pv,
			  double pd,
			  int burnin,
			  int total,
			  int max_partitions,
			  int xsamples,
			  int ysamples,
			  double credible_interval);

resultset1d *
regression_part1d(dataset1d *dataset,
		  double pd,
		  int burnin,
		  int total,
		  int max_partitions,
		  int max_order,
		  int xsamples,
		  int ysamples,
		  double credible_interval);

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
			  double credible_interval);

resultset1d *
regression_part1d_knot(dataset1d *dataset,
		       double pd,
		       int burnin,
		       int total,
		       int max_partitions,
		       int max_order,
		       int xsamples,
		       int ysamples,
		       double credible_interval);

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
		    double credible_interval);

PyObject *pyrjmcmc_make_float_list(const double *data, int size);
PyObject *pyrjmcmc_make_int_list(const int *data, int size);
PyObject *pyrjmcmc_make_int_list_2d(const int **data, int xsize, int ysize);

int
pyrjmcmc_list_to_parameter_list(PyObject *l,
				forwardmodelparameter_t **parameters,
				int *nparameters);

#endif /* rjmcmc_helper_h */
