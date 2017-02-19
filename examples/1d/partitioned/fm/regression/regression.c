
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultset1dfm.h>

#include <rjmcmc/rjmcmc_util.h>

struct regression_data {
  int size;
  double *x;
  double *y;
  double sigma;
};

static struct regression_data *
create_data(int xsamples,
	    double xmin,
	    double xmax,
	    double ymin,
	    double ymax,
	    double n,
	    int nsteps);

static double 
regression_likelihood_hierarchical(void *user_arg,
				   int npartitions,
				   const double *partitions,
				   int nglobalparameters,
				   const double *global_parameters,
				   int hierarchical,
				   int nhierarchical,
				   const double *hierarchicalparameters,
				   part1d_fm_hierarchical_likelihood_state_t *state,
				   part1d_fm_hierarchical_value_at_t value_at,
				   part1d_fm_hierarchical_value_at_t gradient_at,
				   double *logdetce);
						 
static double 
regression_likelihood(void *user_arg,
		      int npartitions,
		      const double *partitions,
		      int nglobalparameters,
		      const double *global_parameters,
		      part1d_fm_likelihood_state_t *state,
		      part1d_fm_value_at_t value_at,
		      part1d_fm_value_at_t gradient_at);

int main(int argc, char *argv[]) 
{
  struct regression_data *data;
  resultset1dfm_t *results;

  int burnin = 1000;
  int total = 20000;
  int min_part = 2;
  int max_part = 10;
  double xmin = 0.0;
  double xmax = 100.0;
  int xsamples = 100;
  int ysamples = 200;
  double confidence = 0.95;
  double pd = 1.0;
  int hierarchical = 1;

  forwardmodelparameter_t local_parameter;
  forwardmodelparameter_t hierarchical_parameter;

  int nproc;
  const double *v;
  const int *iv;
  int i;

  double *xcoords;
  int xcl;

  double sigma_mean = 20.0;

  local_parameter.fmin = -150.0;
  local_parameter.fmax = 150.0;
  local_parameter.fstd_value = 15.0;
  local_parameter.fstd_bd = 15.0;

  hierarchical_parameter.fmin = sigma_mean - 10.0;
  hierarchical_parameter.fmax = sigma_mean + 10.0;
  hierarchical_parameter.fstd_value = 1.0;
  hierarchical_parameter.fstd_bd = 0.0;

  data = create_data(xsamples, xmin, xmax, -100.0, 100.0, sigma_mean, 5);
  if (data == NULL) {
    fprintf(stderr, 
	    "error: failed to create data\n");
    return -1;
  }

  if (rjmcmc_save_coords("data.txt", data->x, data->y, data->size) < 0) {
    fprintf(stderr,
	    "error: failed to save data\n");
    return -1;
  }

  if (hierarchical) {
    results = part1d_forwardmodel_hierarchical(burnin,
					       total,
					       min_part,
					       max_part,
					       xmin,
					       xmax,
					       xsamples,
					       ysamples,
					       confidence,
					       pd,
					       rjmcmc_uniform,
					       rjmcmc_normal,
					       0,
					       NULL,
					       1,
					       &local_parameter,
					       1,
					       &hierarchical_parameter,
					       regression_likelihood_hierarchical,
					       data,
					       RESULTSET1DFM_MEAN);
    
  } else {
    results = part1d_forwardmodel(burnin,
				  total,
				  min_part,
				  max_part,
				  xmin,
				  xmax,
				  xsamples,
				  ysamples,
				  confidence,
				  pd,
				  rjmcmc_uniform,
				  rjmcmc_normal,
				  0,
				  NULL,
				  1,
				  &local_parameter,
				  regression_likelihood,
				  data,
				  RESULTSET1DFM_MEAN);
  }

  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run regression\n");
    return -1;
  }

  if (rjmcmc_save_vector("regression.data", data->y, xsamples) < 0) {
    fprintf(stderr, "error: failed to save data\n");
    return -1;
  }

  xcl = xsamples;
  xcoords = rjmcmc_create_array_1d(xcl);
  if (xcoords == NULL) {
    fprintf(stderr, "error: failed to create array for xsamples\n");
    return -1;
  }
  resultset1dfm_fill_xcoord_vector(results, xcoords, &xcl);

  v = resultset1dfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("regression.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultset1dfm_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("regression.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("regression.partition_hist",
					  2,
					  max_part,
					  iv, total) < 0) {
    fprintf(stderr, "error: failed to save partition histogram data\n");
    return -1;
  }

  iv = resultset1dfm_get_partition_x_histogram(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partition histogram\n");
    return -1;
  }

  if (rjmcmc_save_int_coords("regression.partition_x_hist",
			     xcoords,
			     iv,
			     xsamples) < 0) {
    fprintf(stderr, "error: failed to save partition x histogram\n");
    return -1;
  }

  v = resultset1dfm_get_local_parameter_mean(results, 0);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
    if (rjmcmc_save_coords("regression.mean", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  iv = resultset1dfm_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultset1dfm_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  rjmcmc_destroy_array_1d(xcoords);
  resultset1dfm_destroy(results);

  rjmcmc_destroy_array_1d(data->x);
  rjmcmc_destroy_array_1d(data->y);
  free(data);

  return 0;
}

static struct regression_data *
create_data(int xsamples,
	    double xmin,
	    double xmax,
	    double ymin,
	    double ymax,
	    double n,
	    int nsteps)
{
  int i;
  int j;
  int jwrap;
  struct regression_data *d;
  double l;

  jwrap = xsamples/nsteps;
  if (jwrap <= 1) {
    fprintf(stderr, "create_data: too many steps\n");
    return NULL;
  }

  d = (struct regression_data*)malloc(sizeof(struct regression_data));
  if (d == NULL) {
    return NULL;
  }

  d->x = rjmcmc_create_array_1d(xsamples);
  if (d->x == NULL) {
    return NULL;
  }

  d->y = rjmcmc_create_array_1d(xsamples);
  if (d->y == NULL) {
    return NULL;
  }
  
  l = rjmcmc_random_choose_double(ymin, ymax, rjmcmc_uniform);

  d->size = xsamples;
  d->sigma = n + (n/2.0) * rjmcmc_normal();

  for (i = 0, j = 0; i < xsamples; i ++) {

    d->x[i] = ((double)i + 0.5)/(double)(xsamples) * (xmax - xmin) + xmin;
    d->y[i] = l + n * rjmcmc_normal();

    j ++;
    if (j > jwrap) {
      j = 0;
      l = rjmcmc_random_choose_double(ymin, ymax, rjmcmc_uniform);
    }
  }

  return d;
}

static double 
regression_likelihood_hierarchical(void *user_arg,
				   int npartitions,
				   const double *partitions,
				   int nglobalparameters,
				   const double *global_parameters,
				   int hierarchical,
				   int nhierarchical,
				   const double *hierarchicalparameters,
				   part1d_fm_hierarchical_likelihood_state_t *state,
				   part1d_fm_hierarchical_value_at_t value_at,
				   part1d_fm_hierarchical_value_at_t gradient_at,
				   double *logdetce)
{
  struct regression_data *data = (struct regression_data*)user_arg;
  int i;
  double sum;
  double dy;
  double sigma;
  double denom;

  const double *local_parameters;

  sum = 0;
  
  sigma = hierarchicalparameters[0];
  denom = 2.0 * sigma * sigma;

  for (i = 0; i < data->size; i ++) {

    local_parameters = value_at(state, data->x[i]);

    dy = local_parameters[0] - data->y[i];
    sum += (dy * dy)/denom;

  }

  if (hierarchical) {
    *logdetce = 2.0 * (double)(data->size) * log(sigma);
  }

  return sum;
}

static double 
regression_likelihood(void *user_arg,
		      int npartitions,
		      const double *partitions,
		      int nglobalparameters,
		      const double *global_parameters,
		      part1d_fm_likelihood_state_t *state,
		      part1d_fm_value_at_t value_at,
		      part1d_fm_value_at_t gradient_at)
{
  struct regression_data *data = (struct regression_data*)user_arg;
  double dummy;

  return regression_likelihood_hierarchical(user_arg,
					    npartitions,
					    partitions,
					    nglobalparameters,
					    global_parameters,
					    0, /* Don't calculate logdetce */
					    1, /* nhierarchical parameters */
					    &(data->sigma), /* pass our sigma */
					    (part1d_fm_hierarchical_likelihood_state_t *)state,
					    (part1d_fm_hierarchical_value_at_t)value_at,
					    (part1d_fm_hierarchical_value_at_t)gradient_at,
					    &dummy);
}

