
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultsetfm.h>

#include <rjmcmc/rjmcmc_util.h>

static const double real_sigma = 1.0;
static const double real_value = 0.0;

#define DATASIZE 100
static double data[DATASIZE];

static double unif(void)
{
  return rjmcmc_uniform();
}

static double norm(void)
{
  return rjmcmc_normal();
}

static double likelihood(void *user_arg,
			 int nvalues,
			 const double *values,
			 int hierarchical,
			 int nhierarchical_parameters,
			 const double *hierarchical_values,
			 double *logdetce);

int main(int argc, char *argv[]) 
{

  int burnin = 1000;
  int total = 20000;
  int samples = 100;
  double confidence_interval = 0.95;
  int requested_results = RESULTSETFM_MEAN;

  forwardmodelparameter_t parameter;
  forwardmodelparameter_t hierarchical_parameter;

  int nproc;
  const double *v;
  const int *iv;
  int i;

  resultsetfm_t *results;

  /*
   * Initialize the synthetic data from N(0, sigma);
   */
  for (i = 0; i < DATASIZE; i ++) {
    data[i] = rjmcmc_normal() * real_sigma;
  }

  /*
   * Initialize the search space for the parameters
   */
  
  parameter.fmin = -10.0;
  parameter.fmax = 10.0;
  parameter.fstd_value = 0.5;
  parameter.fstd_bd = 0.0;

  hierarchical_parameter.fmin = 0.1;
  hierarchical_parameter.fmax = 5.0;
  hierarchical_parameter.fstd_value = 0.1;
  hierarchical_parameter.fstd_bd = 0.0;

  /*
   * Run the forward model
   */
  results = single_forwardmodel_hierarchical(burnin,
					     total,
					     unif, /*rjmcmc_uniform,*/
					     norm, /*rjmcmc_normal,*/
					     1,
					     &parameter,
					     1,
					     &hierarchical_parameter,
					     likelihood,
					     NULL,
					     samples,
					     confidence_interval,
					     requested_results);

  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run functionfit\n");
    return -1;
  }

  v = resultsetfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("functionfit.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultsetfm_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultsetfm_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  printf("mean: %f\n", 
	 resultsetfm_get_parameter_mean(results, 0));

  v = resultsetfm_get_parameter_history(results, 0);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get parameter history\n");
    return -1;
  }

  if (rjmcmc_save_vector("single_forwardmodel_hierarchical_c.history",
			 v,
			 total) < 0) {
    fprintf(stderr, "error: failed to save parameter history\n");
    return -1;
  }

  v = resultsetfm_get_hierarchical_parameter_history(results, 0);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get hierarchical parameter history\n");
    return -1;
  }

  if (rjmcmc_save_vector("single_forwardmodel_hierarchical_c.sigma_history",
			 v,
			 total) < 0) {
    fprintf(stderr, "error: failed to save parameter history\n");
    return -1;
  }

  resultsetfm_destroy(results);

  return 0;
}

static double likelihood(void *user_arg,
			 int nvalues,
			 const double *values,
			 int hierarchical,
			 int nhierarchical_parameters,
			 const double *hierarchical_values,
			 double *logdetce)
{
  double dv;
  int i;
  double sum;

  double sigma = hierarchical_values[0];
  
  sum = 0.0;

  for (i = 0; i < DATASIZE; i ++) {
    
    dv = values[0] - data[i];

    sum += (dv*dv);

  }

  if (hierarchical) {

    *logdetce = 2.0 * (double)DATASIZE * log(sigma);

  }

  return sum/(sigma * sigma * 2.0);
}

