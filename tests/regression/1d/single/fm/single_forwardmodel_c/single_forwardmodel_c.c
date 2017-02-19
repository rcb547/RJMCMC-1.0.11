
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultsetfm.h>

#include <rjmcmc/rjmcmc_util.h>

static const double sigma = 1.0;
static const double real_value = 0.0;

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
			 const double *values);

int main(int argc, char *argv[]) 
{

  int burnin = 1000;
  int total = 20000;
  int samples = 100;
  double confidence_interval = 0.95;
  int requested_results = RESULTSETFM_MEAN;

  forwardmodelparameter_t parameter;

  int nproc;
  const double *v;
  const int *iv;
  int i;

  resultsetfm_t *results;

  /*
   * Initialize the search space for the parameters
   */
  
  parameter.fmin = -10.0;
  parameter.fmax = 10.0;
  parameter.fstd_value = 0.5;
  parameter.fstd_bd = 0.0;

  /*
   * Run the forward model
   */

  results = single_forwardmodel(burnin,
				total,
				unif, /*rjmcmc_uniform,*/
				norm, /*rjmcmc_normal,*/
				1,
				&parameter,
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

  if (rjmcmc_save_vector("single_forwardmodel_c.history",
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
			 const double *values)
{
  double dv;

  dv = values[0] - real_value;

  return (dv * dv)/(sigma * sigma * 2.0);
}

