
#include <stdio.h>
#include <math.h>

#include <rjmcmc/regression.h>
#include <rjmcmc/rjmcmc_util.h>

/*
 * Step 1 : Create a random number generator. You can either use the 
 * functions included in the rjmcmc library which use the underlying
 * system random number generators or calls to a mathematical library
 * such as the GNU Scientific library.
 */

/* Option (a) : Using rjmcmc */

static double unif(void)
{
  return rjmcmc_uniform();
}

static double norm(void)
{
  return rjmcmc_normal();
}

/* Option (b) : Using GSL*/

/*
#include <rjmcmc/rjmcmc_random.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
static gsl_rng *rng;

static double unif(void)
{
  return gsl_rng_uniform(rng);
}

static double norm(void)
{
  return gsl_ran_gaussian(rng, 1.0);
}
*/

int main(int argc, char *argv[]) 
{
  dataset1d_t *data;
  resultset1d_t *results;

  /*
   * Step 2: Set the simulation parameters
   */

  int burnin = 1000;               /* Number of iterations to be thrown away */
  int total = 50000;               /* Total iterations */
  int min_part = 2;                /* Minimum number of partitions */
  int max_part = 20;               /* Maximum number of partitions */
  int max_order = 5;               /* Maximum polynomial order */
  int xsamples = 100;              /* Number of points to sample the curve */
  int ysamples = 500;              /* Number of points to sample medians etc */
  double credible_interval = 0.95; /* Credible_Interval interval */
  double pv = 5.0;                 /* Std dev. of value changes */
  double pd = 10.0;	           /* Std dev. of move changes */		

  int nproc;
  const double *v;
  const int *iv;
  int i;

  /* 
   * Step 3: Initialize the random number generators.
   */ 
  
  /* Option (a): */
  rjmcmc_seed(1232131);

  /* Option (b): */
  /*
  rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, 1232131);
  */

  /* Step 4: Load the data */

  /* Option (a): Known errors - the file loaded should have 3 values
   * per line x, y, sigma
   */
  data = dataset1d_load_known("data.txt");
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to load data, "
	    "has it been created with the script?\n");
    return -1;
  }

  /* Option (b): Estimated/Fixed errors - the file loaded should have 2 values
   * per line x, y
   */
  /*
  data = dataset1d_load_fixed("data.txt", 5.0);
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to load data\n");
    return -1;
  }*/

  /*
   * Step 5: Set the data bounds
   */
  data->xmin = 0.0;
  data->xmax = 100.0;
  data->ymin = -50.0;
  data->ymax = 150.0;
  
  /*
   * Step 6: Run the simulation 
   */
  
  /* Option (a): Up to max_order polynomial fitting */

  results = part1d_regression(data,
  			      burnin,
  			      total,
  			      min_part,
  			      max_part,
  			      max_order,
  			      xsamples,
  			      ysamples,
  			      credible_interval,
  			      pd,
  			      unif,
  			      norm,
  			      RESULTSET1D_MEAN,
  			      NULL,
  			      NULL);

  /* Option (b): Zeroth order polynomials */

  /*
  results = part1d_zero_regression(data,
  				   burnin,
  				   total,
  				   min_part,
  				   max_part,
  				   xsamples,
  				   ysamples,
  				   credible_interval,
  				   pd,
  				   unif,
  				   norm,
  				   RESULTSET1D_MEAN,
  				   NULL,
  				   NULL);
  */

  /* Option (d): Joined line segments */
  /*
  results = part1d_natural_regression(data,
				      burnin,
				      total,
				      min_part,
				      max_part,
				      xsamples,
				      ysamples,
				      credible_interval,
				      pv,
				      pd,
				      unif,
				      norm,
				      RESULTSET1D_MEAN,
				      NULL,
				      NULL);
  */
  
  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run regression\n");
    return -1;
  }

  /*
   * Step 7: Query results (get the mean curve)
   */
  v = resultset1d_get_mean(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_vector("part1dtemplate.mean", v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  /* Print out proposal/acceptance information */
  iv = resultset1d_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultset1d_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");
  
  /*
   * Step 8 : Clean up
   */
  dataset1d_destroy(data);
  resultset1d_destroy(results);


  return 0;
}

