
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <getopt.h>

#include <rjmcmc/rjmcmc_config.h>
#include <rjmcmc/regression.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static struct option long_options[] = {
  {"lambda", 0, 0, 'l'},

  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};
static char short_options[] = "lh";

static void usage(const char *pname)
{
  fprintf(stderr, 
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -l|--lambda   Also solve for scaling of error.\n"
	  "\n"
	  " -h|--help     Usage information\n"
	  "\n",
	  pname);
}

int main(int argc, char *argv[]) 
{
  int c;
  int option_index = 0;

  dataset1d_t *data;
  resultset1d_t *results;

  int burnin = 1000;
  int total = 50000;
  int max_order = 9;
  int xsamples = 100;
  int ysamples = 500;
  double credible = 0.95;
  double *xcoords;

  int nproc;
  const double *v;
  const int *iv;
  int i;

  int use_hierarchical = 0;

  while (1) {
    
    c = getopt_long(argc, argv, short_options, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 'l':
      use_hierarchical = -1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  data = dataset1d_load_known("data.txt");
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to load data, "
	    "has it been created with the script?\n");
    return -1;
  }

  /*
   * Set the data bounds
   */
  data->xmin = 0.0;
  data->xmax = 100.0;
  data->ymin = -50.0;
  data->ymax = 150.0;

  /*
   * Set the lambda parameters if hierarchical is desired.
   */
  if (use_hierarchical) {
    data->lambdamin = 0.5;
    data->lambdamax = 2.0;
    data->lambdastd = 0.3;
  }
  
  results = single1d_regression(data,
				burnin,
				total,
				max_order,
				xsamples,
				ysamples,
				credible,
				rjmcmc_uniform,
				rjmcmc_normal,
				RESULTSET1D_MEAN |
				RESULTSET1D_MEDIAN |
				RESULTSET1D_MODE |
				RESULTSET1D_CREDIBLE,
				NULL,
				NULL);

  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run regression\n");
    return -1;
  }

  xcoords = rjmcmc_create_array_1d(xsamples);
  resultset1d_fill_xcoord_vector(results, xcoords);

  v = resultset1d_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("cubic.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  v = resultset1d_get_mean(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_coords("cubic.mean", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  v = resultset1d_get_median(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_coords("cubic.median", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  v = resultset1d_get_mode(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_coords("cubic.mode", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  v = resultset1d_get_credible_min(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_coords("cubic.credible_min", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  v = resultset1d_get_credible_max(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_coords("cubic.credible_max", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save credible_max data\n");
    return -1;
  }

  iv = resultset1d_get_order(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get order data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("cubic.order", iv, max_order + 1) < 0) {
    fprintf(stderr, "error: failed to save order data\n");
    return -1;
  }

  if (use_hierarchical) {
    v = resultset1d_get_lambda(results);
    if (v == NULL) {
      fprintf(stderr, "error: failed to get lambda data\n");
      return -1;
    }
    if (rjmcmc_save_vector("cubic.lambda", v, total) < 0) {
      fprintf(stderr, "error: failed to save lambda data\n");
      return -1;
    }
  }

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

  resultset1d_destroy(results);

  dataset1d_destroy(data);

#if 0//defined(HAVE_GSL)
  gsl_rng_free(rng);
#endif   

  return 0;
}

