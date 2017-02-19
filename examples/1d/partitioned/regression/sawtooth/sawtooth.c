
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <rjmcmc/regression.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

int main(int argc, char *argv[]) 
{
  dataset1d_t *data;
  resultset1d_t *results;

  int burnin = 1000;
  int total = 100000;
  int min_part = 2;
  int max_part = 50;
  int max_order = 5;
  int xsamples = 100;
  int ysamples = 500;
  double confidence = 0.95;
  double pd = 0.5;

  int nproc;
  const double *v;
  const int *iv;

  double sigma = 10.0;
  int i;

  double *xcoords;
  int xcl;

  data = dataset1d_load_fixed("data.txt", sigma);
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
  data->ymin = -100.0;
  data->ymax = 100.0;

  results = part1d_regression(data,
			      burnin,
			      total,
			      min_part,
			      max_part,
			      max_order,
			      xsamples,
			      ysamples,
			      confidence,
			      pd,
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

  xcl = xsamples;
  xcoords = rjmcmc_create_array_1d(xcl);
  if (xcoords == NULL) {
    fprintf(stderr, "error: failed to create array for xsamples\n");
    return -1;
  }
  resultset1d_fill_xcoord_vector(results, xcoords);

  v = resultset1d_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("sawtooth.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  v = resultset1d_get_mean(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_coords("sawtooth.mean", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  v = resultset1d_get_median(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_coords("sawtooth.median", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  v = resultset1d_get_mode(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_coords("sawtooth.mode", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  v = resultset1d_get_credible_min(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_coords("sawtooth.credible_min", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  v = resultset1d_get_credible_max(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_coords("sawtooth.credible_max", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save credible_max data\n");
    return -1;
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

  rjmcmc_destroy_array_1d(xcoords);
  resultset1d_destroy(results);
  
  return 0;
}

