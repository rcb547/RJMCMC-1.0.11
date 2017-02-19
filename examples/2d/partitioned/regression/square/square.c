
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <rjmcmc/regression.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "lh";

static struct option long_options[] = {
  {"lambda", 0, 0, 'l'},
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[]) 
{
  int c;
  int option_index;
  
  dataset2d_t *data;
  resultset2d_t *results;

  int burnin = 1000;
  int total = 20000;
  int min_part = 2;
  int max_part = 1000;
  int xsamples = 100;
  int ysamples = 100;
  int zsamples = 200;
  double credible = 0.95;
  double pv = 5.0;
  double pd = 5.0;

  int nproc;
  const double *v;
  const double **m;
  const int *iv;
  const int **im;

  int i;
  
  double *xcoords;
  int xcl;

  double *ycoords;
  int ycl;

  double sigma = 10.0;
  int use_lambda = 0;
  
  double lambdamin = 0.5;
  double lambdamax = 3.0;
  double lambdastd = 0.05;

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'l':
      use_lambda = -1;
      break;

    default:
      fprintf(stderr, "error: invalid option\n");
      return -1;

    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  /*
   * Set the data bounds
   */

  /*
   * Set the sampling parameters for lambda if enabled.
   */
  if (use_lambda) {
    sigma = 5.0;
  } else {
    sigma = 10.0;
  }
  
  data = dataset2d_load_fixed("data.txt", sigma);
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to load data, "
	    "has it been created with the script?\n");
    return -1;
  }

  printf("Auto xrange: %f %f\n", data->xmin, data->xmax);
  printf("Auto yrange: %f %f\n", data->ymin, data->ymax);
  printf("Auto zrange: %f %f\n", data->zmin, data->zmax);
  data->xmin = -50.0;
  data->xmax = 50.0;
  data->ymin = -50.0;
  data->ymax = 50.0;
  data->zmin = -50.0;
  data->zmax = 50.0;
  printf("Manual xrange: %f %f\n", data->xmin, data->xmax);
  printf("Manual yrange: %f %f\n", data->ymin, data->ymax);
  printf("Manual zrange: %f %f\n", data->zmin, data->zmax);

  if (use_lambda) {
    data->lambdamin = lambdamin;
    data->lambdamax = lambdamax;
    data->lambdastd = lambdastd;
  } else {
    data->lambdastd = 0.0;
  }

  results = part2d_regression(data,
			      burnin,
			      total,
			      min_part,
			      max_part,
			      xsamples,
			      ysamples,
			      zsamples,
			      credible,
			      pv,
			      pd,
			      rjmcmc_uniform,
			      rjmcmc_normal,
			      RESULTSET2D_MEAN |
			      RESULTSET2D_MEDIAN |
			      RESULTSET2D_MODE |
			      RESULTSET2D_CREDIBLE,
			      NULL,
			      NULL);
  
  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run regression\n");
    return -1;
  }

  xcl = xsamples;
  xcoords = rjmcmc_create_array_1d(xcl);
  resultset2d_fill_xcoord_vector(results, xcoords, &xcl);
  
  ycl = ysamples;
  ycoords = rjmcmc_create_array_1d(ycl);
  resultset2d_fill_ycoord_vector(results, ycoords, &ycl);

  v = resultset2d_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("square.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  if (use_lambda) {
    v = resultset2d_get_lambda(results);
    if (v == NULL) {
      fprintf(stderr, "error: failed to get lambda data\n");
      return -1;
    }
    if (rjmcmc_save_vector("square.lambda", v, total) < 0) {
      fprintf(stderr, "error: failed to save lambda data\n");
      return -1;
    }
    if (rjmcmc_save_vector_as_histogram("square.lambda_hist", 
					data->lambdamin,
					data->lambdamax,
					xsamples,
					v, total) < 0) {
      fprintf(stderr, "error: failed to save lambda histogram\n");
      return -1;
    }
  }

  iv = resultset2d_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("square.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("square.partition_hist",
					  0,
					  max_part,
					  iv,
					  total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }

  m = resultset2d_get_mean(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("square.mean", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  m = resultset2d_get_median(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("square.median", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  m = resultset2d_get_mode(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("square.mode", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  m = resultset2d_get_credible_min(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("square.credible_min", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  m = resultset2d_get_credible_max(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("square.credible_max", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_max data\n");
    return -1;
  }

  im = resultset2d_get_centres(results);
  if (im == NULL) {
    fprintf(stderr, "error: failed to get centres\n");
    return -1;
  }
  if (rjmcmc_save_int_matrix("square.centres", im, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save centres\n");
    return -1;
  }
  

  iv = resultset2d_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultset2d_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  dataset2d_destroy(data);
  resultset2d_destroy(results);

  rjmcmc_destroy_array_1d(xcoords);
  rjmcmc_destroy_array_1d(ycoords);

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr, 
	  "usage: %s [options]\n"
	  "where options is on or more of:\n"
	  "\n"
	  " -h|--help               show usage information\n"
	  "\n",
	  pname);
}
