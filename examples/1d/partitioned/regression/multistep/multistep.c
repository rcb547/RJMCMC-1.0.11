
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <rjmcmc/regression.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "nzO:h";

static struct option long_options[] = {
  {"zero", 0, 0, 'z'},
  {"natural", 0, 0, 'n'},
  {"max-order", 1, 0, 'O'},

  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[]) 
{
  int c;
  int option_index;

  dataset1d_t *data;
  resultset1d_t *results;

  int burnin = 1000;
  int total = 50000;
  int min_part = 2;
  int max_part = 100;
  int max_order = 5;
  int xsamples = 100;
  int ysamples = 500;
  double credible = 0.95;
  double pv = 5.0;
  double pd = 1.0;				

  int nproc;
  const double *v;
  const int *iv;
  int i;
  int natural = 0;
  int zero = 0;

  double *xcoords;
  int xcl;

  data = dataset1d_load_known("data.txt");
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to load data, "
	    "has it been created with the script?\n");
    return -1;
  }

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'n':
      natural = 1;
      break;

    case 'z':
      zero = 1;
      break;

    case 'O':
      max_order = atoi(optarg);
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
  data->xmin = 0.0;
  data->xmax = 100.0;
  data->ymin = -30.0;
  data->ymax = 30.0;
  
  if (natural) {
    results = part1d_natural_regression(data,
					burnin,
					total,
					min_part,
					max_part,
					xsamples,
					ysamples,
					credible,
					pv,
					pd,
					rjmcmc_uniform,
					rjmcmc_normal,
					RESULTSET1D_MEAN |
					RESULTSET1D_MEDIAN |
					RESULTSET1D_MODE |
					RESULTSET1D_CREDIBLE,
					NULL,
					NULL);
  } else if (zero) {
    results = part1d_zero_regression(data, 
				     burnin,
				     total,
				     min_part,
				     max_part,
				     xsamples,
				     ysamples,
				     credible,
				     pd,
				     rjmcmc_uniform,
				     rjmcmc_normal,
				     RESULTSET1D_MEAN |
				     RESULTSET1D_MEDIAN |
				     RESULTSET1D_MODE |
				     RESULTSET1D_CREDIBLE,
				     NULL,
				     NULL);
  } else {
    results = part1d_regression(data,
				burnin,
				total,
				min_part,
				max_part,
				max_order,
				xsamples,
				ysamples,
				credible,
				pd,
				rjmcmc_uniform,
				rjmcmc_normal,
				RESULTSET1D_MEAN |
				RESULTSET1D_MEDIAN |
				RESULTSET1D_MODE |
				RESULTSET1D_CREDIBLE,
				NULL,
				NULL);
  }

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
  if (rjmcmc_save_vector("multistep.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultset1d_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("multistep.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("multistep.partition_hist",
					  0,
					  max_part,
					  iv, 
					  total) < 0) {
    fprintf(stderr, "error: failed to save partition histogram.\n");
    return -1;
  }

  iv = resultset1d_get_partition_x_histogram(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partition histogram\n");
    return -1;
  }

  if (rjmcmc_save_int_coords("multistep.partition_x_hist",
			     xcoords,
			     iv,
			     xsamples) < 0) {
    fprintf(stderr, "error: failed to save partition x histogram\n");
    return -1;
  }

  v = resultset1d_get_mean(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_coords("multistep.mean", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  v = resultset1d_get_median(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_coords("multistep.median", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  v = resultset1d_get_mode(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_coords("multistep.mode", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  v = resultset1d_get_credible_min(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_coords("multistep.credible_min", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  v = resultset1d_get_credible_max(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_coords("multistep.credible_max", xcoords, v, xsamples) < 0) {
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

static void usage(const char *pname)
{
  fprintf(stderr, 
	  "usage: %s [options]\n"
	  "where options is on or more of:\n"
	  "\n"
	  " -n|--natural            use natural regression method\n"
	  " -z|--zero               use zero order hold method\n"
	  " -O|--max-order  <int>   maximum order for normal regression\n"
	  "\n"
	  " -h|--help               show usage information\n"
	  "\n",
	  pname);
}
