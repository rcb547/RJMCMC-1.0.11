
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <rjmcmc/regression.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "lsr:t:h";

static struct option long_options[] = {
  {"lambda", 0, 0, 'l'},
  {"random-seed", 1, 0, 'r'},
  {"total", 1, 0, 't'}, 
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static FILE* debug_fp = NULL;
void debug_boundaries(void *state,
		      double *boundaries,
		      int nboundaries,
		      regression1d_value_at_t value_at,
		      double lambda,
		      void *user_arg);


int main(int argc, char *argv[]) 
{
  int c;
  int option_index;

  dataset1d_t *data;
  resultset1d_t *results;

  int burnin = 1000;
  int total = 50000;
  int min_part = 2;
  int max_part = 20;
  int xsamples = 100;
  int ysamples = 500;
  double credible = 0.95;
  double pd = 1.0;

  int nproc;
  const double *v;
  const int *iv;
  int i;

  int use_lambda = 0;

  double *xc;
  int xcl;

  int seed = 0;

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
    case 'l':
      use_lambda = 1;
      break;
      
    case 'r':
      seed = atoi(optarg);
      break;

    case 't':
      total = atoi(optarg);
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

  if (use_lambda) {
    data->lambdamin = 0.5;
    data->lambdamax = 2.0;
    data->lambdastd = 0.1;
  }
  
  rjmcmc_seed(seed);

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
				   debug_boundaries,
				   NULL);

  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run regression\n");
    return -1;
  }

  xc = rjmcmc_create_array_1d(xsamples);
  xcl = xsamples;
  resultset1d_fill_xcoord_vector(results, xc);

  v = resultset1d_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("zeromultistep.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultset1d_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("zeromultistep.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("zeromultistep.parthist",
					  0,
					  max_part,
					  iv, 
					  total) < 0) {
    fprintf(stderr, "error: failed to save partition histogram.\n");
    return -1;
  }

  v = resultset1d_get_mean(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_coords("zeromultistep.mean", xc, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  v = resultset1d_get_median(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_coords("zeromultistep.median", xc, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  v = resultset1d_get_mode(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_coords("zeromultistep.mode", xc, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  v = resultset1d_get_credible_min(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_coords("zeromultistep.credible_min", xc, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  v = resultset1d_get_credible_max(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_coords("zeromultistep.credible_max", xc, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save credible_max data\n");
    return -1;
  }

  iv = resultset1d_get_partition_x_histogram(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partition x histogram\n");
    return -1;
  }
  if (rjmcmc_save_int_coords("zeromultistep.partitionx", xc, iv, xsamples) < 0) {
    fprintf(stderr, "error: failed to save partition x histogram\n");
    return -1;
  }
   

  if (use_lambda) {
    v = resultset1d_get_lambda(results);
    if (v == NULL) {
      fprintf(stderr, "error: failed to get lambda data\n");
      return -1;
    }

    if (rjmcmc_save_vector("zeromultistep.lambda", v, total) < 0) {
      fprintf(stderr, "error: failed to save lambda data\n");
      return -1;
    }

    if (rjmcmc_save_vector_as_histogram("zeromultistep.lambda_hist",
					data->lambdamin,
					data->lambdamax,
					xsamples,
					v,
					total) < 0) {
      fprintf(stderr, "error: failed to save the lambda histogram data\n");
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
  rjmcmc_destroy_array_1d(xc);
  
  if (debug_fp != NULL) {
    fclose(debug_fp);
  }

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr, 
	  "usage: %s [options]\n"
	  "where options is on or more of:\n"
	  "\n"
	  " -l|--lambda             Estimate lambda (error scale parameter)\n"
	  " -h|--help               show usage information\n"
	  "\n",
	  pname);
}

void debug_boundaries(void *state,
		      double *boundaries,
		      int nboundaries,
		      regression1d_value_at_t value_at,
		      double lambda,
		      void *user_arg)
{
  int i;

  if (debug_fp == NULL) {
    debug_fp = fopen("debug.txt", "w");
  }

  for (i = 0; i < nboundaries; i ++) {
    fprintf(debug_fp, "%f ", boundaries[i]);
  }
  fprintf(debug_fp, "\n");
}
