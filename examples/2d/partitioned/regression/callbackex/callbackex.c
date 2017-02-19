
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <rjmcmc/regression.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "h";

static struct option long_options[] = {
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

#define NXSAMPLES 25
#define NYSAMPLES 25

struct user_struct {
  double xc[NXSAMPLES];
  double yc[NXSAMPLES];
  double mean[NXSAMPLES][NYSAMPLES];
  int n;
};

static void regression2d_cb(void *state,
			    double *x,
			    double *y,
			    int ncells,
			    regression2d_value_at_t value_at,
			    double lambda,
			    void *user_arg);


int main(int argc, char *argv[]) 
{
  int c;
  int option_index;
  
  dataset2d_t *data;
  resultset2d_t *results;

  int burnin = 10000;
  int total = 50000;
  int min_part = 2;
  int max_part = 100;
  int xsamples = 100;
  int ysamples = 100;
  int zsamples = 200;
  double credible = 0.95;
  double pd = 1.0;
  double pv = 1.0;

  int nproc;
  const double *v;
  const double **m;
  const int *iv;
  const int **im;

  int i;
  int j;

  double *xcoords;
  int xcl;

  double *ycoords;
  int ycl;
  
  double sigma = 1.0;
  double lambda = 1.0;

  struct user_struct us;

  data = dataset2d_load_fixed("data.txt", sigma);
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to load data, "
	    "has it been created with the script?\n");
    return -1;
  }

  if (lambda == 0.0) {

    data->lambdamin = 0.5;
    data->lambdamax = 2.0;
    data->lambdastd = 0.1;
  }
  

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
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
  printf("Auto zrange: %f %f\n", data->zmin, data->zmax);
  data->xmin = -50.0;
  data->xmax = 50.0;
  data->ymin = -50.0;
  data->ymax = 50.0;
  data->zmin = -50.0;
  data->zmax = 50.0;

  /*
   * Initialize our user struct to zero at set the grid we 
   * want to sample.
   */
  us.n = 0;
  for (j = 0; j < NYSAMPLES; j ++) {
    for (i = 0; i < NXSAMPLES; i ++) {
      us.mean[i][j] = 0.0;
    }
  }

  for (i = 0; i < NYSAMPLES; i ++) {
    us.xc[i] = data->xmin + (data->xmax - data->xmin) * (double)i/(double)(NXSAMPLES - 1);
  }
  for (j = 0; j < NYSAMPLES; j ++) {
    us.yc[j] = data->ymin + (data->ymax - data->ymin) * (double)j/(double)(NYSAMPLES - 1);
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
			      regression2d_cb,
			      (void*)&us);
  
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
  if (rjmcmc_save_vector("callbackex.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  v = resultset2d_get_lambda(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get lambda data\n");
    return -1;
  }
  if (rjmcmc_save_vector("callbackex.lambda", v, total) < 0) {
    fprintf(stderr, "error: failed to save lambda data\n");
    return -1;
  }
  if (rjmcmc_save_vector_as_histogram("callbackex.lambda_hist", 
				      data->lambdamin,
				      data->lambdamax,
				      xsamples,
				      v, total) < 0) {
    fprintf(stderr, "error: failed to save lambda histogram\n");
    return -1;
  }

  iv = resultset2d_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("callbackex.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("callbackex.partition_hist",
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
  if (rjmcmc_save_matrix("callbackex.mean", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  m = resultset2d_get_median(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("callbackex.median", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  m = resultset2d_get_mode(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("callbackex.mode", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  m = resultset2d_get_credible_min(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("callbackex.conf_min", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  m = resultset2d_get_credible_max(results);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("callbackex.conf_max", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_max data\n");
    return -1;
  }

  im = resultset2d_get_centres(results);
  if (im == NULL) {
    fprintf(stderr, "error: failed to get centres\n");
    return -1;
  }
  if (rjmcmc_save_int_matrix("callbackex.centres", im, xsamples, ysamples) < 0) {
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

static void regression2d_cb(void *state,
			    double *x,
			    double *y,
			    int ncells,
			    regression2d_value_at_t value_at,
			    double lambda,
			    void *user_arg)
{
  static const int SAMPLERATE = 5;
  struct user_struct *u = (struct user_struct*)user_arg;
  int i;
  int j;

  char filename[256];
  FILE *fp;

  for (j = 0; j < NYSAMPLES; j ++) {
    for (i = 0; i < NXSAMPLES; i ++) {

      u->mean[i][j] = (u->mean[i][j] * u->n + 
		       value_at(state, u->xc[i], u->yc[j]))/
	(double)(u->n + 1);
    }
  }

  if (u->n % SAMPLERATE == 0 && u->n/SAMPLERATE < 2000) {
    sprintf(filename, "callbackex_mean_%08d.txt", u->n/SAMPLERATE + 1);
    fp = fopen(filename, "w");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to open file for writing\n");
      exit(-1);
    }
    
    
    for (j = 0; j < NYSAMPLES; j ++) {
      for (i = 0; i < NXSAMPLES; i ++) {
	fprintf(fp, "%g ", u->mean[i][j]);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
  }


  u->n ++;
}
