
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <rjmcmc/dataset2d.h>
#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "b:t:T:lh";

static struct option long_options[] = {
  {"burnin", 1, 0, 'b'},
  {"total", 1, 0, 't'},
  {"thin", 1, 0, 'T'},
  {"lambda", 0, 0, 'l'},
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static double misfit(void *user,
		     int nglobalparameters,
		     const double *global_paramters,
		     part2d_fm_likelihood_state_t *state,
		     part2d_fm_value_at_t value_at,
		     part2d_fm_value_at_t gradient_at,
		     const bbox2d_t *bound);

int main(int argc, char *argv[]) 
{
  int c;
  int option_index;
  
  dataset2d_t *data;
  resultset2dfm_t *results;

  int burnin = 10000;
  int total = 50000;
  int thin = 1;
  int min_part = 2;
  int max_part = 100;
  int xsamples = 100;
  int ysamples = 100;
  int zsamples = 200;
  double credible = 0.95;
  double pv = 0.5;
  double pv_bd = 1.0;

  double pd = 0.25;

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

  double sigma = 15.0;

  forwardmodelparameter_t local_parameter;

  data = dataset2d_load_fixed("data.txt", sigma);
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

    case 'b':
      burnin = atoi(optarg);
      break;

    case 't':
      total = atoi(optarg);
      break;

    case 'T':
      thin = atoi(optarg);
      break;

    default:
      fprintf(stderr, "error: invalid option\n");
      return -1;

    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (burnin < 1000) {
    fprintf(stderr, "error: burnin must be greater than 1000\n");
    return -1;
  }

  if (total < 10000 ||
      total < burnin) {
    fprintf(stderr, 
	    "error: total must be greater than 10000 and greater than "
	    "the burnin parameter\n");
    return -1;
  }

  if (thin < 0) {
    fprintf(stderr, 
	    "error: thin must be greater than 0.\n");
    return -1;
  }

  /*
   * Set the data bounds
   */
  printf("%d points\n", data->npoints);
  printf("Auto xrange: %f %f\n", data->xmin, data->xmax);
  printf("Auto yrange: %f %f\n", data->ymin, data->ymax);
  printf("Auto zrange: %f %f\n", data->zmin, data->zmax);
  data->xmin = -50.0;
  data->xmax = 50.0;
  data->ymin = -50.0;
  data->ymax = 50.0;

  /* data->zmin = -50.0; */
  /* data->zmax = 50.0; */
  data->zmin = -5.0;
  //data->zmax = 30.0;
  data->zmax = 60.0;

  local_parameter.fmin = data->zmin;
  local_parameter.fmax = data->zmax;
  local_parameter.fstd_value = pv;
  local_parameter.fstd_bd = pv_bd;

  position_map2d_set_type(0);

  results = part2d_forwardmodel(burnin,
				total,
				thin,
				min_part,
				max_part,
				data->xmin,
				data->xmax,
				data->ymin,
				data->ymax,
				xsamples,
				ysamples,
				zsamples,
				credible,
				pd,
				pd,
				rjmcmc_uniform,
				rjmcmc_normal,
				0,
				NULL,
				1,
				&local_parameter,
				misfit,
				(void*)data,
				RESULTSET2DFM_MEAN |
				RESULTSET2DFM_MEDIAN |
				RESULTSET2DFM_MODE |
				RESULTSET2DFM_CREDIBLE);
  
  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run regression\n");
    return -1;
  }
  xcl = xsamples;
  xcoords = rjmcmc_create_array_1d(xcl);
  resultset2dfm_fill_xcoord_vector(results, xcoords, &xcl);
  
  ycl = ysamples;
  ycoords = rjmcmc_create_array_1d(ycl);
  resultset2dfm_fill_ycoord_vector(results, ycoords, &ycl);

  v = resultset2dfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("gaussian.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }
  

  iv = resultset2dfm_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("gaussian.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("gaussian.partition_hist",
					  0,
					  max_part,
					  iv,
					  total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_mean(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("gaussian.mean", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_median(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("gaussian.median", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_mode(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("gaussian.mode", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_credible_min(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("gaussian.credible_min", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_credible_max(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("gaussian.credible_max", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_max data\n");
    return -1;
  }

  im = resultset2dfm_get_centres(results);
  if (im == NULL) {
    fprintf(stderr, "error: failed to get centres\n");
    return -1;
  }
  if (rjmcmc_save_int_matrix("gaussian.centres", im, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save centres\n");
    return -1;
  }
  

  iv = resultset2dfm_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultset2dfm_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  dataset2d_destroy(data);
  resultset2dfm_destroy(results);

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
	  " -b|--burnin <int>       The number of burnin steps\n"
	  " -t|--total <int>        The total number of steps\n"
	  "\n"
	  " -h|--help               Show usage information and exit\n"
	  "\n",
	  pname);
}

static double misfit(void *user,
		     int nglobalparameters,
		     const double *global_paramters,
		     part2d_fm_likelihood_state_t *state,
		     part2d_fm_value_at_t value_at,
		     part2d_fm_value_at_t gradient_at,
		     const bbox2d_t *bound)
{
  dataset2d_t *data = (dataset2d_t *)user;
  int i;

  double sum;
  double sigma2;
  double dz;
  double n;
  const double *lp;

  sum = 0.0;

  for (i = 0; i < data->npoints; i ++) {
    
    lp = value_at(state, data->points[i].x, data->points[i].y);
    if (lp == NULL) {
      fprintf(stderr, "misfit: failed to determine local value\n");
      exit(-1);
    }

    /* if (lp[0] > 20.0) { */
    /*   printf("lp: %f\n", lp[0]); */
    /* } */
    dz = data->points[i].z - lp[0];

    /* if (data->points[i].z > 20.0) { */
    /*   printf("dz: %02d %f\n", i, data->points[i].z); */
    /* } */
    n = data->points[i].n;
    sum += (dz*dz)/(2.0 * (n*n));
  }

  //  printf("%f %f\n", sum, n);
  return sum;
}
