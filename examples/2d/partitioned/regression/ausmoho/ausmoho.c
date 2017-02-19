
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/dataset2d.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "lt:h";

static struct option long_options[] = {
  {"lambda", 0, 0, 'l'},
  {"total", 1, 0, 't'},
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static double misfit_hierarchical(void *user,
				  int nglobalparameters,
				  const double *global_paramters,
				  int hierarchical,
				  int nhierarchicalparameters,
				  const double *hierarchical_parameters,
				  part2d_fm_likelihood_state_t *state,
				  part2d_fm_value_at_t value_at,
				  part2d_fm_value_at_t gradient_at,
				  const bbox2d_t *bound,
				  double *logdetce);

static double misfit(void *user,
                     int nglobalparameters,
                     const double *global_paramters,
                     part2d_fm_likelihood_state_t *state,
                     part2d_fm_value_at_t value_at,
                     part2d_fm_value_at_t gradient_at,
                     const bbox2d_t *bound);

static double global_lambda = 2.0;

int main(int argc, char *argv[]) 
{
  int c;
  int option_index;
  
  dataset2d_t *data;
  resultset2dfm_t *results;

  int burnin = 1000;
  int total = 400000;
  int thin = 1;
  int min_part = 2;
  int max_part = 500;
  int init_part = 50;
  int xsamples = 64;
  int ysamples = 64;
  int zsamples = 200;
  double credible = 0.95;
  double pv = 2.0;
  double pv_bd = 2.2;
  double pd = 0.5;
  int lambda = 0;

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

  forwardmodelparameter_t local_parameters;
  forwardmodelparameter_t sigma_parameters;

  data = dataset2d_load_known("data.txt");
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
      lambda = -1;
      break;

    case 't':
      total = atoi(optarg);
      if (total < 1000) {
	fprintf(stderr, "error: total must be greater than 1000\n");
	return -1;
      }
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
  printf("%d points\n", data->npoints);
  printf("Auto xrange: %f %f\n", data->xmin, data->xmax);
  printf("Auto yrange: %f %f\n", data->ymin, data->ymax);
  printf("Auto zrange: %f %f\n", data->zmin, data->zmax);
  data->xmin = 100.0;
  data->xmax = 160.0;
  data->ymin = -55.0;
  data->ymax = 0.0;
  data->zmin = 10.0;
  data->zmax = 62.0;

  local_parameters.fmin = data->zmin;
  local_parameters.fmax = data->zmax;
  local_parameters.fstd_value = pv;
  local_parameters.fstd_bd = pv_bd;

  sigma_parameters.fmin = 0.05;
  sigma_parameters.fmax = 4.5;
  sigma_parameters.fstd_value = 0.02;
  sigma_parameters.fstd_bd = 0.0;
    
  position_map2d_set_type(0);

  if (lambda) {
    results = part2d_forwardmodel_hierarchical(burnin,
					       total,
					       thin,
					       min_part,
					       max_part,
					       init_part,
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
					       &local_parameters,
					       1,
					       &sigma_parameters,
					       misfit_hierarchical,
					       (void*)data,
					       RESULTSET2DFM_MEAN |
					       RESULTSET2DFM_MEDIAN |
					       RESULTSET2DFM_MODE |
					       RESULTSET2DFM_CREDIBLE);
  } else {
    results = part2d_forwardmodel(burnin,
				  total,
				  thin,
				  min_part,
				  max_part,
				  init_part,
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
				  &local_parameters,
				  misfit,
				  (void*)data,
				  RESULTSET2DFM_MEAN |
				  RESULTSET2DFM_MEDIAN |
				  RESULTSET2DFM_MODE |
				  RESULTSET2DFM_CREDIBLE);
  }

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
  if (rjmcmc_save_vector("ausmoho.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }
  
  iv = resultset2dfm_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("ausmoho.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("ausmoho.partition_hist",
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
  if (rjmcmc_save_matrix("ausmoho.mean", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_median(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get median data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("ausmoho.median", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save median data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_mode(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get mode data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("ausmoho.mode", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save mode data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_credible_min(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_min data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("ausmoho.credible_min", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_min data\n");
    return -1;
  }

  m = resultset2dfm_get_local_parameter_credible_max(results, 0);
  if (m == NULL) {
    fprintf(stderr, "error: failed to get credible_max data\n");
    return -1;
  }
  if (rjmcmc_save_matrix("ausmoho.credible_max", m, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save credible_max data\n");
    return -1;
  }

  im = resultset2dfm_get_centres(results);
  if (im == NULL) {
    fprintf(stderr, "error: failed to get centres\n");
    return -1;
  }
  if (rjmcmc_save_int_matrix("ausmoho.centres", im, xsamples, ysamples) < 0) {
    fprintf(stderr, "error: failed to save centres\n");
    return -1;
  }
  
  if (lambda) {
    v = resultset2dfm_get_hierarchical_parameter(results, 0);
    if (v == NULL) {
      fprintf(stderr, "error: failed to get hierarchical parameter history\n");
      return -1;
    }
    
    if (rjmcmc_save_vector_as_histogram("ausmoho.sigma_histogram",
					sigma_parameters.fmin,
					sigma_parameters.fmax,
					100,
					v,
					total) < 0) {
      fprintf(stderr, "error: failed to save hierarchical parameter to histogram\n");
      return -1;
    }
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
	  " -h|--help               show usage information\n"
	  "\n",
	  pname);
}

static double misfit_hierarchical(void *user,
				  int nglobalparameters,
				  const double *global_paramters,
				  int hierarchical,
				  int nhierarchicalparameters,
				  const double *hierarchical_parameters,
				  part2d_fm_likelihood_state_t *state,
				  part2d_fm_value_at_t value_at,
				  part2d_fm_value_at_t gradient_at,
				  const bbox2d_t *bound,
				  double *logdetce)
{
  dataset2d_t *data = (dataset2d_t *)user;
  int i;

  double sum;
  double sigma2;
  double dz;
  double n;
  const double *lp;
  double sigma;

  sum = 0.0;
  sigma = hierarchical_parameters[0];

  for (i = 0; i < data->npoints; i ++) {
    
    lp = value_at(state, data->points[i].x, data->points[i].y);
    if (lp == NULL) {
      fprintf(stderr, "misfit: failed to determine local value\n");
      exit(-1);
    }

    dz = data->points[i].z - lp[0];
    //n = sigma * data->points[i].n;
    n = sigma;

    if (n <= 0.0) {
      fprintf(stderr, "misfit: invalid noise at %d: %f\n", i, n);
      exit(-1);
    }
    sum += (dz*dz)/(2.0 * (n*n));
  }

  if (hierarchical) {
    *logdetce = 2.0 * (double)data->npoints * log(sigma);
  }

  return sum;
}

static double misfit(void *user,
                     int nglobalparameters,
                     const double *global_parameters,
                     part2d_fm_likelihood_state_t *state,
                     part2d_fm_value_at_t value_at,
                     part2d_fm_value_at_t gradient_at,
                     const bbox2d_t *bound)
{
  double dummy_hierarchical;
  double dummy_logdetce;

  dummy_hierarchical = global_lambda;
  
  return misfit_hierarchical(user,
			     nglobalparameters,
			     global_parameters,
			     0,
			     1,
			     &dummy_hierarchical,
			     state,
			     value_at,
			     gradient_at,
			     bound,
			     &dummy_logdetce);
}
