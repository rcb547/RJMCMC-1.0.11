
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultset1dfm.h>

#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "t:ls:h";

static struct option long_options[] = {
  {"lambda", 0, 0, 'l'},
  {"total", 1, 0, 't'},
  {"sigmav", 1, 0, 's'},
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

#define XSAMPLES 100

static double xmin = 0.0;
static double xmax = 100.0;
static const int xsamples = XSAMPLES;
static int npoints;
static const double sigma = 1.0;

static double fx(double x);

struct my_data {
  double x[XSAMPLES];
  double y[XSAMPLES];
  double n[XSAMPLES];
};
  
static double natural_misfit(void *user_arg,
			     int npartitions,
			     const double *partitions,
			     int nglobalparameters,
			     const double *global_parameters,
			     part1d_fm_likelihood_state_t *state,
			     part1d_fm_value_at_t value_at,
			     part1d_fm_value_at_t gradient_at);

static double natural_misfit_hierarchical(void *user_arg,
					  int npartitions,
					  const double *partitions,
					  int nglobalparameters,
					  const double *global_parameters,
					  int hierarchical,
					  int nhierarchical,
					  const double *hierarchicalparameters,
					  part1d_fm_hierarchical_likelihood_state_t *state,
					  part1d_fm_hierarchical_value_at_t value_at,
					  part1d_fm_hierarchical_value_at_t gradient_at,
					  double *logdetce);

static void usage(const char *pname);

int main(int argc, char *argv[]) 
{
  int burnin = 10000;
  int total = 100000;
  int thin = 1;
  int min_part = 2;
  int max_part = 60;
  int ysamples = 200;
  double confidence = 0.95;
  double pd = 0.8;
  FILE *fp;

  forwardmodelparameter_t local_parameter;

  int nproc;
  const double *v;
  const int *iv;
  int i;

  double *xcoords;
  double *y;
  int xcl;

  resultset1dfm_t *results;

  struct my_data data;

  int hierarchical = 0;

  forwardmodelparameter_t hierarchical_parameter;

  double sigmav = 2.0;

  int c;
  int option_index;

  double tsigma;
  double x;

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'l':
      hierarchical = 1;
      break;

    case 't':
      total = atoi(optarg);
      if (total < 20000) {
	fprintf(stderr, "error: invalid total\n");
	return -1;
      }
      break;

    case 's':
      sigmav = atof(optarg);
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
   * Load our data
   */

  fp = fopen("data.txt", "r");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to open data file\n");
    return -1;
  }

  i = 0;
  local_parameter.fmin = 1e37;
  local_parameter.fmax = -1e37;
    
  while (!feof(fp)) {
    if (fscanf(fp, "%lf %lf %lf\n", &(data.x[i]), &(data.y[i]), &(data.n[i])) != 3) {
      if (feof(fp)) {
	continue;
      } else {
	fprintf(stderr, "error: failed to read data point\n");
	return -1;
      }
    }

    if (data.y[i] < local_parameter.fmin) {
      local_parameter.fmin = data.y[i];
    }

    if (data.y[i] > local_parameter.fmax) {
      local_parameter.fmax = data.y[i];
    }

    i ++;
  }
  fclose(fp);

  npoints = i;

  xmin = data.x[0];
  xmax = data.x[i - 1];

  printf("%d points\n", npoints);
  printf("x %f %f\n", xmin, xmax);
  printf("y %f %f %f\n", local_parameter.fmin, local_parameter.fmax, local_parameter.fmax - local_parameter.fmin);

  xmin = 0.0;
  xmax = 100.0;

  local_parameter.fstd_value = sigmav;
  local_parameter.fstd_bd = sigmav;

  if (hierarchical) {
    
    hierarchical_parameter.fmin = 0.1;
    hierarchical_parameter.fmax = 1.5;
    hierarchical_parameter.fstd_value = 0.1;
    hierarchical_parameter.fstd_bd = 0.0;
    
    results = part1d_forwardmodel_natural_hierarchical(burnin,
						       total,
						       min_part,
						       max_part,
						       xmin,
						       xmax,
						       xsamples,
						       ysamples,
						       confidence,
						       pd,
						       rjmcmc_uniform,
						       rjmcmc_normal,
						       0,
						       NULL,
						       1,
						       &local_parameter,
						       1,
						       &hierarchical_parameter,
						       natural_misfit_hierarchical,
						       &data,
						       RESULTSET1DFM_MEAN);

  } else {
    
    results = part1d_forwardmodel_natural(burnin,
					  total,
					  min_part,
					  max_part,
					  xmin,
					  xmax,
					  xsamples,
					  ysamples,
					  confidence,
					  pd,
					  rjmcmc_uniform,
					  rjmcmc_normal,
					  0,
					  NULL,
					  1,
					  &local_parameter,
					  natural_misfit,
					  &data,
					  RESULTSET1DFM_MEAN);
  }
  
  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run natural\n");
    return -1;
  }

  xcl = xsamples;
  xcoords = rjmcmc_create_array_1d(xcl);
  if (xcoords == NULL) {
    fprintf(stderr, "error: failed to create array for xsamples\n");
    return -1;
  }
  y = rjmcmc_create_array_1d(xcl);
  if (y == NULL) {
    fprintf(stderr, "error: failed to create array for y data\n");
    return -1;
  }
  
  resultset1dfm_fill_xcoord_vector(results, xcoords, &xcl);

  v = resultset1dfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("natural.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultset1dfm_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("natural.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("natural.partition_hist",
					  2,
					  max_part,
					  iv, total) < 0) {
    fprintf(stderr, "error: failed to save partition histogram data\n");
    return -1;
  }

  iv = resultset1dfm_get_partition_x_histogram(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partition histogram\n");
    return -1;
  }

  if (rjmcmc_save_int_coords("natural.partition_x_hist",
			     xcoords,
			     iv,
			     xsamples) < 0) {
    fprintf(stderr, "error: failed to save partition x histogram\n");
    return -1;
  }

  v = resultset1dfm_get_local_parameter_mean(results, 0);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get mean data\n");
    return -1;
  }
  if (rjmcmc_save_coords("natural.mean", xcoords, v, xsamples) < 0) {
    fprintf(stderr, "error: failed to save mean data\n");
    return -1;
  }

  iv = resultset1dfm_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultset1dfm_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  if (hierarchical) {
    v = resultset1dfm_get_hierarchical(results, 0);
    
    if (v == NULL) {
      fprintf(stderr, "error: failed to get hierarchical parameter history\n");
      return -1;
    }

    if (rjmcmc_save_vector_as_histogram("natural.sigma_histogram",
					hierarchical_parameter.fmin,
					hierarchical_parameter.fmax,
					100,
					v,
					total) < 0) {
      fprintf(stderr, "error: failed to save hierarchical parameter to histogram\n");
      return -1;
    }

  }

  rjmcmc_destroy_array_1d(xcoords);
  resultset1dfm_destroy(results);

  return 0;
}

static double 
natural_misfit(void *user_arg,
	       int npartitions,
	       const double *partitions,
	       int nglobalparameters,
	       const double *global_parameters,
	       part1d_fm_likelihood_state_t *state,
	       part1d_fm_value_at_t value_at,
	       part1d_fm_value_at_t gradient_at)
{
  double t;
  
  return natural_misfit_hierarchical(user_arg,
				     npartitions,
				     partitions,
				     nglobalparameters,
				     global_parameters,
				     0,
				     1,
				     &sigma,
				     (part1d_fm_hierarchical_likelihood_state_t *)state,
				     (part1d_fm_hierarchical_value_at_t)value_at,
				     (part1d_fm_hierarchical_value_at_t)gradient_at,
				     &t);
}

static double
natural_misfit_hierarchical(void *user_arg,
			    int npartitions,
			    const double *partitions,
			    int nglobalparameters,
			    const double *global_parameters,
			    int hierarchical,
			    int nhierarchical,
			    const double *hierarchicalparameters,
			    part1d_fm_hierarchical_likelihood_state_t *state,
			    part1d_fm_hierarchical_value_at_t value_at,
			    part1d_fm_hierarchical_value_at_t gradient_at,
			    double *logdetce)
  
{
  struct my_data *data = (struct my_data*)user_arg;
  int i;
  double dv;
  double sum;
  double ss;
  double sssum;

  const double *local_parameters;

  sum = 0.0;
  sssum = 0.0;

  for (i = 0; i < npoints; i ++) {

    local_parameters = value_at(state, data->x[i]);

    dv = data->y[i] - local_parameters[0];

    ss = data->n[i] * hierarchicalparameters[0];
    ss = ss*ss;

    sssum = sssum + ss;

    sum += dv*dv/(2.0 * ss);
  }

  if (hierarchical) {
    *logdetce = log(sssum);
  }

  return sum;
}

static void usage(const char *pname)
{
  fprintf(stderr, 
	  "usage: %s [options]\n"
	  "where options is on or more of:\n"
	  "\n"
	  " -l|--lambda             use hierarchical\n"
	  "\n"
	  " -h|--help               show usage information\n"
	  "\n",
	  pname);
}
