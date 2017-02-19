
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultset1dfm.h>

#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "ls:t:h";

static struct option long_options[] = {
  {"lambda", 0, 0, 'l'},
  {"sigmav", 1, 0, 's'},
  {"total", 1, 0, 't'},
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

#define DATASIZE 100
#define XSAMPLES 100

static const double xmin = 0.0;
static const double xmax = 100.0;
static const int xsamples = XSAMPLES;
static int npoints = -1;

static double fx(double x);

struct my_data {
  double x[DATASIZE];
  double y[DATASIZE];
  double n[DATASIZE];
};
  
static const double *test_value_at(part1d_fm_likelihood_state_t *state,
				   double x)
{
#define NPT 9
  double points[NPT][2] = {
    {13.43642,   347.43374},
    {76.37746,  -244.93097},
    {49.54351,   -50.50894},
    {65.15930,   288.72335},
    { 9.38596,  -471.65252},
    {83.57651,   -67.23293},
    {76.22801,  -497.89395},
    {44.53872,   221.54003},
    {22.87622,   445.27070}};
  int i;
  double dist;
  double mindist;
  int mindisti;
  static double result;

  mindist = 1e37;
  mindisti = -1;

  for (i = 0; i < NPT; i ++) {
    dist = fabs(x - points[i][0]);
    if (dist < mindist) {
      mindist = dist;
      mindisti = i;
    }
  }

  if (mindisti < 0) {
    result = 0.0;
  } else {
    result = points[mindisti][1];
  }
  
  return &result;
}
    
static double synthdata_misfit(void *user_arg,
			       int npartitions,
			       const double *partitions,
			       int nglobalparameters,
			       const double *global_parameters,
			       part1d_fm_likelihood_state_t *state,
			       part1d_fm_value_at_t value_at,
			       part1d_fm_value_at_t gradient_at);

static double synthdata_misfit_hierarchical(void *user_arg,
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
  int total = 200000;
  int thin = 1;
  int min_part = 2;
  int max_part = 40;
  int ysamples = 200;
  double confidence = 0.95;
  double pd = 2.5;
  double sigmav = 5.0;

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

  int c;
  int option_index;

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
	fprintf(stderr, "error: total must be greater than 20000\n");
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
   * Load our synthetic data
   */
  fp = fopen("data.txt", "r");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to open data file\n");
    return -1;
  }

  i = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf %lf %lf\n", &(data.x[i]), &(data.y[i]), &(data.n[i])) != 3) {
      if (feof(fp)) {
	continue;
      } else {
	fprintf(stderr, "error: failed to read data file\n");
	return -1;
      }
    }
    i ++;
  }
  npoints = i;

  fclose(fp);

  local_parameter.fmin = -500.0;
  local_parameter.fmax = 500.0;
  local_parameter.fstd_value = sigmav;
  local_parameter.fstd_bd = sigmav;

  /* printf("test: %f\n", *(test_value_at(NULL, 10.0))); */
  /* printf("test: %f\n", *(test_value_at(NULL, 15.0))); */
  /* printf("test: %f\n", *(test_value_at(NULL, 25.0))); */
  /* printf("test: %f\n", *(test_value_at(NULL, 70.707071))); */
  /* printf("npoints: %d\n", npoints); */
  /* printf("like: %f\n", synthdata_misfit(&data, */
  /* 					9, */
  /* 					NULL, */
  /* 					0, */
  /* 					NULL, */
  /* 					NULL, */
  /* 					test_value_at, */
  /* 					NULL)); */


  if (hierarchical) {
    
    hierarchical_parameter.fmin = 0.5;
    hierarchical_parameter.fmax = 2.0;
    hierarchical_parameter.fstd_value = 0.1;
    hierarchical_parameter.fstd_bd = 0.0;
    
    results = part1d_forwardmodel_hierarchical(burnin,
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
					       synthdata_misfit_hierarchical,
					       &data,
					       RESULTSET1DFM_MEAN);

  } else {
    
    results = part1d_forwardmodel(burnin,
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
				  synthdata_misfit,
				  &data,
				  RESULTSET1DFM_MEAN);
  }

  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run synthdata\n");
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
  for (i = 0; i < xsamples; i ++) {
    y[i] = fx(xcoords[i]);
  }

  v = resultset1dfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("synthdata.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultset1dfm_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector("synthdata.partitions", iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }
  if (rjmcmc_save_int_vector_as_histogram("synthdata.partition_hist",
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

  if (rjmcmc_save_int_coords("synthdata.partition_x_hist",
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
  if (rjmcmc_save_coords("synthdata.mean", xcoords, v, xsamples) < 0) {
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

    if (rjmcmc_save_vector_as_histogram("synthdata.sigma_histogram",
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

static double fx(double x)
{
  double y;

  if (x < 50.0) {
    y = 25.0;
  } else {
    y = 0.0;
  }

  return y;
}

static double synthdata_misfit(void *user_arg,
			  int npartitions,
			  const double *partitions,
			  int nglobalparameters,
			  const double *global_parameters,
			  part1d_fm_likelihood_state_t *state,
			  part1d_fm_value_at_t value_at,
			  part1d_fm_value_at_t gradient_at)
{
  double t;
  double fakesigma = 1.0;

  return synthdata_misfit_hierarchical(user_arg,
				  npartitions,
				  partitions,
				  nglobalparameters,
				  global_parameters,
				  0,
				  1,
				  &fakesigma,
				  (part1d_fm_hierarchical_likelihood_state_t *)state,
				  (part1d_fm_hierarchical_value_at_t)value_at,
				  (part1d_fm_hierarchical_value_at_t)gradient_at,
				  &t);
}

static double
synthdata_misfit_hierarchical(void *user_arg,
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
  double ds;
  double sum;

  FILE *fp;
  const double *local_parameters;

  sum = 0.0;

  for (i = 0; i < npoints; i ++) {

    local_parameters = value_at(state, data->x[i]);

    dv = data->y[i] - local_parameters[0];

    ds = hierarchicalparameters[0] * data->n[i];
    sum += (dv*dv)/(2.0 * ds*ds);

  }

  if (hierarchical) {
    *logdetce = 2.0 * xsamples * log(hierarchicalparameters[0]);
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
