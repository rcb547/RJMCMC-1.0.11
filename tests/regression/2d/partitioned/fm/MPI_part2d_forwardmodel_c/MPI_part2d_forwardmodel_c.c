
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <rjmcmc/dataset2d.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/forwardmodel_mpi.h>

#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "h";

static struct option long_options[] = {
  {"help", 0, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static double fxy(double x, double y);

static double likelihood(void *user,
			 int nglobalparameters,
			 const double *global_paramters,
			 part2d_fm_likelihood_state_t *state,
			 part2d_fm_value_at_t value_at,
			 part2d_fm_value_at_t gradient_at);

int main(int argc, char *argv[]) 
{
  int c;
  int option_index;
  
  dataset2d_t *data;
  resultset2dfm_t *results;

  int burnin = 10000;
  int total = 50000;
  int min_part = 2;
  int max_part = 100;
  int xsamples = 100;
  int ysamples = 100;
  int zsamples = 200;
  double confidence = 0.95;
  double pd = 0.5;
  int npoints = 2500;

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

  double pv = 1.0;
  double sigma = 5.0;
  double xmin = -5.0;
  double xmax = 5.0;
  double ymin = -5.0;
  double ymax = 5.0;

  forwardmodelparameter_t local_parameter;

  int mpi_size;
  int mpi_rank;

  char filename[256];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  /*
   * Initialize the rng
   */
  rjmcmc_seed(mpi_rank * 913);

  data = dataset2d_allocate(npoints);
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to create data\n");
    return -1;
  }

  data->xmin = xmin;
  data->xmax = xmax;
  data->ymin = ymin;
  data->ymax = ymax;
  data->zmin = -50.0;
  data->zmax = 50.0;

  if (mpi_rank == 0) {

    for (i = 0; i < npoints; i ++) {
      data->points[i].x = rjmcmc_random_choose_double(xmin, xmax, rjmcmc_uniform);
      data->points[i].y = rjmcmc_random_choose_double(ymin, ymax, rjmcmc_uniform);
      
      data->points[i].z = 
	fxy(data->points[i].x, data->points[i].y) + 
	sigma * rjmcmc_normal();
      
      data->points[i].n = sigma;
    }
  }

  MPI_Bcast(data->points, npoints * 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  printf("Data range: %f %f\n", data->zmin, data->zmax);
  printf("%d %f %f %f %f\n", mpi_rank, 
	 data->points[0].x,
	 data->points[0].y,
	 data->points[0].z,
	 data->points[0].n);

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
  local_parameter.fmin = data->zmin;
  local_parameter.fmax = data->zmax;
  local_parameter.fstd_value = pv;
  local_parameter.fstd_bd = pv;

  results = MPI_part2d_forwardmodel(burnin,
				    total,
				    min_part,
				    max_part,
				    data->xmin,
				    data->xmax,
				    data->ymin,
				    data->ymax,
				    xsamples,
				    ysamples,
				    zsamples,
				    confidence,
				    pd,
				    pd,
				    rjmcmc_uniform,
				    rjmcmc_normal,
				    0,
				    NULL,
				    1,
				    &local_parameter,
				    likelihood,
				    (void*)data,
				    RESULTSET2DFM_MEAN |
				    RESULTSET2DFM_MEDIAN |
				    RESULTSET2DFM_MODE |
				    RESULTSET2DFM_CREDIBLE,
				    mpi_size,
				    mpi_rank,
				    0,
				    MPI_COMM_WORLD);
  
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
  sprintf(filename, "MPI_part2d_forwardmodel_c.misfit_%d", mpi_rank);
  if (rjmcmc_save_vector(filename, v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultset2dfm_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  sprintf(filename, "MPI_part2d_forwardmodel_c.partitions_%d", mpi_rank);
  if (rjmcmc_save_int_vector(filename, iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }

  if (mpi_rank == 0) {
    if (rjmcmc_save_int_vector_as_histogram("MPI_part2d_forwardmodel_c.partition_hist",
					    0,
					    max_part,
					    iv,
					    total) < 0) {
      fprintf(stderr, "error: failed to save partitions data\n");
      return -1;
    }
    
    /*
     * Mean
     */
    m = resultset2dfm_get_local_parameter_mean(results, 0);
    if (m == NULL) {
      fprintf(stderr, "error: failed to get mean data\n");
      return -1;
    }
    if (rjmcmc_save_matrix("MPI_part2d_forwardmodel_c.mean", m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save mean data\n");
      return -1;
    }
    
    /*
     * Mode 
     */
    m = resultset2dfm_get_local_parameter_mode(results, 0);
    if (m == NULL) {
      fprintf(stderr, "error: failed to get mode data\n");
      return -1;
    }
    if (rjmcmc_save_matrix("MPI_part2d_forwardmodel_c.mode", m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save mode data\n");
      return -1;
    }
    
    /*
     * Median
     */
    m = resultset2dfm_get_local_parameter_median(results, 0);
    if (m == NULL) {
      fprintf(stderr, "error: failed to get median data\n");
      return -1;
    }
    if (rjmcmc_save_matrix("MPI_part2d_forwardmodel_c.median", m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save median data\n");
      return -1;
    }
    
    /*
     * Credible intervals
     */
    m = resultset2dfm_get_local_parameter_credible_min(results, 0);
    if (m == NULL) {
      fprintf(stderr, "error: failed to get credible_min data\n");
      return -1;
    }
    if (rjmcmc_save_matrix("MPI_part2d_forwardmodel_c.credible_min", m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save credible_min data\n");
      return -1;
    }

    m = resultset2dfm_get_local_parameter_credible_max(results, 0);
    if (m == NULL) {
      fprintf(stderr, "error: failed to get credible_max data\n");
      return -1;
    }
    if (rjmcmc_save_matrix("MPI_part2d_forwardmodel_c.credible_max", m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save credible_max data\n");
      return -1;
    }
    
    
    im = resultset2dfm_get_centres(results);
    if (im == NULL) {
      fprintf(stderr, "error: failed to get centres\n");
      return -1;
    }
    if (rjmcmc_save_int_matrix("MPI_part2d_forwardmodel_c.centres", im, xsamples, ysamples) < 0) {
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
  }

  if (mpi_rank != 0) {
    /*
     * Mean
     */
    m = resultset2dfm_get_local_parameter_mean(results, 0);
    if (m == NULL) {
      fprintf(stderr, "error: failed to get mean data\n");
      return -1;
    }
    sprintf(filename, "MPI_part2d_forwardmodel_c.mean_%d", mpi_rank);
    if (rjmcmc_save_matrix(filename, m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save mean data\n");
      return -1;
    }
  }

  dataset2d_destroy(data);
  resultset2dfm_destroy(results);
  
  rjmcmc_destroy_array_1d(xcoords);
  rjmcmc_destroy_array_1d(ycoords);

  MPI_Finalize();

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

static double fxy(double x, double y)
{
  // Function is a raised disc area at 0,0
  double r = sqrt(x*x + y*y);
  if (r < 2.5) {
    return 25.0;
  } else {
    return 0.0;
  }
}

static double likelihood(void *user,
			 int nglobalparameters,
			 const double *global_paramters,
			 part2d_fm_likelihood_state_t *state,
			 part2d_fm_value_at_t value_at,
			 part2d_fm_value_at_t gradient_at)
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
      fprintf(stderr, "likelihood: failed to determine local value\n");
      exit(-1);
    }

    dz = data->points[i].z - lp[0];
    n = data->points[i].n;
    sum += (dz*dz)/(2.0 * (n*n));
  }

  return sum;
}
