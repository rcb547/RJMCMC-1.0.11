
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <mpi.h>

#include <rjmcmc/dataset2d.h>

#include <rjmcmc/forwardmodel_mpi.h>

#if !defined(HAVE_MPI_H)
#error "No MPI."
#endif

#include <rjmcmc/position_map2d.h>

#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

static char short_options[] = "hm:p:";

static struct option long_options[] = {
  {"help", 0, 0, 'h'},
  {"method", 1, 0, 'm'},
  {"partitions", 1, 0, 'p'},
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
  int min_part = 2;
  int max_part = 100;
  int xsamples = 100;
  int ysamples = 100;
  int zsamples = 200;
  double confidence = 0.95;
  double pd = 5.0;
  int method = 0;
  int seed_base = 101;
  int seed_scale = 983;

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

  double pv = 5.0;
  double sigma = 10.0;

  forwardmodelparameter_t local_parameter;

  int mpi_size;
  int mpi_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  /* printf("%d/%d\n", mpi_rank, mpi_size); */

  rjmcmc_seed(seed_base + mpi_rank * seed_scale);

  data = dataset2d_load_fixed("data.txt", sigma);
  if (data == NULL) {
    fprintf(stderr, 
	    "error: unable to load data, "
	    "has it been created with the script?\n");
    return -1;
  }

  if (mpi_rank == 0) {
    printf("Data range: %f %f\n", data->zmin, data->zmax);
  }

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'm':
      method = atoi(optarg);
      break;

    case 'p':
      max_part = atoi(optarg);
      if (max_part < 50) {
	fprintf(stderr, "error: no. partitions must be greater than 50\n");
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

  if (method < 0 || method > 2) {
    fprintf(stderr, "error: invalid method\n");
    return -1;
  }
  
  position_map2d_set_type(method);

  /*
   * Set the data bounds
   */
  if (mpi_rank == 0) {
    printf("Auto zrange: %f %f\n", data->zmin, data->zmax);
  }
  data->xmin = -50.0;
  data->xmax = 50.0;
  data->ymin = -50.0;
  data->ymax = 50.0;
  data->zmin = -50.0;
  data->zmax = 50.0;

  local_parameter.fmin = data->zmin;
  local_parameter.fmax = data->zmax;
  local_parameter.fstd_value = pv;
  local_parameter.fstd_bd = pv;

  results = MPI_part2d_forwardmodel(burnin,
				    total,
				    0, /* Thin */
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
				    misfit,
				    (void*)data,
				    RESULTSET2DFM_MEAN |
				    RESULTSET2DFM_MEDIAN |
				    RESULTSET2DFM_MODE |
				    RESULTSET2DFM_CREDIBLE,
				    mpi_size,
				    mpi_rank,
				    0, /* Set the root to be 0 */
				    MPI_COMM_WORLD);
  
  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run regression\n");
    return -1;
  }

  if (mpi_rank == 0) {
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
    if (rjmcmc_save_vector("regression.misfit", v, total) < 0) {
      fprintf(stderr, "error: failed to save misfit data\n");
      return -1;
    }

    iv = resultset2dfm_get_partitions(results);
    if (iv == NULL) {
      fprintf(stderr, "error: failed to get partitions data\n");
      return -1;
    }
    if (rjmcmc_save_int_vector("regression.partitions", iv, total) < 0) {
      fprintf(stderr, "error: failed to save partitions data\n");
      return -1;
    }
    if (rjmcmc_save_int_vector_as_histogram("regression.partition_hist",
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
    if (rjmcmc_save_matrix("regression.mean", m, xsamples, ysamples) < 0) {
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
    if (rjmcmc_save_matrix("regression.mode", m, xsamples, ysamples) < 0) {
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
    if (rjmcmc_save_matrix("regression.median", m, xsamples, ysamples) < 0) {
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
    if (rjmcmc_save_matrix("regression.credible_min", m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save credible_min data\n");
      return -1;
    }
    
    m = resultset2dfm_get_local_parameter_credible_max(results, 0);
    if (m == NULL) {
      fprintf(stderr, "error: failed to get credible_max data\n");
      return -1;
    }
    if (rjmcmc_save_matrix("regression.credible_max", m, xsamples, ysamples) < 0) {
      fprintf(stderr, "error: failed to save credible_max data\n");
      return -1;
    }
    
    
    im = resultset2dfm_get_centres(results);
    if (im == NULL) {
      fprintf(stderr, "error: failed to get centres\n");
      return -1;
    }
    if (rjmcmc_save_int_matrix("regression.centres", im, xsamples, ysamples) < 0) {
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

    rjmcmc_destroy_array_1d(xcoords);
    rjmcmc_destroy_array_1d(ycoords);

  }
    
  dataset2d_destroy(data);
  resultset2dfm_destroy(results);

  MPI_Finalize();

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr, 
	  "usage: %s [options]\n"
	  "where options is on or more of:\n"
	  "\n"
	  " -p|--partitions <int>   Max paritions\n"
	  " -m|--method <int>       2D Nearest neighbour search method:\n"
	  "                          0 - Linear (default)\n"
	  "                          1 - Delaunay Triangulation\n"
	  "                          2 - Quadtree\n"
	  "\n"
	  " -h|--help               show usage information\n"
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

    dz = data->points[i].z - lp[0];
    n = data->points[i].n;
    sum += (dz*dz)/(2.0 * (n*n));
  }

  return sum;
}
