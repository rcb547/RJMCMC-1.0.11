
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <mpi.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/forwardmodel_mpi.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultset1dfm.h>

#include <rjmcmc/rjmcmc_util.h>

#define XSAMPLES 100

static const double xmin = 0.0;
static const double xmax = 100.0;
static const int xsamples = XSAMPLES;
static const double base_sigma = 5.0;

static double fx(double x);

struct my_data {
  double x[XSAMPLES];
  double ry[XSAMPLES];
  double sy[XSAMPLES];
};
  
static double likelihood(void *user_arg,
			 int npartitions,
			 const double *partitions,
			 int nglobalparameters,
			 const double *global_parameters,
			 int hierarchical,
			 int nhierarchicalparameters,
			 const double *hierarchical_parameters,
			 part1d_fm_hierarchical_likelihood_state_t *state,
			 part1d_fm_hierarchical_value_at_t value_at,
			 part1d_fm_hierarchical_value_at_t gradient_at,
			 double *logdetce);

int main(int argc, char *argv[]) 
{
  int burnin = 1000;
  int total = 20000;
  int min_part = 2;
  int max_part = 10;
  int ysamples = 200;
  double confidence = 0.95;
  double pd = 3.0;

  forwardmodelparameter_t local_parameter;
  forwardmodelparameter_t hierarchical_parameter;

  int nproc;
  const double *v;
  const int *iv;
  int i;

  double *xcoords;
  double *y;
  int xcl;

  resultset1dfm_t *results;

  struct my_data data;

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

  /*
   * Create our data based upon sampling our function
   */
  if (mpi_rank == 0) {

    for (i = 0; i < xsamples; i ++) {
      data.x[i] = (double)i/(double)(xsamples - 1) * (xmax - xmin) + xmin;
      data.ry[i] = fx(data.x[i]);
      data.sy[i] = data.ry[i] + rjmcmc_normal() * base_sigma;
    }

  }
    
  /*
   * Make sure all processes have the same data.
   */
  MPI_Bcast(data.x, xsamples, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(data.ry, xsamples, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(data.sy, xsamples, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  local_parameter.fmin = -150.0;
  local_parameter.fmax = 150.0;
  local_parameter.fstd_value = 10.0;
  local_parameter.fstd_bd = 10.0;

  hierarchical_parameter.fmin = base_sigma * 0.5;
  hierarchical_parameter.fmax = base_sigma * 2.0;
  hierarchical_parameter.fstd_value = 0.1;
  hierarchical_parameter.fstd_bd = 0.0;


  results = MPI_part1d_forwardmodel_natural_hierarchical(burnin,
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
							 likelihood,
							 &data,
							 RESULTSET1DFM_MEAN,
							 mpi_size,
							 mpi_rank,
							 0,
							 MPI_COMM_WORLD);
  
  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run part1d_forwardmodel_natural_hierarchical_c\n");
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
  if (mpi_rank == 0) {
    for (i = 0; i < xsamples; i ++) {
      y[i] = fx(xcoords[i]);
    }
    if (rjmcmc_save_vector("part1d_forwardmodel_natural_hierarchical_c.data", data.sy, xsamples) < 0) {
      fprintf(stderr, "error: failed to save data\n");
      return -1;
    }
  }

  v = resultset1dfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  sprintf(filename, "part1d_forwardmodel_natural_hierarchical_c.misfit_%d", mpi_rank);
  if (rjmcmc_save_vector(filename, v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultset1dfm_get_partitions(results);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get partitions data\n");
    return -1;
  }
  sprintf(filename, "part1d_forwardmodel_natural_hierarchical_c.partitions_%d", mpi_rank);
  if (rjmcmc_save_int_vector(filename, iv, total) < 0) {
    fprintf(stderr, "error: failed to save partitions data\n");
    return -1;
  }

  if (mpi_rank == 0) {
    
    if (rjmcmc_save_int_vector_as_histogram("part1d_forwardmodel_natural_hierarchical_c.partition_hist",
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
    
    if (rjmcmc_save_int_coords("part1d_forwardmodel_natural_hierarchical_c.partition_x_hist",
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
    if (rjmcmc_save_vector("part1d_forwardmodel_natural_hierarchical_c.mean", v, xsamples) < 0) {
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
  }

  rjmcmc_destroy_array_1d(xcoords);
  rjmcmc_destroy_array_1d(y);

  resultset1dfm_destroy(results);

  MPI_Finalize();

  return 0;
}

static double fx(double x)
{
  double y;

  if (x < 25.0) {
    y = 30.0 + 1.0*x;
  } else if (x < 50.0) {
    y = 100.0 - 2.0*x;
  } else if (x < 75.0) {
    y = -140.0 + 3.0*x;
  } else {
    y = 125.0 - 0.5*x;
  }

  return y;
}

static double likelihood(void *user_arg,
			 int npartitions,
			 const double *partitions,
			 int nglobalparameters,
			 const double *global_parameters,
			 int hierarchical,
			 int nhierarchicalparameters,
			 const double *hierarchical_parameters,
			 part1d_fm_hierarchical_likelihood_state_t *state,
			 part1d_fm_hierarchical_value_at_t value_at,
			 part1d_fm_hierarchical_value_at_t gradient_at,
			 double *logdetce)
{
  struct my_data *data = (struct my_data*)user_arg;
  int i;
  double dv;
  double sum;

  const double *local_parameters;
  double sigma;

  sum = 0.0;

  sigma = hierarchical_parameters[0];

  for (i = 0; i < xsamples; i ++) {

    local_parameters = value_at(state, data->x[i]);

    dv = data->sy[i] - local_parameters[0];

    sum += dv*dv;
  }

  if (hierarchical) {
    *logdetce = 2.0 * (double)xsamples * log(sigma);
  }

  return sum/(2.0 * sigma * sigma);
}

