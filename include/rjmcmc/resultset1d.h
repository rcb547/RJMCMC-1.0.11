#ifndef resultset1d_h
#define resultset1d_h

#include <rjmcmc/rjmcmc_config.h>

#if defined(HAVE_MPI_H)
#include <mpi.h>
#endif 

typedef enum {
  RESULTSET1D_MEAN       = 0x01,
  RESULTSET1D_MEDIAN     = 0x02,
  RESULTSET1D_MODE       = 0x04,
  RESULTSET1D_CREDIBLE = 0x08,

  RESULTSET1D_GRADIENT   = 0x10
} resultset1d_result_t;

typedef struct _resultset1d resultset1d_t;

resultset1d_t *
resultset1d_create(int burnin,
		   int total,
		   int xsamples,
		   int ysamples,
		   int nprocesses,
		   int maxpartitions,
		   int maxorder,
		   double xmin,
		   double xmax,
		   double ymin,
		   double ymax,
		   double credible_interval,
		   int results);

void 
resultset1d_destroy(resultset1d_t *r);


/*
 * Sampling methods
 */

void
resultset1d_propose(resultset1d_t *r,
		    int p);

void
resultset1d_accept(resultset1d_t *r,
		   int p);


void
resultset1d_sample(resultset1d_t *r,
		   int i,
		   const double *v);

void
resultset1d_sample_gradient(resultset1d_t *r,
			    int i,
			    const double *dv);

void 
resultset1d_sample_misfit(resultset1d_t *r,
			  int i,
			  double misfit);

void 
resultset1d_sample_order(resultset1d_t *r,
			 int i,
			 int order);

void 
resultset1d_sample_lambda(resultset1d_t *r,
			 int i,
			 double lambda);

void
resultset1d_sample_npartitions(resultset1d_t *r,
			       int i,
			       int npartitions);

void
resultset1d_sample_partition_x(resultset1d_t *r,
			       int i,
			       double x);

/*
 * Assemble after simulation (needed for histogram derived values, eg median)
 */

void
resultset1d_assemble_results(resultset1d_t *r);

#if defined(HAVE_MPI_H)
void
MPI_resultset1d_assemble_results(resultset1d_t *r,
				 int mpisize,
				 int mpirank,
				 int root,
				 MPI_Comm comm);
#endif /* HAVE_OPENMPI_MPI_H */

/*
 * Getting results
 */

const int *
resultset1d_get_propose(resultset1d_t *r,
			int *nprocesses);

const int *
resultset1d_get_accept(resultset1d_t *r,
		       int *nprocesses);

const double *
resultset1d_get_misfit(resultset1d_t *r);

const int *
resultset1d_get_order(resultset1d_t *r);

const double *
resultset1d_get_lambda(resultset1d_t *r);

const int *
resultset1d_get_partitions(resultset1d_t *r);

const int *
resultset1d_get_partition_x_histogram(resultset1d_t *r);

const double *
resultset1d_get_mean(resultset1d_t *r);

const double *
resultset1d_get_median(resultset1d_t *r);

const double *
resultset1d_get_mode(resultset1d_t *r);

const double *
resultset1d_get_credible_min(resultset1d_t *r);
		     
const double *
resultset1d_get_credible_max(resultset1d_t *r);

const int **
resultset1d_get_histogram(resultset1d_t *r);

const double *
resultset1d_get_gradient(resultset1d_t *r);

const double *
resultset1d_get_gradient_credible_min(resultset1d_t *r);

const double *
resultset1d_get_gradient_credible_max(resultset1d_t *r);

int
resultset1d_get_total(resultset1d_t *r);

int 
resultset1d_get_max_partitions(resultset1d_t *r);

int 
resultset1d_get_xsamples(resultset1d_t *r);

int
resultset1d_get_ysamples(resultset1d_t *r);

int 
resultset1d_get_max_order(resultset1d_t *r);

void
resultset1d_fill_xcoord_vector(resultset1d_t *r,
			       double *x);

void
resultset1d_fill_ycoord_vector(resultset1d_t *r,
			       double *y);

#endif /* resultset_h */
