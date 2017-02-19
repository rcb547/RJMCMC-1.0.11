#ifndef resultset1dfm_h
#define resultset1dfm_h

#include <rjmcmc/rjmcmc.h>
#include <rjmcmc/forwardmodelparameter.h>

#if defined(HAVE_MPI_H)
#include <mpi.h>
#endif /* HAVE_MPI_H */

typedef enum {
  RESULTSET1DFM_MEAN       = 0x01,
  RESULTSET1DFM_MEDIAN     = 0x02,
  RESULTSET1DFM_MODE       = 0x04,
  RESULTSET1DFM_CREDIBLE   = 0x08,
} resultset1dfm_result_t;

typedef struct _resultset1dfm resultset1dfm_t;

resultset1dfm_t *
resultset1dfm_create(int burnin,
		     int total,
		     int nglobalparameters,
		     const forwardmodelparameter_t *global_parameters,
		     int nlocalparameters,
		     const forwardmodelparameter_t *local_parameters,
		     int nhierarchicalparameters,
		     int xsamples,
		     int ysamples,
		     int max_partitions,
		     double xmin,
		     double xmax,
		     int nprocesses,
		     double credible_interval,
		     int results);

void
resultset1dfm_destroy(resultset1dfm_t *r);

/*
 * Sampling methods
 */

void 
resultset1dfm_propose(resultset1dfm_t *r,
		      int p);

void 
resultset1dfm_accept(resultset1dfm_t *r,
		     int p);

void
resultset1dfm_propose_local_value(resultset1dfm_t *r,
				  int li);

void
resultset1dfm_accept_local_value(resultset1dfm_t *r,
				 int li);



void
resultset1dfm_sample_global_parameter(resultset1dfm_t *r,
				      int i,
				      int gi,
				      double g);

void
resultset1dfm_sample_local_parameter(resultset1dfm_t *r,
				     int i,
				     int li,
				     double *l);

void
resultset1dfm_sample_npartitions(resultset1dfm_t *r,
				 int i,
				 int np);

void
resultset1dfm_sample_partition_x(resultset1dfm_t *r,
				 double x);

void 
resultset1dfm_sample_misfit(resultset1dfm_t *r,
			    int i,
			    double misfit);

void 
resultset1dfm_sample_hierarchical(resultset1dfm_t *r,
			   int i,
			   int si,
			   double hierarchical);

/*
 * Assembling results
 */

void
resultset1dfm_assemble_results(resultset1dfm_t *r);

#if defined(HAVE_MPI_H)
void
MPI_resultset1dfm_assemble_results(resultset1dfm_t *r,
				   int mpisize,
				   int mpirank,
				   int root,
				   MPI_Comm comm);
#endif /* HAVE_MPI_H */

/*
 * Getting results
 */

int 
resultset1dfm_get_max_partitions(resultset1dfm_t *r);

int
resultset1dfm_get_total(resultset1dfm_t *r);

int
resultset1dfm_get_xsamples(resultset1dfm_t *r);

int
resultset1dfm_get_nparameters(resultset1dfm_t *r);

const int *
resultset1dfm_get_propose(resultset1dfm_t *r,
			  int *nprocesses);

int
resultset1dfm_get_propose_f(resultset1dfm_t *r,
			    int maxsize,
			    int *propose);

const int *
resultset1dfm_get_accept(resultset1dfm_t *r,
			 int *nprocesses);

int
resultset1dfm_get_accept_f(resultset1dfm_t *r,
			   int maxsize,
			   int *accept);

const double *
resultset1dfm_get_misfit(resultset1dfm_t *r);

int
resultset1dfm_get_misfit_f(resultset1dfm_t *r,
			   int maxsize,
			   double *misfit);

const double *
resultset1dfm_get_global_parameter(resultset1dfm_t *r, int gi);

int
resultset1dfm_get_global_parameter_f(resultset1dfm_t *r, 
				     int gi,
				     int maxsize,
				     double *global);

const double *
resultset1dfm_get_local_parameter_mean(resultset1dfm_t *r, int li);

int resultset1dfm_get_local_parameter_mean_f(resultset1dfm_t *r,
					     int li,
					     int maxsize,
					     double *mean);

const double *
resultset1dfm_get_local_parameter_median(resultset1dfm_t *r, int li);

int resultset1dfm_get_local_parameter_median_f(resultset1dfm_t *r,
					       int li,
					       int maxsize,
					       double *median);

const double *
resultset1dfm_get_local_parameter_mode(resultset1dfm_t *r, int li);

int resultset1dfm_get_local_parameter_mode_f(resultset1dfm_t *r,
					     int li,
					     int maxsize,
					     double *mode);

const double *
resultset1dfm_get_local_parameter_credible_min(resultset1dfm_t *r, int li);

int resultset1dfm_get_local_parameter_credible_min_f(resultset1dfm_t *r,
						     int li,
						     int maxsize,
						     double *credible_min);

const double *
resultset1dfm_get_local_parameter_credible_max(resultset1dfm_t *r, int li);

int resultset1dfm_get_local_parameter_credible_max_f(resultset1dfm_t *r,
						     int li,
						     int maxsize,
						     double *credible_max);

const int **
resultset1dfm_get_local_parameter_histogram(resultset1dfm_t *r, int li);

int
resultset1dfm_get_local_parameter_histogram_f(resultset1dfm_t *r, 
					      int li,
					      int xsamples,
					      int ysamples,
					      int *histogram);

const int *
resultset1dfm_get_partitions(resultset1dfm_t *r);

int
resultset1dfm_get_partitions_f(resultset1dfm_t *r,
			       int maxsize,
			       int *partitions);

const int *
resultset1dfm_get_partition_x_histogram(resultset1dfm_t *r);

int
resultset1dfm_get_partition_x_histogram_f(resultset1dfm_t *r,
					  int maxsize,
					  int *histogram);

const double *
resultset1dfm_get_hierarchical(resultset1dfm_t *r,
			       int si);

int 
resultset1dfm_get_hierarchical_f(resultset1dfm_t *r,
				 int hi,
				 int maxsize,
				 double *hierarchical);

const int **
resultset1dfm_get_local_parameter_histogram(resultset1dfm_t *r,
					    int li);

int
resultset1dfm_get_local_parameter_histogram_f(resultset1dfm_t *r,
					      int li,
					      int xsamples,
					      int ysamples,
					      int *histogram);

void
resultset1dfm_fill_xcoord_vector(resultset1dfm_t *r,
				 double *x,
				 int *l);

int
resultset1dfm_get_xcoord_vector_f(resultset1dfm_t *r,
				  int maxsize,
				  double *x);

void 
resultset1dfm_fill_ycoord_vector(resultset1dfm_t *r,
				 int li,
				 double *y,
				 int *l);

int
resultset1dfm_get_ycoord_vector_f(resultset1dfm_t *r,
				  int li,
				  int maxsize,
				  double *y);


#endif /* resultset1dfm_h */
