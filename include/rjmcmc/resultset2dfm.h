#ifndef resultset2dfm_h
#define resultset2dfm_h

#include <rjmcmc/rjmcmc.h>
#include <rjmcmc/forwardmodelparameter.h>

#if defined(HAVE_MPI_H)
#include <mpi.h>
#endif /* HAVE_MPI_H */

typedef enum {
  RESULTSET2DFM_MEAN       = 0x01,
  RESULTSET2DFM_MEDIAN     = 0x02,
  RESULTSET2DFM_MODE       = 0x04,
  RESULTSET2DFM_CREDIBLE   = 0x08,
} resultset2dfm_result_t;

typedef enum {
  RESULTSET2DFM_BINARY     = 0x01,
  RESULTSET2DFM_JSON       = 0x02
} resultset2dfm_format_t;

typedef struct _resultset2dfm resultset2dfm_t;

resultset2dfm_t *
resultset2dfm_create(int burnin,
		     int total,
		     int thin,
		     int nglobalparameters,
		     const forwardmodelparameter_t *global_parameters,
		     int nlocalparameters,
		     const forwardmodelparameter_t *local_parameters,
		     int nhierarchicalparameters,
		     int xsamples,
		     int ysamples,
		     int zsamples,
		     int max_partitions,
		     double xmin,
		     double xmax,
		     double ymin,
		     double ymax,
		     int nprocesses,
		     double credible_interval,
		     int results);

void
resultset2dfm_destroy(resultset2dfm_t *r);

/*
 * Save/Restore
 */
int resultset2dfm_save(const resultset2dfm_t *r,
		       const char *filename,
		       int fmt);

resultset2dfm_t *resultset2dfm_load(const char *filename,
				    int fmt,
				    int *results,
				    int *burnin,
				    int *total,
				    int *thin,
				    int *xsamples,
				    int *ysamples,
				    int *zsamples,
				    int *nprocesses,
				    int *maxpartitions,
				    double *xmin,
				    double *xmax,
				    double *ymin,
				    double *ymax,
				    double *credibleinterval,
				    int *nhierarchical,
				    int *nglobal,
				    int *nlocal,
				    forwardmodelparameter_t **global_parameters,
				    forwardmodelparameter_t **local_parameters);
		       
		       
/*
 * Sampling methods
 */

void 
resultset2dfm_propose(resultset2dfm_t *r,
		      int p);

void 
resultset2dfm_accept(resultset2dfm_t *r,
		     int p);

void
resultset2dfm_sample_global_parameter(resultset2dfm_t *r,
				      int i,
				      int gi,
				      double g);

void
resultset2dfm_sample_hierarchical_parameter(resultset2dfm_t *r,
					    int i,
					    int hi,
					    double v);

void
resultset2dfm_sample_local_parameter(resultset2dfm_t *r,
				     int i,
				     int li,
				     double **l);

void 
resultset2dfm_sample_misfit(resultset2dfm_t *r,
			    int i,
			    double misfit);

void
resultset2dfm_sample_npartitions(resultset2dfm_t *r,
				 int i,
				 int npartitions);

void
resultset2dfm_sample_centre(resultset2dfm_t *r,
			    double x,
			    double y);

/*
 * Assembling results
 */

void
resultset2dfm_assemble_results(resultset2dfm_t *r);

#if defined(HAVE_MPI_H)
void
MPI_resultset2dfm_assemble_results(resultset2dfm_t *r,
				   int mpisize,
				   int mpirank,
				   int root,
				   MPI_Comm comm);
#endif /* HAVE_MPI_H */

/*
 * Getting results
 */

int
resultset2dfm_get_nparameters(resultset2dfm_t *r);

const int *
resultset2dfm_get_propose(resultset2dfm_t *r,
			  int *nprocesses);

int
resultset2dfm_get_propose_f(resultset2dfm_t *r,
			    int maxsize,
			    int *propose);

const int *
resultset2dfm_get_accept(resultset2dfm_t *r,
			 int *nprocesses);

const double *
resultset2dfm_get_misfit(resultset2dfm_t *r);

int
resultset2dfm_get_misfit_f(resultset2dfm_t *r,
			   int maxsize,
			   double *misfit);

const double *
resultset2dfm_get_global_parameter(resultset2dfm_t *r, int gi);

int
resultset2dfm_get_global_parameter_f(resultset2dfm_t *r,
				     int gi,
				     int maxsize,
				     double *global_parameter);

const double *
resultset2dfm_get_hierarchical_parameter(resultset2dfm_t *r, int hi);

int
resultset2dfm_get_hierarchical_parameter_f(resultset2dfm_t *r,
					   int hi,
					   int maxsize,
					   double *hierarchical_parameter);

const double **
resultset2dfm_get_local_parameter_mean(resultset2dfm_t *r, int li);

int
resultset2dfm_get_local_parameter_mean_f(resultset2dfm_t *r,
					 int li,
					 int xsamples,
					 int ysamples,
					 double *mean);

const double **
resultset2dfm_get_local_parameter_variance(resultset2dfm_t *r, int li);

int
resultset2dfm_get_local_parameter_variance_f(resultset2dfm_t *r,
					     int li,
					     int xsamples,
					     int ysamples,
					     double *variance);

const double **
resultset2dfm_get_local_parameter_mode(resultset2dfm_t *r, int li);

int
resultset2dfm_get_local_parameter_mode_f(resultset2dfm_t *r,
					 int li,
					 int xsamples,
					 int ysamples,
					 double *mode);
const double **
resultset2dfm_get_local_parameter_median(resultset2dfm_t *r, int li);

int
resultset2dfm_get_local_parameter_median_f(resultset2dfm_t *r,
					   int li,
					   int xsamples,
					   int ysamples,
					   double *median);

const double **
resultset2dfm_get_local_parameter_credible_min(resultset2dfm_t *r, int li);

int
resultset2dfm_get_local_parameter_credible_min_f(resultset2dfm_t *r,
						 int li,
						 int xsamples,
						 int ysamples,
						 double *credible_min);

const double **
resultset2dfm_get_local_parameter_credible_max(resultset2dfm_t *r, int li);

int
resultset2dfm_get_local_parameter_credible_max_f(resultset2dfm_t *r,
						 int li,
						 int xsamples,
						 int ysamples,
						 double *credible_max);

const int *
resultset2dfm_get_partitions(resultset2dfm_t *r);

int
resultset2dfm_get_partitions_f(resultset2dfm_t *r,
			       int maxsize,
			       int *npartitions);

const int **
resultset2dfm_get_centres(resultset2dfm_t *r);

void
resultset2dfm_fill_xcoord_vector(resultset2dfm_t *r,
				 double *x,
				 int *l);

int
resultset2dfm_fill_xcoord_vector_f(resultset2dfm_t *r,
				   int maxsize,
				   double *x);

void
resultset2dfm_fill_ycoord_vector(resultset2dfm_t *r,
				 double *y,
				 int *l);

int
resultset2dfm_fill_ycoord_vector_f(resultset2dfm_t *r,
				   int maxsize,
				   double *y);

#endif /* resultset2dfm_h */
