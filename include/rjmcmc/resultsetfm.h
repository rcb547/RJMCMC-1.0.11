#ifndef resultsetfm_h
#define resultsetfm_h

#include <rjmcmc/rjmcmc_config.h>
#include <rjmcmc/forwardmodelparameter.h>

#if defined(HAVE_MPI_H)
#include <mpi.h>
#endif 

typedef enum {
  RESULTSETFM_MEAN       = 0x01,
  RESULTSETFM_MEDIAN     = 0x02,
  RESULTSETFM_MODE       = 0x04,
  RESULTSETFM_CREDIBLE   = 0x08,
} resultsetfm_result_t;


typedef struct _resultsetfm resultsetfm_t;

resultsetfm_t *
resultsetfm_create(int burnin,
		   int total,
		   int nparameters,
		   const forwardmodelparameter_t *parameters,
		   int nhierarchicalparameters,
		   const forwardmodelparameter_t *hierarchical_parameters,
		   int samples,
		   double credible_interval,
		   int results);

void 
resultsetfm_destroy(resultsetfm_t *r);


/*
 * Sampling methods
 */

void
resultsetfm_propose(resultsetfm_t *r,
		    int p);

void
resultsetfm_accept(resultsetfm_t *r,
		   int p);

void 
resultsetfm_sample_misfit(resultsetfm_t *r,
			  int i,
			  double misfit);

void
resultsetfm_sample(resultsetfm_t *r,
		   int i,
		   const double *v);

void
resultsetfm_sample_hierarchical(resultsetfm_t *r,
				int i,
				double *hierarchical);

/*
 * Assemble after simulation (needed for histogram derived values, eg median)
 */

void
resultsetfm_assemble_results(resultsetfm_t *r);

#if defined(HAVE_MPI_H)

void
MPI_resultsetfm_assemble_results(resultsetfm_t *r,
				 int mpisize,
				 int mpirank,
				 int root,
				 MPI_Comm comm);

#endif /* HAVE_MPI_H */

/*
 * Getting results
 */

int
resultsetfm_get_nprocesses(resultsetfm_t *r);

const int *
resultsetfm_get_propose(resultsetfm_t *r,
			int *nprocesses);

int
resultsetfm_get_propose_f(resultsetfm_t *r,
			  int maxsize,
			  int *propose);

int 
resultsetfm_get_propose_n(resultsetfm_t *r,
			  int process);

const int *
resultsetfm_get_accept(resultsetfm_t *r,
		       int *nprocesses);

int
resultsetfm_get_accept_f(resultsetfm_t *r,
			 int maxsize,
			 int *accept);

int
resultsetfm_get_accept_n(resultsetfm_t *r,
			 int process);

const double *
resultsetfm_get_misfit(resultsetfm_t *r);

int
resultsetfm_get_misfit_f(resultsetfm_t *r,
			 int maxsize,
			 double *misfit);


const double *
resultsetfm_get_parameter_history(resultsetfm_t *r, int p);

int
resultsetfm_get_parameter_history_f(resultsetfm_t *r,
				    int p,
				    int maxsize,
				    double *phistory);

const double *
resultsetfm_get_hierarchical_parameter_history(resultsetfm_t *r,
					       int p);

int
resultsetfm_get_hierarchical_parameter_history_f(resultsetfm_t *r,
						 int p,
						 int maxsize,
						 double *hphistory);

double 
resultsetfm_get_parameter_mean(resultsetfm_t *r, int p);

double
resultsetfm_get_parameter_variance(resultsetfm_t *r, int p);

double 
resultsetfm_get_parameter_mode(resultsetfm_t *r, int p);

double 
resultsetfm_get_parameter_median(resultsetfm_t *r, int p);

double 
resultsetfm_get_parameter_credible_min(resultsetfm_t *r, int p);

double 
resultsetfm_get_parameter_credible_max(resultsetfm_t *r, int p);

#endif /* resultset_h */
