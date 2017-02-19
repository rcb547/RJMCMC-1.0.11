
#include <stdio.h>
#include <stdlib.h>

#include <rjmcmc/rjmcmc.h>

#include <rjmcmc/resultsetfm.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_debug.h>
#include <rjmcmc/rjmcmc_defines.h>

#define RESULTSETFM_2DARRAY(el, a, b) \
  el = rjmcmc_create_array_2d(a, b); \
  if (el == NULL) {		     \
    return NULL;		     \
  }

#define RESULTSETFM_CONDITIONAL_2DARRAY(cond, el, a, b)	\
  el = NULL;						\
  if (cond) {						\
    el = rjmcmc_create_array_2d(a, b);			\
    if (el == NULL) {					\
      return NULL;					\
    }							\
  }

#define RESULTSETFM_1DARRAY(el, a)		\
  el = rjmcmc_create_array_1d(a);		\
  if (el == NULL) {				\
    return NULL;				\
  }
  
#define RESULTSETFM_CONDITIONAL_1DARRAY(cond, el, a)	\
  el = NULL;						\
  if (cond) {						\
    el = rjmcmc_create_array_1d(a);			\
    if (el == NULL) {					\
      return NULL;					\
    }							\
  }

struct _resultsetfm {
  
  int burnin;
  int total;
  int samples;

  int nparameters;
  const forwardmodelparameter_t *parameters;

  int nhierarchical_parameters;
  const forwardmodelparameter_t *hierarchical_parameters;

  int *propose;
  int *accept;

  /*
   * Misfit [total];
   */
  double *misfit;

  /*
   * Parameter values history [nparameters][total]
   */
  double **parameter_hist;

  /*
   * Hierachical parameter values history [nhiearchical_parameters][total]
   */
  double **hierarchical_parameter_hist;

  /*
   * mean of parameters [nparameters]
   */
  double *mean;

  /*
   * mode of parameters [nparameters]
   */
  double *mode;

  /*
   * median of parameters [nparameters]
   */
  double *median;

  /*
   * Credible interval of parameters [nparameters]
   */
  double credible_interval;
  double *credible_min;
  double *credible_max;

};
    
resultsetfm_t *
resultsetfm_create(int burnin,
		   int total,
		   int nparameters,
		   const forwardmodelparameter_t *parameters,
		   int nhierarchical_parameters,
		   const forwardmodelparameter_t *hierarchical_parameters,
		   int samples,
		   double credible_interval,
		   int results)
{
  resultsetfm_t *r;

  r = (resultsetfm_t*)malloc(sizeof(resultsetfm_t));
  if (r == NULL) {
    return NULL;
  }

  r->burnin = burnin;
  r->total = total;
  r->samples = samples;

  r->nparameters = nparameters;
  r->parameters = parameters;

  r->nhierarchical_parameters = nhierarchical_parameters;
  r->hierarchical_parameters = hierarchical_parameters;

  r->credible_interval = credible_interval;

  r->propose = rjmcmc_create_int_array_1d(nparameters + nhierarchical_parameters);
  if (r->propose == NULL) {
    return NULL;
  }

  r->accept = rjmcmc_create_int_array_1d(nparameters + nhierarchical_parameters);
  if (r->accept == NULL) {
    return NULL;
  }

  RESULTSETFM_1DARRAY(r->misfit, total);
  RESULTSETFM_2DARRAY(r->parameter_hist, nparameters, total);


  RESULTSETFM_CONDITIONAL_2DARRAY(nhierarchical_parameters > 0, 
				  r->hierarchical_parameter_hist, 
				  nhierarchical_parameters, total);

  RESULTSETFM_CONDITIONAL_1DARRAY((results & RESULTSETFM_MEAN) > 0,
				  r->mean, nparameters);

  RESULTSETFM_CONDITIONAL_1DARRAY((results & RESULTSETFM_MEDIAN) > 0,
				  r->median, nparameters);
  RESULTSETFM_CONDITIONAL_1DARRAY((results & RESULTSETFM_MODE) > 0,
				  r->mode, nparameters);

  RESULTSETFM_CONDITIONAL_1DARRAY((results & RESULTSETFM_CREDIBLE) > 0,
				  r->credible_min, nparameters);
  RESULTSETFM_CONDITIONAL_1DARRAY((results & RESULTSETFM_CREDIBLE) > 0,
				  r->credible_max, nparameters);

  return r;
      
}

void 
resultsetfm_destroy(resultsetfm_t *r)
{
  if (r != NULL) {

    rjmcmc_destroy_array_1d(r->credible_max);
    rjmcmc_destroy_array_1d(r->credible_min);
    rjmcmc_destroy_array_1d(r->mode);
    rjmcmc_destroy_array_1d(r->median);
    rjmcmc_destroy_array_1d(r->mean);

    rjmcmc_destroy_array_2d(r->nhierarchical_parameters, r->hierarchical_parameter_hist);

    rjmcmc_destroy_array_2d(r->nparameters, r->parameter_hist);
    rjmcmc_destroy_array_1d(r->misfit);

    rjmcmc_destroy_int_array_1d(r->accept);
    rjmcmc_destroy_int_array_1d(r->propose);
    
    free(r);
  }
}


/*
 * Sampling methods
 */

void
resultsetfm_propose(resultsetfm_t *r,
		    int p)
{
  RJMCMC_NULLCHECKVOID(r, "resultsetfm_propose: NULL results\n");

  RJMCMC_INDEXCHECKVOID(p, (r->nparameters + r->nhierarchical_parameters),
			"resultsetfm_propose: invalid index\n");

  r->propose[p] ++;
}

void
resultsetfm_accept(resultsetfm_t *r,
		   int p)
{
  RJMCMC_NULLCHECKVOID(r, "resultsetfm_accept: NULL results\n");

  RJMCMC_INDEXCHECKVOID(p, (r->nparameters + r->nhierarchical_parameters),
			"resultsetfm_accept: invalid index\n");

  r->accept[p] ++;
}

void
resultsetfm_sample(resultsetfm_t *r,
		   int i,
		   const double *v)
{
  int p;

  if (r == NULL) {
    rjmcmc_error("resultsetfm_sample: NULL results\n");
    return;
  }

  if (i < 0 || i >= r->total) {
    rjmcmc_error("resultsetfm_sample: invalid index\n");
    return;
  }

  for (p = 0; p < r->nparameters; p ++) {
    r->parameter_hist[p][i] = v[p];
  }
}

void 
resultsetfm_sample_misfit(resultsetfm_t *r,
			  int i,
			  double misfit)
{
  if (r == NULL) {
    rjmcmc_error("resultsetfm_sample_misfit: NULL results\n");
    return;
  }

  if (i < 0 || i >= r->total) {
    rjmcmc_error("resultsetfm_sample_misfit: invalid index\n");
    return;
  }

  r->misfit[i] = misfit;
}

void
resultsetfm_sample_hierarchical(resultsetfm_t *r,
				int i,
				double *hierarchical)
{
  int p;

  if (r == NULL) {
    rjmcmc_error("resultsetfm_sample_hierarchical: NULL results\n");
    return;
  }

  if (i < 0 || i >= r->total) {
    rjmcmc_error("resultsetfm_sample_hierarchical: invalid index\n");
    return;
  }

  for (p = 0; p < r->nhierarchical_parameters; p ++) {
    r->hierarchical_parameter_hist[p][i] = hierarchical[p];
  }
}



/*
 * Assemble after simulation (needed for histogram derived values, eg median)
 */

void
resultsetfm_assemble_results(resultsetfm_t *r)
{
  int i;
  int j;
  double t;

  int **hist;
  int conf_drop;

  if (r->mean != NULL) {
    for (i = 0; i < r->nparameters; i ++) {
      if (rjmcmc_mean_variance(r->parameter_hist[i] + r->burnin, 
			       r->total - r->burnin,
			       &(r->mean[i]),
			       &t) < 0) {
	rjmcmc_error("resultsetfm_assemble_results: failed to compute mean and variance\n");
	return;
      }
    }
  }

  hist = NULL;
  if (r->mode != NULL ||
      r->median != NULL ||
      r->credible_min != NULL ||
      r->credible_max != NULL) {
    
    /* We need to create do the histograms for each of the parameters */
    hist = rjmcmc_create_int_array_2d(r->nparameters, r->samples);
    if (hist == NULL) {
      rjmcmc_error("resultsetfm_assemble_results: failed to allocate memory for histogram\n");
      return;
    }

    for (i = 0; i < r->nparameters; i ++) {

      for (j = r->burnin; j < r->total; j ++) {

	hist[i][rjmcmc_map_to_index(r->parameter_hist[i][j],
				    r->parameters[i].fmin,
				    r->parameters[i].fmax,
				    r->samples)] ++;

      }
    }						      
  }

  if (r->mode != NULL) {
    for (i = 0; i < r->nparameters; i ++) {
      
      r->mode[i] = rjmcmc_mode_from_histogram(hist[i],
					      r->parameters[i].fmin,
					      r->parameters[i].fmax,
					      r->samples);
    }
  }

  if (r->median != NULL) {
    for (i = 0; i < r->nparameters; i ++) {
      
      r->median[i] = rjmcmc_median_from_histogram(hist[i],
						  r->parameters[i].fmin,
						  r->parameters[i].fmax,
						  r->samples);
    }
  }

  conf_drop = (int)((double)(r->total - r->burnin) * 
		    (1.0 - r->credible_interval)/2.0);

  if (r->credible_min != NULL) {
    for (i = 0; i < r->nparameters; i ++) {

      r->credible_min[i] = rjmcmc_head_from_histogram(hist[i],
							r->parameters[i].fmin,
							r->parameters[i].fmax,
							r->samples,
							conf_drop);

    }
  }

  if (r->credible_max != NULL) {
    for (i = 0; i < r->nparameters; i ++) {

      r->credible_max[i] = rjmcmc_tail_from_histogram(hist[i],
							r->parameters[i].fmin,
							r->parameters[i].fmax,
							r->samples,
							conf_drop);
    }
  }

  if (hist != NULL) {
    rjmcmc_destroy_int_array_2d(r->nparameters, hist);
  }
}

#if defined(HAVE_MPI_H)

void
MPI_resultsetfm_assemble_results(resultsetfm_t *r,
				 int mpisize,
				 int mpirank,
				 int root,
				 MPI_Comm comm)
{
  double *v;
  int vsize;
  int *iv;

  int i;
  int j;
  int nprocesses;

  double t;

  vsize = r->nparameters;
  v = rjmcmc_create_array_1d(vsize);
  if (v == NULL) {
    rjmcmc_error("MPI_resultsetfm_assemble_results: failed to allocate memory\n");
    return;
  }

  nprocesses = r->nparameters + r->nhierarchical_parameters;
  iv = rjmcmc_create_int_array_1d(nprocesses);
  if (iv == NULL) {
    rjmcmc_error("MPI_resultsetfm_assemble_results: failed to allocate memory\n");
    return;
  }

  /*
   * Aggregate the total acceptance/proposal counts
   */
  MPI_Reduce(r->propose, iv, nprocesses, MPI_INT, MPI_SUM, root, comm);
  if (mpirank == root) {
    for (j = 0; j < nprocesses; j ++) {
      r->propose[j] = iv[j];
    }
  }
  
  MPI_Reduce(r->accept, iv, nprocesses, MPI_INT, MPI_SUM, root, comm);
  if (mpirank == root) {
    for (j = 0; j < nprocesses; j ++) {
      r->accept[j] = iv[j];
    }
  }


  /*
   * Assemble the mean. First compute the mean then reduce by sum and 
   * divide by the number of processes.
   */
  if (r->mean != NULL) {
    for (i = 0; i < r->nparameters; i ++) {
      if (rjmcmc_mean_variance(r->parameter_hist[i] + r->burnin, 
			       r->total - r->burnin,
			       &(r->mean[i]),
			       &t) < 0) {
	rjmcmc_error("resultsetfm_assemble_results: failed to compute mean and variance\n");
	return;
      }
    }

    MPI_Reduce(r->mean, v, r->nparameters, MPI_DOUBLE, MPI_SUM, root, comm);
    if (mpirank == root) {
      for (j = 0; j < r->nparameters; j ++) {
	r->mean[j] = v[j]/(double)mpisize;
      }
    }
  }

}

#endif /* HAVE_MPI_H */

/*
 * Getting results
 */

int
resultsetfm_get_nparameters(resultsetfm_t *r)
{
  return r->nparameters + r->nhierarchical_parameters;
}
			    
const int *
resultsetfm_get_propose(resultsetfm_t *r,
			int *nprocesses)
{
  *nprocesses = resultsetfm_get_nparameters(r);
  return r->propose;
}

int
resultsetfm_get_propose_f(resultsetfm_t *r,
			  int maxsize,
			  int *propose)
{
  int n;
  int i;

  n = resultsetfm_get_nparameters(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    propose[i] = r->propose[i];
  }
  
  return n;
}

int
resultsetfm_get_propose_n(resultsetfm_t *r,
			  int process)
{
  int np = resultsetfm_get_nparameters(r);

  if (process >= 0 && process < np) {
    return r->propose[process];
  }

  return -1;
}

const int *
resultsetfm_get_accept(resultsetfm_t *r,
		       int *nprocesses)
{
  *nprocesses = resultsetfm_get_nparameters(r);
  return r->accept;
}

int
resultsetfm_get_accept_n(resultsetfm_t *r,
			 int process)
{
  int np = resultsetfm_get_nparameters(r);

  if (process >= 0 && process < np) {
    return r->accept[process];
  }

  return -1;
}

int
resultsetfm_get_accept_f(resultsetfm_t *r,
			 int maxsize,
			 int *accept)
{
  int n;
  int i;

  n = resultsetfm_get_nparameters(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    accept[i] = r->accept[i];
  }
  
  return n;
}

const double *
resultsetfm_get_misfit(resultsetfm_t *r)
{
  return r->misfit;
}

int
resultsetfm_get_misfit_f(resultsetfm_t *r,
			 int maxsize,
			 double *misfit)
{
  int n;
  int i;

  n = r->total;
  if (n > maxsize) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    misfit[i] = r->misfit[i];
  }

  return n;
}

const double *
resultsetfm_get_parameter_history(resultsetfm_t *r, int p)
{
  if (r == NULL) {
    rjmcmc_error("resultsetfm_get_parameter_history: NULL results\n");
    return NULL;
  }

  if (p < 0 || p >= r->nparameters) {
    rjmcmc_error("resultsetfm_get_parameter_history: invalid index\n");
    return NULL;
  }

  return r->parameter_hist[p];
}

int
resultsetfm_get_parameter_history_f(resultsetfm_t *r, int p, int maxsize, double *phistory)
{
  int n;
  int i;

  if (r == NULL) {
    rjmcmc_error("resultsetfm_get_parameter_history_f: NULL results\n");
    return -1;
  }

  if (p < 0 || p >= r->nparameters) {
    rjmcmc_error("resultsetfm_get_parameter_history_f: invalid index\n");
    return -1;
  }

  n = r->total;
  if (n > maxsize) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    phistory[i] = r->parameter_hist[p][i];
  }

  return n;
}

const double *
resultsetfm_get_hierarchical_parameter_history(resultsetfm_t *r, int p)
{
  return r->hierarchical_parameter_hist[p];
}

int
resultsetfm_get_hierarchical_parameter_history_f(resultsetfm_t *r,
						 int p,
						 int maxsize,
						 double *hphistory)
{
  int n;
  int i;

  if (r == NULL) {
    rjmcmc_error("resultsetfm_get_hierarchical_parameter_history_f: NULL results\n");
    return -1;
  }

  if (p < 0 || p >= r->nhierarchical_parameters) {
    rjmcmc_error("resultsetfm_get_hierarchical_parameter_history_f: invalid index\n");
    return -1;
  }

  n = r->total;
  if (n > maxsize) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    hphistory[i] = r->hierarchical_parameter_hist[p][i];
  }

  return n;
}

double 
resultsetfm_get_parameter_mean(resultsetfm_t *r, int p)
{
  if (r->mean != NULL) {
    return r->mean[p];
  } 

  return 0.0;
}

double 
resultsetfm_get_parameter_mode(resultsetfm_t *r, int p)
{
  if (r->mode != NULL) {
    return r->mode[p];
  }
  
  return 0.0;
}

double 
resultsetfm_get_parameter_median(resultsetfm_t *r, int p)
{
  if (r->median != NULL) {
    return r->median[p];
  }

  return 0.0;
}

double 
resultsetfm_get_parameter_credible_min(resultsetfm_t *r, int p)
{
  if (r->credible_min != NULL) {
    return r->credible_min[p];
  }

  return 0.0;
}

double 
resultsetfm_get_parameter_credible_max(resultsetfm_t *r, int p)
{
  if (r->credible_max != NULL) {
    return r->credible_max[p];
  }

  return 0.0;
}

