

#include <stdio.h>
#include <stdlib.h>

#include <rjmcmc/resultset1dfm.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_defines.h>
#include <rjmcmc/rjmcmc_debug.h>

struct _resultset1dfm {

  int results;

  int burnin;
  int total;

  int xsamples;
  int ysamples;

  int nglobalparameters;
  const forwardmodelparameter_t *global_parameters;

  int nlocalparameters;
  const forwardmodelparameter_t *local_parameters;

  int maxpartitions;

  double xmin;
  double xmax;

  int nprocesses;
  int *propose;
  int *accept;

  int *propose_local;
  int *accept_local;

  /* 
   * Misfit [total]
   */
  double *misfit;

  /*
   * N Partitions [total];
   */
  int *partitions;

  /*
   * Partition x histogram [xsamples]
   */
  int *partition_x_hist;

  /*
   * Global Parameter History [nglobalparameters][total];
   */
  double **global;

  /*
   * Local Parameter Mean(s) [nlocalparameters][xsamples]
   */
  double **local_mean;

  /*
   * Hierarchical parameter histories [nhierarchical][total]
   */
  int nhierarchical;
  double **hierarchical;

  /*
   * Local Parameter Histogram [nlocalparameters][xsamples][ysamples]
   * (only used if mode/median/credible requested)
   */
  int ***histogram;

  /*
   * Local median [nlocalparameters][xsamples]
   */
  double **local_median;
  
  /*
   * Local mode [nlocalparameters][xsamples]
   */
  double **local_mode;

  /*
   * Credible intervals [nlocalparameters][xsamples]
   */
  double credible_interval;
  double **local_cred_min;
  double **local_cred_max;
};

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
		     int maxpartitions,
		     double xmin,
		     double xmax,
		     int nprocesses,
		     double credible_interval,
		     int results)
{
  resultset1dfm_t *r;

  r = (resultset1dfm_t *)malloc(sizeof(resultset1dfm_t));
  if (r == NULL) {
    return NULL;
  }

  r->results = results;
  r->burnin = burnin;
  r->total = total;
  
  r->xsamples = xsamples;
  r->ysamples = ysamples;

  r->nglobalparameters = nglobalparameters;
  r->global_parameters = global_parameters;

  r->nlocalparameters = nlocalparameters;
  r->local_parameters = local_parameters;

  r->maxpartitions = maxpartitions;

  r->xmin = xmin;
  r->xmax = xmax;

  r->nprocesses = nprocesses;

  r->propose = rjmcmc_create_int_array_1d(nprocesses);
  if (r->propose == NULL) {
    return NULL;
  }

  r->accept = rjmcmc_create_int_array_1d(nprocesses);
  if (r->accept == NULL) {
    return NULL;
  }

  r->propose_local = rjmcmc_create_int_array_1d(nlocalparameters);
  if (r->propose_local == NULL) {
    return NULL;
  }

  r->accept_local = rjmcmc_create_int_array_1d(nlocalparameters);
  if (r->accept_local == NULL) {
    return NULL;
  }

  r->misfit = rjmcmc_create_array_1d(total);
  if (r->misfit == NULL) {
    return NULL;
  }

  r->partitions = NULL;
  if (maxpartitions > 0) {
    r->partitions = rjmcmc_create_int_array_1d(total);
    if (r->partitions == NULL) {
      return NULL;
    }
  }

  r->partition_x_hist = NULL;
  if (maxpartitions > 0) {
    r->partition_x_hist = rjmcmc_create_int_array_1d(xsamples);
    if (r->partition_x_hist == NULL) {
      return NULL;
    }
  }

  r->global = NULL;
  if (nglobalparameters > 0) {
    r->global = rjmcmc_create_array_2d(nglobalparameters, total);
    if (r->global == NULL) {
      return NULL;
    }
  }

  r->local_mean = rjmcmc_create_array_2d(nlocalparameters, xsamples);
  if (r->local_mean == NULL) {
    return NULL;
  }

  r->nhierarchical = nhierarchicalparameters;
  r->hierarchical = NULL; 
  if (nhierarchicalparameters > 0) {
    r->hierarchical = rjmcmc_create_array_2d(nhierarchicalparameters, total);
    if (r->hierarchical == NULL) {
      return NULL;
    }
  }

  r->histogram = NULL;
  if ((results & (RESULTSET1DFM_MODE | 
		  RESULTSET1DFM_MEDIAN | 
		  RESULTSET1DFM_CREDIBLE)) > 0) {
    r->histogram = rjmcmc_create_int_array_3d(nlocalparameters, 
					      xsamples, 
					      ysamples);
    if (r->histogram == NULL) {
      return NULL;
    }
  }

  r->local_median = NULL;
  if ((results & RESULTSET1DFM_MEDIAN) > 0) {
    r->local_median = rjmcmc_create_array_2d(nlocalparameters, xsamples);
    if (r->local_median == NULL) {
      return NULL;
    }
  }

  r->local_mode = NULL;
  if ((results & RESULTSET1DFM_MODE) > 0) {
    r->local_mode = rjmcmc_create_array_2d(nlocalparameters, xsamples);
    if (r->local_mode == NULL) {
      return NULL;
    }
  }

  r->credible_interval = credible_interval;
  r->local_cred_min = NULL;
  r->local_cred_max = NULL;
  if ((results & RESULTSET1DFM_CREDIBLE) > 0) {
    r->local_cred_min = rjmcmc_create_array_2d(nlocalparameters, xsamples);
    if (r->local_cred_min == NULL) {
      return NULL;
    }
    r->local_cred_max = rjmcmc_create_array_2d(nlocalparameters, xsamples);
    if (r->local_cred_max == NULL) {
      return NULL;
    }
  }

  return r;
}

void 
resultset1dfm_destroy(resultset1dfm_t *r)
{
  if (r != NULL) {

    rjmcmc_destroy_int_array_1d(r->propose);
    rjmcmc_destroy_int_array_1d(r->accept);
    
    rjmcmc_destroy_int_array_1d(r->propose_local);
    rjmcmc_destroy_int_array_1d(r->accept_local);

    rjmcmc_destroy_array_1d(r->misfit);
    rjmcmc_destroy_int_array_1d(r->partitions);
    rjmcmc_destroy_int_array_1d(r->partition_x_hist);
    
    rjmcmc_destroy_array_2d(r->nglobalparameters, r->global);

    rjmcmc_destroy_array_2d(r->nlocalparameters, r->local_mean);

    rjmcmc_destroy_array_2d(r->nhierarchical, r->hierarchical);

    rjmcmc_destroy_int_array_3d(r->nlocalparameters, r->xsamples, r->histogram);

    rjmcmc_destroy_array_2d(r->nlocalparameters, r->local_median);
    rjmcmc_destroy_array_2d(r->nlocalparameters, r->local_mode);
    rjmcmc_destroy_array_2d(r->nlocalparameters, r->local_cred_min);
    rjmcmc_destroy_array_2d(r->nlocalparameters, r->local_cred_max);

    free(r);
  }
}

void
resultset1dfm_propose(resultset1dfm_t *r,
		      int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses, 
			"resultset1dfm_propose: invalid index\n");

  r->propose[p] ++;
}

void
resultset1dfm_accept(resultset1dfm_t *r,
		     int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses, 
			"resultset1dfm_accept: invalid index\n");

  r->accept[p] ++;
}

void
resultset1dfm_propose_local_value(resultset1dfm_t *r,
				  int li)
{
  RJMCMC_INDEXCHECKVOID(li, 
			r->nlocalparameters, 
			"resultset1dfm_propose_local_value: invalid index\n");

  r->propose_local[li] ++;
}


void
resultset1dfm_accept_local_value(resultset1dfm_t *r,
				 int li)
{
  RJMCMC_INDEXCHECKVOID(li, 
			r->nlocalparameters, 
			"resultset1dfm_accept_local_value: invalid index\n");

  r->accept_local[li] ++;
}


void 
resultset1dfm_sample_global_parameter(resultset1dfm_t *r,
				      int i,
				      int gi,
				      double g)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset1dfm_sample_global_parameter: invalid index\n");

  RJMCMC_INDEXCHECKVOID(gi, 
			r->nglobalparameters, 
			"resultset1dfm_sample_global_parameter: invalid index\n");

  r->global[gi][i] = g;
}

void
resultset1dfm_sample_local_parameter(resultset1dfm_t *r,
				     int i,
				     int li,
				     double *l)
{
  int j;

  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset1dfm_sample_local_parameter: invalid index\n");

  RJMCMC_INDEXCHECKVOID(li, 
			r->nlocalparameters, 
			"resultset1dfm_sample_local_parameter: invalid index\n");

  if (i >= r->burnin) {
    for (j = 0; j < r->xsamples; j ++) {
      r->local_mean[li][j] += l[j];
    }

    if (r->histogram != NULL) {
      for (j = 0; j < r->xsamples; j ++) {
	r->histogram[li][j][rjmcmc_map_to_index(l[j], 
						r->local_parameters[li].fmin,
						r->local_parameters[li].fmax,
						r->ysamples)] ++;
      }
    }
  }
}

void 
resultset1dfm_sample_misfit(resultset1dfm_t *r,
			    int i,
			    double misfit)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset1dfm_sample_misfit: invalid index\n");

  r->misfit[i] = misfit;
}

void
resultset1dfm_sample_npartitions(resultset1dfm_t *r,
				 int i,
				 int npartitions)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset1dfm_sample_npartitions: invalid index\n");

  r->partitions[i] = npartitions;
}

void
resultset1dfm_sample_partition_x(resultset1dfm_t *r,
				 double x)
{
  if (x > r->xmin && x < r->xmax) {
    /*
     * The range check here ensures that we don't count the end partitions
     * in the final histogram.
     */

    r->partition_x_hist[rjmcmc_map_to_index(x, 
					    r->xmin, 
					    r->xmax, 
					    r->xsamples)] ++;
  }
}

void
resultset1dfm_sample_hierarchical(resultset1dfm_t *r,
				  int i,
				  int si,
				  double hierarchical)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset1dfm_sample_hierarchical: invalid index\n");
  
  RJMCMC_INDEXCHECKVOID(si, 
			r->nhierarchical, 
			"resultset1dfm_sample_hierarchical: invalid index\n");

  r->hierarchical[si][i] = hierarchical;
}


/*
 * Assembling
 */
void
resultset1dfm_assemble_results(resultset1dfm_t *r)
{
  int i;
  int j;
  double denom;
  int cred_samples;

  denom = r->total - r->burnin;

  if (r->local_mean) {
    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) { 
	r->local_mean[i][j] /= denom;
      }
    }
  }

  if (r->local_median) {
    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) {
	r->local_median[i][j] = 
	  rjmcmc_median_from_histogram(r->histogram[i][j],
				       r->local_parameters[i].fmin,
				       r->local_parameters[i].fmax,
				       r->ysamples);
      }
    }
  }

  if (r->local_mode) {
    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) {
	r->local_mode[i][j] = 
	  rjmcmc_mode_from_histogram(r->histogram[i][j],
				     r->local_parameters[i].fmin,
				     r->local_parameters[i].fmax,
				     r->ysamples);
      }
    }
  }

  if (r->local_cred_min &&
      r->local_cred_max) {

    cred_samples = (double)(r->total - r->burnin) * (1.0 - r->credible_interval)/2.0;
    
    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) {

	r->local_cred_min[i][j] = 
	  rjmcmc_head_from_histogram(r->histogram[i][j],
				     r->local_parameters[i].fmin,
				     r->local_parameters[i].fmax,
				     r->ysamples,
				     cred_samples);
							     
	r->local_cred_max[i][j] = 
	  rjmcmc_tail_from_histogram(r->histogram[i][j],
				     r->local_parameters[i].fmin,
				     r->local_parameters[i].fmax,
				     r->ysamples,
				     cred_samples);

      }
    }
  }
}

#if defined(HAVE_MPI_H)
void
MPI_resultset1dfm_assemble_results(resultset1dfm_t *r,
				   int mpisize,
				   int mpirank,
				   int root,
				   MPI_Comm comm)
{
  int ivsize;
  int *iv;

  int vsize;
  double *v;

  int cred_samples;

  int i;
  int j;

  double denom;

  vsize = r->xsamples;
  v = rjmcmc_create_array_1d(vsize);
  if (v == NULL) {
    rjmcmc_error("MPI_resultset1dfm_assemble_results: "
		 "failed to create temporary array\n");
    return;
  }

  ivsize = r->nprocesses;
  if (r->histogram && ivsize < r->ysamples) {
    ivsize = r->ysamples;
  }
  iv = rjmcmc_create_int_array_1d(ivsize);
  if (iv == NULL) {
    rjmcmc_error("MPI_resultset1dfm_assemble_results: "
		 "failed to create temporary int array\n");
    return;
  }

  /*
   * Aggregate the total acceptance/proposal counts
   */
  MPI_Reduce(r->propose, iv, r->nprocesses, MPI_INT, MPI_SUM, root, comm);
  if (mpirank == root) {
    for (j = 0; j < r->nprocesses; j ++) {
      r->propose[j] = iv[j];
    }
  }

  MPI_Reduce(r->accept, iv, r->nprocesses, MPI_INT, MPI_SUM, root, comm);
  if (mpirank == root) {
    for (j = 0; j < r->nprocesses; j ++) {
      r->accept[j] = iv[j];
    }
  }

  /*
   * Aggregate the means. First compute locally then distribute.
   */
  if (r->local_mean) {

    denom = r->total - r->burnin;

    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) { 
	r->local_mean[i][j] /= denom;
      }

      MPI_Reduce(r->local_mean[i], v, r->xsamples, MPI_DOUBLE, MPI_SUM, root, comm);

      if (mpirank == root) {
	for (j = 0; j < r->xsamples; j ++) {
	  r->local_mean[i][j] = v[j]/(double)mpisize;
	}
      }
    }
  }

  /*
   * Aggregate the histograms
   */
  if (r->histogram) {

    cred_samples = ((double)(r->total - r->burnin) * (1.0 - r->credible_interval)/2.0) * mpisize;

    for (i = 0; i < r->nlocalparameters; i ++) {
      
      for (j = 0; j < r->xsamples; j ++) {

	MPI_Reduce(r->histogram[i][j], iv, r->ysamples, 
		   MPI_INT, MPI_SUM, root, comm);

	

	if (mpirank == root) {

	  if (r->local_median) {
	    r->local_median[i][j] = 
	      rjmcmc_median_from_histogram(iv,
					   r->local_parameters[i].fmin,
					   r->local_parameters[i].fmax,
					   r->ysamples);
	  }

	  if (r->local_mode) {
	    r->local_mode[i][j] = 
	      rjmcmc_mode_from_histogram(iv,
					 r->local_parameters[i].fmin,
					 r->local_parameters[i].fmax,
					 r->ysamples);
	  }

	  if (r->local_cred_min &&
	      r->local_cred_max) {

	    r->local_cred_min[i][j] = 
	      rjmcmc_head_from_histogram(iv,
					 r->local_parameters[i].fmin,
					 r->local_parameters[i].fmax,
					 r->ysamples,
					 cred_samples);
	    
	    r->local_cred_max[i][j] = 
	      rjmcmc_tail_from_histogram(iv,
					 r->local_parameters[i].fmin,
					 r->local_parameters[i].fmax,
					 r->ysamples,
					 cred_samples);

	  }
	
	}
      }

    }
  }      

  rjmcmc_destroy_int_array_1d(iv);
  rjmcmc_destroy_array_1d(v);
}
#endif /* HAVE_MPI_H */
/*
 * Getting results
 */

int 
resultset1dfm_get_max_partitions(resultset1dfm_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->maxpartitions;
}

int
resultset1dfm_get_total(resultset1dfm_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->total;
}

int
resultset1dfm_get_xsamples(resultset1dfm_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->xsamples;
}

int
resultset1dfm_get_nparameters(resultset1dfm_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->nprocesses;
}

const int *
resultset1dfm_get_propose(resultset1dfm_t *r,
			int *nprocesses)
{
  if (r == NULL) {
    return NULL;
  }

  if (nprocesses != NULL) {
    *nprocesses = r->nprocesses;
  }

  return r->propose;
}

int
resultset1dfm_get_propose_f(resultset1dfm_t *r,
			    int maxsize,
			    int *propose)
{
  int n;
  int i;

  if (r == NULL) {
    return -1;
  }

  n = resultset1dfm_get_nparameters(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    propose[i] = r->propose[i];
  }
  
  return n;
}

const int *
resultset1dfm_get_accept(resultset1dfm_t *r,
		       int *nprocesses)
{
  if (r == NULL) {
    return NULL;
  }

  if (nprocesses != NULL) {
    *nprocesses = r->nprocesses;
  }

  return r->accept;
}

int
resultset1dfm_get_accept_f(resultset1dfm_t *r,
			    int maxsize,
			    int *accept)
{
  int n;
  int i;

  if (r == NULL) {
    return -1;
  }

  n = resultset1dfm_get_nparameters(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    accept[i] = r->accept[i];
  }
  
  return n;
}

const double *
resultset1dfm_get_misfit(resultset1dfm_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->misfit;
}

int
resultset1dfm_get_misfit_f(resultset1dfm_t *r,
			   int maxsize,
			   double *misfit)
{
  int i;
  int n;

  if (r == NULL) {
    return -1;
  }

  n = maxsize;
  if (n > r->total) {
    n = r->total;
  }

  for (i = 0; i < n; i ++) {
    misfit[i] = r->misfit[i];
  }

  return n;
}


const int *
resultset1dfm_get_partitions(resultset1dfm_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->partitions;
}

int
resultset1dfm_get_partitions_f(resultset1dfm_t *r,
			       int maxsize,
			       int *partitions)
{
  int i;
  int n;

  n = maxsize;
  if (n > r->total) {
    n = r->total;
  }

  for (i = 0; i < n; i ++) {
    partitions[i] = r->partitions[i];
  }

  return n;
}


const int *
resultset1dfm_get_partition_x_histogram(resultset1dfm_t *r)
{
  return r->partition_x_hist;
}

int
resultset1dfm_get_partition_x_histogram_f(resultset1dfm_t *r,
					  int maxsize,
					  int *histogram)
{
  int i;
  int n;

  n = maxsize;
  if (n > r->xsamples) {
    n = r->xsamples;
  }

  for (i = 0; i < n; i ++) {
    histogram[i] = r->partition_x_hist[i];
  }

  return n;
}


const double *
resultset1dfm_get_global_parameter(resultset1dfm_t *r, int gi)
{
  
  if (gi >= 0 && gi < r->nglobalparameters) {
    return r->global[gi];
  }

  return NULL;
}

int
resultset1dfm_get_global_parameter_f(resultset1dfm_t *r, 
				     int gi,
				     int maxsize,
				     double *global)
{
  int i;
  int n;
  
  if (r == NULL) {
    return -1;
  }

  if (r->global == NULL) {
    return -1;
  }

  if (gi < 0 || gi >= r->nglobalparameters) {
    return -1;
  }

  n = maxsize;
  if (n > r->total) {
    n = r->total;
  }

  for (i = 0; i < n; i ++) {
    global[i] = r->global[gi][i];
  }

  return n;
}

const double *
resultset1dfm_get_local_parameter_mean(resultset1dfm_t *r, int li)
{
  if (li >= 0 && li < r->nlocalparameters) {
    return r->local_mean[li];
  }

  return NULL;
}

int resultset1dfm_get_local_parameter_mean_f(resultset1dfm_t *r,
					     int li,
					     int maxsize,
					     double *mean)
{
  int i;
  int n;

  n = maxsize;
  if (n > r->xsamples) {
    n = r->xsamples;
  }

  for (i = 0; i < n; i ++) {
    mean[i] = r->local_mean[li][i];
  }

  return n;
}

const double *
resultset1dfm_get_local_parameter_median(resultset1dfm_t *r, int li)
{
  if (r->local_median != NULL) {
    return r->local_median[li];
  }

  return NULL;
}

int resultset1dfm_get_local_parameter_median_f(resultset1dfm_t *r,
					       int li,
					       int maxsize,
					       double *median)
{
  int i;
  int n;

  if (r->local_median == NULL) {
    return 0;
  }

  n = maxsize;
  if (n > r->xsamples) {
    n = r->xsamples;
  }

  for (i = 0; i < n; i ++) {
    median[i] = r->local_mean[li][i];
  }

  return n;
}


const double *
resultset1dfm_get_local_parameter_mode(resultset1dfm_t *r, int li)
{
  if (r->local_mode != NULL) {
    return r->local_mode[li];
  }

  return NULL;
}

int resultset1dfm_get_local_parameter_mode_f(resultset1dfm_t *r,
					     int li,
					     int maxsize,
					     double *mode)
{
  int i;
  int n;

  if (r->local_mode == NULL) {
    return 0;
  }

  n = maxsize;
  if (n > r->xsamples) {
    n = r->xsamples;
  }

  for (i = 0; i < n; i ++) {
    mode[i] = r->local_mode[li][i];
  }

  return n;
}

const double *
resultset1dfm_get_local_parameter_credible_min(resultset1dfm_t *r, int li)
{
  if (r->local_cred_min != NULL) {
    return r->local_cred_min[li];
  }

  return NULL;
}

int resultset1dfm_get_local_parameter_credible_min_f(resultset1dfm_t *r,
						     int li,
						     int maxsize,
						     double *credible_min)
{
  int i;
  int n;

  if (r->local_cred_min == NULL) {
    return 0;
  }

  n = maxsize;
  if (n > r->xsamples) {
    n = r->xsamples;
  }

  for (i = 0; i < n; i ++) {
    credible_min[i] = r->local_cred_min[li][i];
  }

  return n;
}

const double *
resultset1dfm_get_local_parameter_credible_max(resultset1dfm_t *r, int li)
{
  if (r->local_cred_max != NULL) {
    return r->local_cred_max[li];
  }
  
  return NULL;
}

int resultset1dfm_get_local_parameter_credible_max_f(resultset1dfm_t *r,
						     int li,
						     int maxsize,
						     double *credible_max)
{
  int i;
  int n;

  if (r->local_cred_max == NULL) {
    return 0;
  }

  n = maxsize;
  if (n > r->xsamples) {
    n = r->xsamples;
  }

  for (i = 0; i < n; i ++) {
    credible_max[i] = r->local_cred_max[li][i];
  }

  return n;
}

const int **
resultset1dfm_get_local_parameter_histogram(resultset1dfm_t *r, int li)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset1dfm_get_local_parameter_histogram: "
		      "null results.");
  RJMCMC_NULLCHECKPTR(r->histogram,
		      "resultset1dfm_get_local_parameter_histogram: "
		      "NULL histogram.\n");
  RJMCMC_INDEXCHECKPTR(li, r->nlocalparameters, 
		       "resultset1dfm_get_local_parameter_histogram: "
		       "invalid index.");

  return (const int **)(r->histogram[li]);
}

int
resultset1dfm_get_local_parameter_histogram_f(resultset1dfm_t *r, 
					      int li,
					      int xsamples,
					      int ysamples,
					      int *histogram)
{
  const int **lhist;
  int x;
  int y;

  lhist = resultset1dfm_get_local_parameter_histogram(r, li);
  if (lhist == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset1dfm_get_local_parameter_histogram_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_local_parameter_histogram_f: "
			   "invalid ysamples\n");
  /*
   * This is a test to ensure that the fortran array is not transposed.
   * Run with a fortran forward model with different xsamples and ysamples
   * and ensure that the 2d histogram has a circle at its centre.
   */

  //#define CIRCLE_TEST
#ifdef CIRCLE_TEST
  
  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {
      int r2 = 
	(y - ysamples/2)*(y - ysamples/2) + 
	(x - xsamples/2)*(x - xsamples/2);

      if (r2 > 16) {
	histogram[y * xsamples + x] = 0;
      } else {
	histogram[y * xsamples + x] = 1;
      }
	
    }
  }

#else
  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      histogram[y * xsamples + x] = lhist[x][y];

    }
  }
#endif

  return 0;
}

const double *
resultset1dfm_get_hierarchical(resultset1dfm_t *r,
			       int si)
{
  if (si >= 0 && si < r->nhierarchical) {
    return r->hierarchical[si];
  } 

  return NULL;
}

int
resultset1dfm_get_hierarchical_f(resultset1dfm_t *r, 
				 int hi,
				 int maxsize,
				 double *hierarchical)
{
  int i;
  int n;
  
  if (r == NULL) {
    return -1;
  }

  if (r->hierarchical == NULL) {
    return -1;
  }

  if (hi < 0 || hi >= r->nhierarchical) {
    return -1;
  }

  n = maxsize;
  if (n > r->total) {
    n = r->total;
  }

  for (i = 0; i < n; i ++) {
    hierarchical[i] = r->hierarchical[hi][i];
  }

  return n;
}

void
resultset1dfm_fill_xcoord_vector(resultset1dfm_t *r,
				 double *x,
				 int *l)
{
  int i;
  int e;

  e = *l;
  if (e > r->xsamples) {
    e = r->xsamples;
    *l = e;
  }

  rjmcmc_fill_coord_vector(r->xmin,
			   r->xmax,
			   e,
			   x);
}

int
resultset1dfm_get_xcoord_vector_f(resultset1dfm_t *r,
				  int maxsize,
				  double *x)
{
  int s;

  s = maxsize;
  resultset1dfm_fill_xcoord_vector(r, x, &s);

  return s;
}

void 
resultset1dfm_fill_ycoord_vector(resultset1dfm_t *r,
				 int li,
				 double *y,
				 int *l)
{
  int i;
  int e;

  
  e = *l;
  if (e > r->xsamples) {
    e = r->xsamples;
    *l = e;
  }

  
  rjmcmc_fill_coord_vector(r->local_parameters[li].fmin,
			   r->local_parameters[li].fmax,
			   r->ysamples,
			   y);

}

int
resultset1dfm_get_ycoord_vector_f(resultset1dfm_t *r,
				  int li,
				  int maxsize,
				  double *y)
{
  int s;

  s = maxsize;
  resultset1dfm_fill_ycoord_vector(r, li, y, &s);

  return s;
}
