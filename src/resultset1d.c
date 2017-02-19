

#include <stdio.h>
#include <stdlib.h>

#include <rjmcmc/resultset1d.h>
#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_defines.h>

struct _resultset1d {

  int results;

  int burnin;
  int total;

  int xsamples;
  int ysamples;

  int nprocesses;

  int maxpartitions;
  int maxorder;

  double xmin;
  double xmax;
  double ymin;
  double ymax;

  double gradmin;
  double gradmax;

  int *propose;
  int *accept;

  /* 
   * Misfit [total]
   */
  double *misfit;

  /*
   * Lambda [total]
   */
  double *lambda;

  /*
   * Order [maxorder + 1];
   */
  int *order;

  /*
   * N Partitions [total];
   */
  int *partitions;

  /*
   * Partition x histogram [xsamples]
   */
  int *partition_x_hist;

  /*
   * Mean [xsamples]
   */
  double *mean;

  /*
   * Hist [xsamples][ysamples]
   */
  int **hist;

  /*
   * Mode [xsamples]
   */
  double *mode;

  /*
   * Median [xsamples]
   */
  double *median;

  /*
   * Credible [xsamples] x2
   */
  double conf_interval;
  double *conf_min;
  double *conf_max;

  /*
   * Gradient [xsamples]
   */
  double *gradient;

  /*
   * Gradient Hist [xsamples][ysamples]
   */
  int **gradient_hist;

  double *gradient_conf_min;
  double *gradient_conf_max;

};

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
		   int results)
{
  resultset1d_t *r;

  r = (resultset1d_t *)malloc(sizeof(resultset1d_t));
  if (r == NULL) {
    return NULL;
  }

  r->results = results;
  r->burnin = burnin;
  r->total = total;
  
  r->xsamples = xsamples;
  r->ysamples = ysamples;

  r->nprocesses = nprocesses;

  r->maxpartitions = maxpartitions;
  r->maxorder = maxorder;

  r->xmin = xmin;
  r->xmax = xmax;
  r->ymin = ymin;
  r->ymax = ymax;

  r->gradmax = (r->ymax - r->ymin)/
    ((r->xmax - r->xmin) / (double)(maxpartitions + 1));
  r->gradmin = -(r->gradmax);

  r->propose = rjmcmc_create_int_array_1d(nprocesses);
  if (r->propose == NULL) {
    return NULL;
  }

  r->accept = rjmcmc_create_int_array_1d(nprocesses);
  if (r->accept == NULL) {
    return NULL;
  }

  if (total > 0) {
    r->misfit = rjmcmc_create_array_1d(total);
    if (r->misfit == NULL) {
      return NULL;
    }

    r->lambda = rjmcmc_create_array_1d(total);
    if (r->lambda == NULL) {
      return NULL;
    }
  } else {
    r->misfit = NULL;
    r->lambda = NULL;
  }

  r->order = NULL;
  if (r->maxorder >= 0) {
    r->order = rjmcmc_create_int_array_1d(r->maxorder + 1);
    if (r->order == NULL) {
      return NULL;
    }
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

  r->mean = NULL;
  if (results & RESULTSET1D_MEAN) {
    r->mean = rjmcmc_create_array_1d(xsamples);
    if (r->mean == NULL) {
      return NULL;
    }
  }
   
  r->hist = NULL;
  if (results & (RESULTSET1D_MEDIAN | 
		 RESULTSET1D_MODE | 
		 RESULTSET1D_CREDIBLE)) {
    r->hist = rjmcmc_create_int_array_2d(xsamples, ysamples);
    if (r->hist == NULL) {
      return NULL;
    }
  }

  r->mode = NULL;
  if (results & RESULTSET1D_MODE) {
    r->mode = rjmcmc_create_array_1d(xsamples);
    if (r->mode == NULL) {
      return NULL;
    }
  }

  r->median = NULL;
  if (results & RESULTSET1D_MEDIAN) {
    r->median = rjmcmc_create_array_1d(xsamples);
    if (r->median == NULL) {
      return NULL;
    }
  }

  r->conf_interval = credible_interval;
  r->conf_min = NULL;
  r->conf_max = NULL;
  if (results & RESULTSET1D_CREDIBLE) {
    r->conf_min = rjmcmc_create_array_1d(xsamples);
    r->conf_max = rjmcmc_create_array_1d(xsamples);
    if (r->conf_min == NULL ||
	r->conf_max == NULL) {
      return NULL;
    }
  }

  r->gradient = NULL;
  if (results & RESULTSET1D_GRADIENT) {
    r->gradient = rjmcmc_create_array_1d(xsamples);
    if (r->gradient == NULL) {
      return NULL;
    }
  }

  r->gradient_hist = NULL;
  r->gradient_conf_min = NULL;
  r->gradient_conf_max = NULL;
  if ((results & RESULTSET1D_GRADIENT) &&
      (results & RESULTSET1D_CREDIBLE)) {
    r->gradient_hist = rjmcmc_create_int_array_2d(xsamples, ysamples);
    if (r->gradient_hist == NULL) {
      return NULL;
    }

    r->gradient_conf_min = rjmcmc_create_array_1d(xsamples);
    r->gradient_conf_max = rjmcmc_create_array_1d(xsamples);
    if (r->gradient_conf_min == NULL ||
	r->gradient_conf_max == NULL) {
      return NULL;
    }
  }

  return r;
}

void 
resultset1d_destroy(resultset1d_t *r)
{
  if (r != NULL) {

    rjmcmc_destroy_int_array_1d(r->propose);
    rjmcmc_destroy_int_array_1d(r->accept);
    
    rjmcmc_destroy_array_1d(r->misfit);
    rjmcmc_destroy_array_1d(r->lambda);
    rjmcmc_destroy_int_array_1d(r->order);
    rjmcmc_destroy_int_array_1d(r->partitions);
    rjmcmc_destroy_int_array_1d(r->partition_x_hist);

    rjmcmc_destroy_array_1d(r->mean);

    rjmcmc_destroy_int_array_2d(r->xsamples, r->hist);

    rjmcmc_destroy_array_1d(r->mode);
    rjmcmc_destroy_array_1d(r->median);
    rjmcmc_destroy_array_1d(r->conf_min);
    rjmcmc_destroy_array_1d(r->conf_max);

    rjmcmc_destroy_int_array_2d(r->xsamples, r->gradient_hist);

    rjmcmc_destroy_array_1d(r->gradient);
    rjmcmc_destroy_array_1d(r->gradient_conf_min);
    rjmcmc_destroy_array_1d(r->gradient_conf_max);

    free(r);
  }
}

void
resultset1d_propose(resultset1d_t *r,
		    int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses, 
			"resultset1d_propose: invalid index\n");
  r->propose[p] ++;
}

void
resultset1d_accept(resultset1d_t *r,
		   int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses, 
			"resultset1d_accept: invalid index\n");
  r->accept[p] ++;
}

void
resultset1d_sample(resultset1d_t *r,
		   int i,
		   const double *v)
{
  int j;

  RJMCMC_INDEXCHECKVOID(i,
			r->total,
			"resulset1d_sample: invalid index\n");

  if (i >= r->burnin) {
    
    /*
     * Sample for mean
     */
    if (r->mean) {
      for (j = 0; j < r->xsamples; j ++) {
	r->mean[j] += v[j];
      }
    }

    if (r->hist != NULL) {
      /*
       * Sample for mode, median, credible
       */
      for (j = 0; j < r->xsamples; j ++) {
	r->hist[j][rjmcmc_map_to_index(v[j],
				       r->ymin,
				       r->ymax,
				       r->ysamples)] ++;
      }
    }
  }
}

void
resultset1d_sample_gradient(resultset1d_t *r,
			    int i,
			    const double *dv)
{
  int j;

  RJMCMC_INDEXCHECKVOID(i,
			r->total,
			"resulset1d_sample_gradient: invalid index\n");

  if (i >= r->burnin) {
    
    /*
     * Sample for mean
     */
    if (r->gradient) {
      for (j = 0; j < r->xsamples; j ++) {
	r->gradient[j] += dv[j];
      }
    }

    if (r->gradient_hist != NULL) {
      /*
       * Sample for mode, median, credible
       */
      for (j = 0; j < r->xsamples; j ++) {
	r->gradient_hist[j][rjmcmc_map_to_index(dv[j],
						r->gradmin,
						r->gradmax,
						r->ysamples)] ++;
      }
    }
  }
}

void 
resultset1d_sample_misfit(resultset1d_t *r,
			  int i,
			  double misfit)
{
  RJMCMC_INDEXCHECKVOID(i,
			r->total,
			"resulset1d_sample_misfit: invalid index\n");

  r->misfit[i] = misfit;
}

void 
resultset1d_sample_order(resultset1d_t *r,
			 int i,
			 int order)
{
  RJMCMC_INDEXCHECKVOID(order,
			r->maxorder + 1,
			"resultset1d_sample_order: invalid order\n");
  r->order[order] ++;
}

void 
resultset1d_sample_lambda(resultset1d_t *r,
			 int i,
			 double lambda)
{
  RJMCMC_INDEXCHECKVOID(i,
			r->total,
			"resulset1d_sample_lambda: invalid index\n");

  r->lambda[i] = lambda;
}

void
resultset1d_sample_npartitions(resultset1d_t *r,
			       int i,
			       int npartitions)
{
  RJMCMC_INDEXCHECKVOID(i,
			r->total,
			"resulset1d_sample_npartitions: invalid index\n");

  r->partitions[i] = npartitions;
}

void
resultset1d_sample_partition_x(resultset1d_t *r,
			       int i,
			       double x)
{
  int j;

  if (x > r->xmin && x < r->xmax) {
    j = rjmcmc_map_to_index(x, r->xmin, r->xmax, r->xsamples);
    r->partition_x_hist[j] ++;
  }
}

void
resultset1d_assemble_results(resultset1d_t *r)
{
  int j;
  double denom;
  int conf_samples;

  denom = r->total - r->burnin;

  if (r->mean) {
    for (j = 0; j < r->xsamples; j ++) { 
      r->mean[j] /= denom;
    }
  }

  if (r->gradient) {
    for (j = 0; j < r->xsamples; j ++) { 
      r->gradient[j] /= denom;
    }
  }
    

  if (r->median) {
    for (j = 0; j < r->xsamples; j ++) {
      r->median[j] = rjmcmc_median_from_histogram(r->hist[j],
						  r->ymin,
						  r->ymax,
						  r->ysamples);
    }
  }

  if (r->mode) {
    for (j = 0; j < r->xsamples; j ++) {
      r->mode[j] = rjmcmc_mode_from_histogram(r->hist[j],
					      r->ymin,
					      r->ymax,
					      r->ysamples);
    }
  }

  if (r->conf_min && r->conf_max) {

    conf_samples = (double)(r->total - r->burnin) * (1.0 - r->conf_interval)/2.0;

    for (j = 0; j < r->xsamples; j ++) {

      r->conf_min[j] = rjmcmc_head_from_histogram(r->hist[j],
						  r->ymin,
						  r->ymax,
						  r->ysamples,
						  conf_samples);
      r->conf_max[j] = rjmcmc_tail_from_histogram(r->hist[j],
						  r->ymin,
						  r->ymax,
						  r->ysamples,
						  conf_samples);
    }
  }

  if (r->gradient_conf_min && r->gradient_conf_max) {
    conf_samples = (double)(r->total - r->burnin) * (1.0 - r->conf_interval)/2.0;

    for (j = 0; j < r->xsamples; j ++) {
      r->gradient_conf_min[j] = rjmcmc_head_from_histogram(r->gradient_hist[j],
							   r->gradmin,
							   r->gradmax,
							   r->ysamples,
							   conf_samples);
      r->gradient_conf_max[j] = rjmcmc_tail_from_histogram(r->gradient_hist[j],
							   r->gradmin,
							   r->gradmax,
							   r->ysamples,
							   conf_samples);
    }
  }    
}

#if defined(HAVE_MPI_H)
void
MPI_resultset1d_assemble_results(resultset1d_t *r,
				 int mpisize,
				 int mpirank,
				 int root,
				 MPI_Comm comm)
{
  int j;
  double denom;
  int conf_samples;

  double *v;
  int *iv;

  denom = (r->total - r->burnin) * (double)mpisize;

  v = rjmcmc_create_array_1d(r->xsamples);
  if (v == NULL) {
    rjmcmc_error("MPI_resultset1d_assemble_results: "
		 "failed to create temporary array\n");
    return;
  }

  j = r->ysamples;
  if (j < r->nprocesses) {
    j = r->nprocesses;
  }
  iv = rjmcmc_create_int_array_1d(j);
  if (iv == NULL) {
    rjmcmc_error("MPI_resultset1d_assemble_results: "
		 "failed to create temporary histogram array\n");
    return;
  }

  if (r->mean) {
    MPI_Reduce(r->mean, v, r->xsamples, MPI_DOUBLE, MPI_SUM, root, comm);
    if (mpirank == root) {
      for (j = 0; j < r->xsamples; j ++) { 
	r->mean[j] = v[j]/denom;
      }
    }
  }

  if (r->gradient) {
    MPI_Reduce(r->gradient, v, r->xsamples, MPI_DOUBLE, MPI_SUM, root, comm);
    if (mpirank == root) {
      for (j = 0; j < r->xsamples; j ++) { 
	r->gradient[j] = v[j]/denom;
      }
    }
  }
    
  if (r->hist) {

    for (j = 0; j < r->xsamples; j ++) {

      MPI_Reduce(r->hist[j], iv, r->ysamples, MPI_INT, MPI_SUM, root, comm);

      if (mpirank == root) {

	if (r->median) {
	  r->median[j] = rjmcmc_median_from_histogram(iv,
						      r->ymin,
						      r->ymax,
						      r->ysamples);
	}

	if (r->mode) {
	  r->mode[j] = rjmcmc_mode_from_histogram(iv,
						  r->ymin,
						  r->ymax,
						  r->ysamples);
	}

	if (r->conf_min && r->conf_max) {
	  conf_samples = (double)(mpisize * (r->total - r->burnin)) * (1.0 - r->conf_interval)/2.0;

	  r->conf_min[j] = rjmcmc_head_from_histogram(iv,
						      r->ymin,
						      r->ymax,
						      r->ysamples,
						      conf_samples);
	  r->conf_max[j] = rjmcmc_tail_from_histogram(iv,
						      r->ymin,
						      r->ymax,
						      r->ysamples,
						      conf_samples);
	}
      }
    }
  }

  if (r->gradient_hist) {
    conf_samples = (double)(mpisize * (r->total - r->burnin)) * (1.0 - r->conf_interval)/2.0;

    for (j = 0; j < r->xsamples; j ++) {

      MPI_Reduce(r->gradient_hist[j], iv, r->ysamples, 
		 MPI_INT, MPI_SUM, root, comm);

      

      if (mpirank == root) {
	if (r->gradient_conf_min && r->gradient_conf_max) {

	  r->gradient_conf_min[j] = rjmcmc_head_from_histogram(iv,
							       r->gradmin,
							       r->gradmax,
							       r->ysamples,
							       conf_samples);
	  r->gradient_conf_max[j] = rjmcmc_tail_from_histogram(iv,
							       r->gradmin,
							       r->gradmax,
							       r->ysamples,
							       conf_samples);
	}
      }
    }
  }    

  /*
   * In the MPI version we also need to aggregate the acceptance and
   * proposal counts. 
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

}
#endif /* HAVE_MPI_H */

/*
 * Getting results
 */

const int *
resultset1d_get_propose(resultset1d_t *r,
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

const int *
resultset1d_get_accept(resultset1d_t *r,
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

const double *
resultset1d_get_misfit(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->misfit;
}

const int *
resultset1d_get_order(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->order;
}

const double *
resultset1d_get_lambda(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->lambda;
}

const int *
resultset1d_get_partitions(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->partitions;
}

const int *
resultset1d_get_partition_x_histogram(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->partition_x_hist;
}

const double *
resultset1d_get_mean(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->mean;
}

const double *
resultset1d_get_median(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->median;
}

const double *
resultset1d_get_mode(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->mode;
}

const double *
resultset1d_get_credible_min(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->conf_min;
}
		     
const double *
resultset1d_get_credible_max(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->conf_max;
}

const int **
resultset1d_get_histogram(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

return (const int **)r->hist;
}

const double *
resultset1d_get_gradient(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->gradient;
}

const double *
resultset1d_get_gradient_credible_min(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->gradient_conf_min;
}

const double *
resultset1d_get_gradient_credible_max(resultset1d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->gradient_conf_max;
}

int
resultset1d_get_total(resultset1d_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->total;
}

int 
resultset1d_get_max_partitions(resultset1d_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->total;
}

int 
resultset1d_get_xsamples(resultset1d_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->xsamples;
}

int
resultset1d_get_ysamples(resultset1d_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->ysamples;
}

int 
resultset1d_get_max_order(resultset1d_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->maxorder;
}

void
resultset1d_fill_xcoord_vector(resultset1d_t *r,
			       double *x)
{
  rjmcmc_fill_coord_vector(r->xmin,
			   r->xmax,
			   r->xsamples,
			   x);
}

void
resultset1d_fill_ycoord_vector(resultset1d_t *r,
			       double *y)
{
  rjmcmc_fill_coord_vector(r->ymin,
			   r->ymax,
			   r->ysamples,
			   y);
}
