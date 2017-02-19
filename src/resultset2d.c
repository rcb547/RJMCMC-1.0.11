
#include <stdio.h>
#include <stdlib.h>
#include "rjmcmc/resultset2d.h"

#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_defines.h"

struct _resultset2d {

  int results;

  int burnin;
  int total;

  int xsamples;
  int ysamples;
  int zsamples;

  int nprocesses;

  int max_partitions;

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;

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
   * N Partitions [total]
   */
  int *partitions;

  /*
   * Mean [xsamples][ysamples]
   */
  double **mean;

  /*
   * Centres [xsamples][ysamples]
   */
  int **centres;

  /*
   * Hist [xsamples][ysamples][zsamples]
   */
  int ***hist;

  /*
   * Mode [xsamples][ysamples]
   */
  double **mode;

  /*
   * Median [xsamples][ysamples]
   */
  double **median;

  /*
   * Credible [xsamples][ysamples] x2
   */
  double conf_interval;
  double **conf_min;
  double **conf_max;

};

resultset2d_t *
resultset2d_create(int burnin,
		   int total,
		   int xsamples,
		   int ysamples,
		   int zsamples,
		   int nprocesses,
		   int maxpartitions,
		   double xmin,
		   double xmax,
		   double ymin,
		   double ymax,
		   double zmin,
		   double zmax,
		   double credible_interval,
		   int results)
{
  resultset2d_t *r;
  int i;

  r = (resultset2d_t *)malloc(sizeof(resultset2d_t));
  if (r == NULL) {
    return NULL;
  }

  r->results = results;
  r->burnin = burnin;
  r->total = total;

  r->xsamples = xsamples;
  r->ysamples = ysamples;
  r->zsamples = zsamples;

  r->nprocesses = nprocesses;

  r->max_partitions = maxpartitions;

  r->xmin = xmin;
  r->xmax = xmax;
  r->ymin = ymin;
  r->ymax = ymax;
  r->zmin = zmin;
  r->zmax = zmax;

  r->propose = rjmcmc_create_int_array_1d(nprocesses);
  if (r->propose == NULL) {
    return NULL;
  }
  
  r->accept = rjmcmc_create_int_array_1d(nprocesses);
  if (r->accept == NULL) {
    return NULL;
  }
  
  r->misfit = rjmcmc_create_array_1d(total);
  if (r->misfit == NULL) {
    return NULL;
  }

  r->lambda = rjmcmc_create_array_1d(total);
  if (r->lambda == NULL) {
    return NULL;
  }

  r->partitions = rjmcmc_create_int_array_1d(total);
  if (r->partitions == NULL) {
    return NULL;
  }

  r->mean = NULL;
  if (results & RESULTSET2D_MEAN) {
    r->mean = rjmcmc_create_array_2d(xsamples, ysamples);
    if (r->mean == NULL) {
      return NULL;
    }
  }
   
  r->centres = rjmcmc_create_int_array_2d(xsamples, ysamples);
  if (r->centres == NULL) {
    return NULL;
  }
    
  r->hist = NULL;
  if (results & (RESULTSET2D_MEDIAN | 
		 RESULTSET2D_MODE | 
		 RESULTSET2D_CREDIBLE)) {
    r->hist = rjmcmc_create_int_array_3d(xsamples, ysamples, zsamples);
    if (r->hist == NULL) {
      return NULL;
    }
  }

  r->mode = NULL;
  if (results & RESULTSET2D_MODE) {
    r->mode = rjmcmc_create_array_2d(xsamples, ysamples);
    if (r->mode == NULL) {
      return NULL;
    }
  }

  r->median = NULL;
  if (results & RESULTSET2D_MEDIAN) {
    r->median = rjmcmc_create_array_2d(xsamples, ysamples);
    if (r->median == NULL) {
      return NULL;
    }
  }

  r->conf_interval = credible_interval;
  r->conf_min = NULL;
  r->conf_max = NULL;
  if (results & RESULTSET2D_CREDIBLE) {
    r->conf_min = rjmcmc_create_array_2d(xsamples, ysamples);
    if (r->conf_min == NULL) {
      return NULL;
    }
    r->conf_max = rjmcmc_create_array_2d(xsamples, ysamples);
    if (r->conf_max == NULL) {
      return NULL;
    }
  }
    

  return r;
}

void 
resultset2d_destroy(resultset2d_t *r)
{
  if (r != NULL) {

    rjmcmc_destroy_int_array_1d(r->propose);
    rjmcmc_destroy_int_array_1d(r->accept);

    rjmcmc_destroy_array_1d(r->misfit);
    rjmcmc_destroy_array_1d(r->lambda);
    rjmcmc_destroy_int_array_1d(r->partitions);

    rjmcmc_destroy_array_2d(r->xsamples, r->mean);

    rjmcmc_destroy_int_array_2d(r->xsamples, r->centres);

    rjmcmc_destroy_int_array_3d(r->xsamples, r->ysamples, r->hist);

    rjmcmc_destroy_array_2d(r->xsamples, r->mode);
    rjmcmc_destroy_array_2d(r->xsamples, r->median);
    rjmcmc_destroy_array_2d(r->xsamples, r->conf_min);
    rjmcmc_destroy_array_2d(r->xsamples, r->conf_max);

    free(r);
  }
   
}

void resultset2d_propose(resultset2d_t *r,
			 int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses,
			"resultset2d_propose: invalid process\n");

  r->propose[p] ++;
}

void
resultset2d_accept(resultset2d_t *r,
		   int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses,
			"resultset2d_accept: invalid process\n");

  r->accept[p] ++;
}

void
resultset2d_sample(resultset2d_t *r,
		   int i,
		   const double **v)
{
  int j;
  int k;

  RJMCMC_INDEXCHECKVOID(i, 
			r->total,
			"resultset2d_sample: invalid index\n");

  if (i >= r->burnin) {
    
    /*
     * Sample for mean
     */
    if (r->mean) {
      for (j = 0; j < r->xsamples; j ++) {
	for (k = 0; k < r->ysamples; k ++) {
	  r->mean[j][k] += v[j][k];
	}
      }
    }

    if (r->hist != NULL) {
      /*
       * Sample for mode, median, credible
       */
      for (j = 0; j < r->xsamples; j ++) {
	for (k = 0; k < r->ysamples; k ++) {
	  r->hist[j][k][rjmcmc_map_to_index(v[j][k],
					    r->zmin,
					    r->zmax,
					    r->zsamples)] ++;
	}
      }
    }
  }  
}

void 
resultset2d_sample_misfit(resultset2d_t *r,
			  int i,
			  double misfit)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total,
			"resultset2d_sample_misfit: invalid index\n");
  r->misfit[i] = misfit;
}

void
resultset2d_sample_centre(resultset2d_t *r,
			  double x,
			  double y)
{
  int ix;
  int iy;

  ix = rjmcmc_map_to_index(x, r->xmin, r->xmax, r->xsamples);
  iy = rjmcmc_map_to_index(y, r->ymin, r->ymax, r->ysamples);
  r->centres[ix][iy] ++;
}

void 
resultset2d_sample_lambda(resultset2d_t *r,
			 int i,
			 double lambda)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total,
			"resultset2d_sample_lambda: invalid index\n");
  r->lambda[i] = lambda;
}

void
resultset2d_sample_npartitions(resultset2d_t *r,
			       int i,
			       int npartitions)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total,
			"resultset2d_sample_npartitions: invalid index\n");
  r->partitions[i] = npartitions;
}

void
resultset2d_assemble_results(resultset2d_t *r)
{
  int j;
  int k;
  double denom;
  int conf_samples;

  denom = r->total - r->burnin;

  if (r->mean) {
    for (j = 0; j < r->xsamples; j ++) {
      for (k = 0; k < r->ysamples; k ++) {
	r->mean[j][k] /= denom;
      }
    }
  }

  if (r->median) {
    for (j = 0; j < r->xsamples; j ++) {
      for (k = 0; k < r->ysamples; k ++) {
	r->median[j][k] = rjmcmc_median_from_histogram(r->hist[j][k],
						       r->zmin,
						       r->zmax,
						       r->zsamples);
      }
    }
  }

  if (r->mode) {
    for (j = 0; j < r->xsamples; j ++) {
      for (k = 0; k < r->ysamples; k ++) {
	r->mode[j][k] = rjmcmc_mode_from_histogram(r->hist[j][k],
						   r->zmin,
						   r->zmax,
						   r->zsamples);
      }
    }
  }
  
  if (r->conf_min && r->conf_max) {

    conf_samples = (double)(r->total - r->burnin) * (1.0 - r->conf_interval)/2.0;

    for (j = 0; j < r->xsamples; j ++) {
      for (k = 0; k < r->ysamples; k ++) {
	r->conf_min[j][k] = rjmcmc_head_from_histogram(r->hist[j][k],
						       r->zmin,
						       r->zmax,
						       r->zsamples,
						       conf_samples);
	r->conf_max[j][k] = rjmcmc_tail_from_histogram(r->hist[j][k],
						       r->zmin,
						       r->zmax,
						       r->zsamples,
						       conf_samples);
      }
    }
  }
}

/*
 * Getting results
 */

const int *
resultset2d_get_propose(resultset2d_t *r,
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
resultset2d_get_accept(resultset2d_t *r,
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
resultset2d_get_misfit(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->misfit;
}

const double *
resultset2d_get_lambda(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->lambda;
}

const int *
resultset2d_get_partitions(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->partitions;
}

const double **
resultset2d_get_mean(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return (const double **)r->mean;
}

const int **
resultset2d_get_centres(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return (const int**)r->centres;
}

const double **
resultset2d_get_median(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return (const double **)r->median;
}

const double **
resultset2d_get_mode(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return (const double **)r->mode;
}

const double **
resultset2d_get_credible_min(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return (const double **)r->conf_min;
}
		     
const double **
resultset2d_get_credible_max(resultset2d_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return (const double **)r->conf_max;
}

void
resultset2d_fill_xcoord_vector(resultset2d_t *r,
			       double *x,
			       int *l)
{
  int i;
  int e;

  if (r == NULL) {
    return;
  }

  e = *l;
  if (e > r->xsamples) {
    e = r->xsamples;
    *l = e;
  }

  for (i = 0; i < e; i ++) {
    x[i] = r->xmin + (r->xmax - r->xmin) * ((double)i + 0.5)/(double)(r->xsamples - 1);
  }
}

void
resultset2d_fill_ycoord_vector(resultset2d_t *r,
			       double *y,
			       int *l)
{
  int i;
  int e;

  if (r == NULL) {
    return;
  }

  e = *l;
  if (e > r->ysamples) {
    e = r->ysamples;
    *l = e;
  }

  for (i = 0; i < e; i ++) {
    y[i] = r->ymin + (r->ymax - r->ymin) * ((double)i + 0.5)/(double)(r->ysamples - 1);
  }
}

void
resultset2d_fill_zcoord_vector(resultset2d_t *r,
			       double *z,
			       int *l)
{
  int i;
  int e;

  e = *l;
  if (e > r->zsamples) {
    e = r->zsamples;
    *l = e;
  }

  for (i = 0; i < e; i ++) {
    z[i] = r->zmin + (r->zmax - r->zmin) * ((double)i + 0.5)/(double)(r->zsamples - 1);
  }

}

