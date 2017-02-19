
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "rjmcmc/part1d_regression_rj.h"

#include "rjmcmc/position_map1d.h"
#include "rjmcmc/curvefit.h"

#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_defines.h"

static const double DEFAULT_AUTO_Z = 3.0;

struct _model {
  
  double **a;
  int *order;

  double lambda;

  double **pk;
  double **kcdf;

  double **prior_product;
  double *ppratio;
};

typedef struct _model model_t;
  
struct _part1d_regression_rj {

  /*
   * Constant parameters
   */
  int min_partitions;
  int max_partitions;
  int max_order;
  int ndatasets;
  double xmin;
  double xmax;
  double pd;
  double auto_z;

  /*
   * Varying parameter
   */
  int npartitions;

  /*
   * Coordinates of the partition boundaries
   */
  position_map1d_t *p;

  /*
   * Models for each dataset
   */
  model_t *models;

  /*
   * Curve fitting data
   */
  curvefit_result_t *cf;

  /*
   * detCm, mean misfit etc scratch arrays for recalculating distribution
   * of k.
   */
  double **mean;
  double **sigma;
  double *autoprior;
  double *detCm;
  double *mean_misfit;
  double **S;

};

static model_t *models_create(int max_part,
			      int max_order,
			      int ndatasets);
static void models_destroy(int max_part,
			   int max_order,
			   int ndatasets,
			   model_t *m);
static void models_clone(int max_partitions,
			 int max_order,
			 int ndatasets,
			 const model_t *src,
			 model_t *dst);

static void models_delete(int max_partitions,
			  int max_order,
			  int ndatasets,
			  int del_iy,
			  int npart,
			  model_t *m);

static int update_partition(part1d_regression_rj_t *proposed,
			    const dataset1d_t *d,
			    int di,
			    int pi,
			    int xi,
			    int xj,
			    double xl,
			    double xr,
			    rjmcmc_uniform_rand_t random,
			    rjmcmc_normal_rand_t normal,
			    double *prob);

static double value_at(const part1d_regression_rj_t *current,
		       int di,
		       double x);

static int partitions_valid(part1d_regression_rj_t *p,
			    const dataset1d_t **datasets,
			    int ndatasets);

part1d_regression_rj_t *
part1d_regression_rj_create(int min_partitions,
			    int max_partitions,
			    int max_order,
			    int ndatasets,
			    double xmin, 
			    double xmax,
			    double pd)
{
  part1d_regression_rj_t *r;
  
  r = (part1d_regression_rj_t*)malloc(sizeof(part1d_regression_rj_t));
  if (r == NULL) {
    return NULL;
  }

  r->min_partitions = min_partitions;
  r->max_partitions = max_partitions;
  r->max_order = max_order;
  r->ndatasets = ndatasets;

  r->xmin = xmin;
  r->xmax = xmax;
  r->pd = pd;
  r->auto_z = DEFAULT_AUTO_Z; 

  r->npartitions = 0;

  r->p = position_map1d_create(max_partitions, xmin, xmax);
  if (r->p == NULL) {
    return NULL;
  }

  r->models = models_create(max_partitions,
			    max_order,
			    ndatasets);
  if (r->models == NULL) {
    return NULL;
  }

  r->cf = curvefit_create(max_order);
  if (r->cf == NULL) {
    return NULL;
  }

  r->mean = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (r->mean == NULL) {
    return NULL;
  }

  r->sigma = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (r->sigma == NULL) {
    return NULL;
  }

  r->autoprior = rjmcmc_create_array_1d(max_order + 1);
  if (r->autoprior == NULL) {
    return NULL;
  }

  r->detCm = rjmcmc_create_array_1d(max_order + 1);
  if (r->detCm == NULL) {
    return NULL;
  }

  r->mean_misfit = rjmcmc_create_array_1d(max_order + 1);
  if (r->mean_misfit == NULL) {
    return NULL;
  }

  r->S = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (r->S == NULL) {
    return NULL;
  }

  return r;
}

void
part1d_regression_rj_destroy(part1d_regression_rj_t *p)
{
  if (p != NULL) {

    position_map1d_destroy(p->p);

    models_destroy(p->max_partitions,
		   p->max_order,
		   p->ndatasets,
		   p->models);

    curvefit_destroy(p->cf);

    rjmcmc_destroy_array_2d(p->max_order + 1, p->mean);
    rjmcmc_destroy_array_2d(p->max_order + 1, p->sigma);

    rjmcmc_destroy_array_1d(p->autoprior);
    rjmcmc_destroy_array_1d(p->detCm);
    rjmcmc_destroy_array_1d(p->mean_misfit);
    rjmcmc_destroy_array_2d(p->max_order + 1, p->S);

    free(p);

  }
}

void
part1d_regression_rj_clone(const part1d_regression_rj_t *src,
			   part1d_regression_rj_t *dst)
{
  RJMCMC_NULLCHECKVOID(src, "part1d_regression_rj_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "part1d_regression_rj_clone: null dst\n");

  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions, 
			    "part1d_regression_rj_clone: size mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_order != dst->max_order, 
			    "part1d_regression_rj_clone: order mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->ndatasets != dst->ndatasets, 
			    "part1d_regression_rj_clone: count mismatch\n");

  position_map1d_clone(src->p, dst->p);
  models_clone(src->max_partitions,
	       src->max_order,
	       src->ndatasets,
	       src->models, dst->models);
  dst->npartitions = src->npartitions;
}

int 
part1d_regression_rj_initialize(part1d_regression_rj_t *p,
				const dataset1d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal)
{
  int npart;
  int pi;
  int di;

  double x;
  double xl;
  double xr;

  int xi;
  int xj;

  int order;

  const dataset1d_t *data;

  double curve_prob;

  int n;
  int i;

  npart = 2;

  for (pi = 2; pi < npart; pi ++) {
    x = p->xmin + (p->xmax - p->xmin)*random();
    position_map1d_insert(p->p, 
			x,
			pi);

    /*
     * Ensure that we create a valid set of initial partitions 
     * by deleting any that are too close to another partition
     * boundary.
     */

    if (!partitions_valid(p, 
			  datasets,
			  ndatasets)) {
      position_map1d_delete(p->p,
			  x,
			  pi);
      npart --;
      pi --;
    }
  }

  p->npartitions = npart;
  npart --;
  xl = position_map1d_position_of_index(p->p, 0);


  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    if (data->lambdastd > 0.0) {
      p->models[di].lambda = 
	rjmcmc_random_choose_double(data->lambdamin,
				    data->lambdamax,
				    random);
    }
  }

  for (i = 0; i < npart; i ++) {

    pi = position_map1d_predecessor_of_point(p->p, xl);

    /* Create the initial models */
    xr = position_map1d_next_position(p->p, xl);
    if (xr < xl) {
      rjmcmc_error(
	      "part1d_regression_rj_initialize: "
	      "failed to get next point\n");
      return -1;
    }

    for (di = 0; di < ndatasets; di ++) {
      data = datasets[di];

      order = rjmcmc_random_choose_int(0, p->max_order, random);
      p->models[di].order[pi] = order;
      n = dataset1d_range(data, xl, xr, &xi, &xj);

      if (update_partition(p,
			   data,
			   di,
			   pi,
			   xi,
			   xj,
			   xl,
			   xr,
			   random,
			   normal,
			   &curve_prob) < 0) {
	rjmcmc_error("part1d_regression_rj_initialize: "
		     "failed to update partition\n");
	return -1;
      }

    }
  }

  return 0;
}

static double
value_at(const part1d_regression_rj_t *current,
	 int di,
	 double x)
{
  int iy;
  int order;
  const double *a;

  double y;

  iy = position_map1d_predecessor_of_point(current->p, x);
  if (iy < 0) {
    return -DBL_MAX;
  }

  if (iy == 1) {
    iy = position_map1d_predecessor_of_index(current->p, iy);
    if (iy < 0 || iy == 1) {
      return -DBL_MAX;
    }
  }


  a = current->models[di].a[iy];
  order = current->models[di].order[iy];

  y = rjmcmc_polynomial_value(a, order + 1, x);

  return y;
}

double
part1d_regression_rj_value_at(const part1d_regression_rj_t *current,
			      double x)
{
  return value_at(current, 0, x);
}

int
part1d_regression_rj_evaluate(const part1d_regression_rj_t *current,
			      int di,
			      double xmin,
			      double xmax,
			      int nsamples,
			      double *samples)
{
  int i;

  double xi;

  for (i = 0; i < nsamples; i ++) {
    xi = xmin + (double)i * (xmax - xmin)/(double)(nsamples - 1);
    samples[i] = value_at(current, di, xi);
  }
   
  return 0;  
}

double
part1d_regression_rj_misfit(part1d_regression_rj_t *p,
			    const dataset1d_t **datasets,
			    int ndatasets)
{
  int i;
  int di;

  double y;
  double n;
  double dy;
  double sum;
  double dsum;
  double l2;

  const dataset1d_t *dataset;

  sum = 0.0;
  
  for (di = 0; di < ndatasets; di ++) {
    
    dsum = 0.0;
    dataset = datasets[di];

    l2 = 1.0;
    if (dataset->lambdastd > 0.0) {
      l2 = p->models[di].lambda * p->models[di].lambda;
    }
    
    for (i = 0; i < dataset->npoints; i ++) {
      
      y = value_at(p, di, dataset->points[i].x);
      if (y > dataset->ymax || y < dataset->ymin) {
	return DBL_MAX;
      }
      
      dy = y - dataset->points[i].y;
      n = dataset->points[i].n;
      dsum += (dy*dy)/(2.0*n*n*l2);
    }

    sum += dsum;
  }

  return sum;

}

int
part1d_regression_rj_propose_birth(const part1d_regression_rj_t *current,
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *birth_prob)
{
  double new_x;
  int new_iy;

  int prev_iy;
  double prev_x;

  double prob;
  int di;

  int xi;
  int xj;

  double new_x_right;

  const dataset1d_t *data;

  int n;

  double curve_prob;

  if (current->npartitions == current->max_partitions) {
    /* rjmcmc_error( */
    /* 	    "part1d_regression_rj_propose_birth: " */
    /* 	    "%d %d\n",  */
    /* 	    current->npartitions, */
    /* 	    current->max_partitions); */
    return 0;
  }

  part1d_regression_rj_clone(current, proposed);

  new_x = random() * (proposed->xmax - proposed->xmin) + proposed->xmin;
  new_iy = proposed->npartitions;

  if (position_map1d_insert(proposed->p, new_x, new_iy) < 0) {
    rjmcmc_error(
	    "part1d_regression_rj_propose_birth: "
	    "failed to add new point\n");
    return 0;
  }

  if (!partitions_valid(proposed,
			datasets,
			ndatasets)) {
    return 0;
  }
  
  proposed->npartitions ++;

  new_x_right = position_map1d_next_position(proposed->p, new_x);
  if (new_x_right < new_x) {
    rjmcmc_error("part1d_regression_rj_propose_birth: "
		 "failed to find right extent of new point\n");
    return 0;
  }
    
  prev_iy = position_map1d_predecessor_of_index(proposed->p, new_iy);
  if (prev_iy < 0) {
    rjmcmc_error("part1d_regression_rj_propose_birth: "
		 "failed to find predecessor\n");
    return 0;
  }

  /*
   * Store the prior/proposal ratio as it will be overwritten when
   * the curves are resampled.
   */
  prob = 1.0;
  for (di = 0; di < ndatasets; di ++) {
    prob /= proposed->models[di].ppratio[prev_iy];
  }

  /*
   * Update the partitions for the new partition (b)
   */

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    n = dataset1d_range(data, new_x, new_x_right, &xi, &xj);
    if (n <= 1) {
      return 0;
    }

    if (update_partition(proposed,
			 data,
			 di,
			 new_iy,
			 xi,
			 xj,
			 new_x,
			 new_x_right,
			 random,
			 normal,
			 &curve_prob) < 0) {
      rjmcmc_error(
	      "part1d_regression_rj_propose_birth: "
	      "failed to update new partition\n");
      return 0;
    }

    prob *= curve_prob;
  }
  
  /*
   * Update the partition that the new partition affected (c)
   */

  prev_x = position_map1d_position_of_index(proposed->p, prev_iy);

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    n = dataset1d_range(data, prev_x, new_x, &xi, &xj);
    if (n <= 1) {
      return 0;
    }

    if (update_partition(proposed,
			 data,
			 di,
			 prev_iy,
			 xi,
			 xj,
			 prev_x,
			 new_x,
			 random,
			 normal,
			 &curve_prob) < 0) {
      rjmcmc_error(
	      "part1d_regression_rj_propose_birth: "
	      "failed to update new partition\n");
      return 0;
    }

    prob *= curve_prob;
  }

  *birth_prob = prob;

  return -1;
}

int 
part1d_regression_rj_propose_death(const part1d_regression_rj_t *current,
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *death_prob)
{
  int del_iy;
  double deleted_pos;
  int new_iy;
  int di;
  
  double prob;
  double curve_prob;
  
  const dataset1d_t *data;
  
  int order;
  int n;
  
  double xl;
  double xr;

  int xi;
  int xj;

  part1d_regression_rj_clone(current, proposed);

  if (proposed->npartitions <= 2 || proposed->npartitions <= proposed->min_partitions) {
    /* Can't remove any more points */
    /* rjmcmc_error( */
    /* 	    "part1d_regression_rj_propose_death:" */
    /* 	    "too few partitions %d %d\n", */
    /* 	    proposed->npartitions, */
    /* 	    proposed->min_partitions); */

    return 0;
  }
  
  /* Remove one value at random, note that we can't remove the two endpoints
   * which are always located and iy indices of 0 and 1, hence the 2's in 
   * the random choice of index to remove
   */
  del_iy = rjmcmc_random_choose_int(2, proposed->npartitions - 1,
				    random);

  /*
   * Record the prior/proposal ratio for the deleted partition
   */
  prob = 1.0;
  for (di = 0; di < proposed->ndatasets; di ++) {
    prob /= proposed->models[di].ppratio[del_iy];
  }
    
  deleted_pos = position_map1d_position_of_index(proposed->p, del_iy);
  
  if (position_map1d_delete(proposed->p, deleted_pos, del_iy) < 0) {
    rjmcmc_error("part1d_regression_rj_propose_death: failed to delete position\n");
    return 0;
  }

  models_delete(proposed->max_partitions,
		proposed->max_order,
		proposed->ndatasets,
		del_iy,
		proposed->npartitions,
		proposed->models);
  proposed->npartitions --;

  new_iy = position_map1d_predecessor_of_point(proposed->p, deleted_pos);
  if (new_iy < 0) {
    rjmcmc_error("part1d_regression_rj_propose_death: failed to find predecessor\n");
    return 0;
  }

  /*
   * Record the prior/proposal ratio for the left neighbour partition
   */
  prob = 1.0;
  for (di = 0; di < proposed->ndatasets; di ++) {
    prob /= proposed->models[di].ppratio[new_iy];
  }


  prob = 1.0;

  xl = position_map1d_position_of_index(proposed->p, new_iy);
  xr = position_map1d_next_position(proposed->p, xl);


  /*
   * Update the partition that contained the deleted partition
   */

  for (di = 0; di < ndatasets; di ++) {

    data = datasets[di];
    order = proposed->models[di].order[new_iy];

    n = dataset1d_range(data, xl, xr, &xi, &xj);

    if (update_partition(proposed,
			 data,
			 di,
			 new_iy,
			 xi,
			 xj,
			 xl,
			 xr,
			 random,
			 normal,
			 &curve_prob) < 0) {
      rjmcmc_error(
	      "part1d_regression_rj_propose_death: "
	      "failed to update partition\n");
      return 0;
    }

    prob *= curve_prob;
  }

  *death_prob = prob;
  /* printf("death: %f\n", *death_prob); */
  /* *death_prob = 1.0; */
  return -1;
}

int 
part1d_regression_rj_propose_value(const part1d_regression_rj_t *current, 
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *value_prob)
{
  int di;
  int iy;

  double xl;
  double xr;

  int n;
  int xi;
  int xj;

  const dataset1d_t *data;
  int order;

  double prob;
  double curve_prob;

  part1d_regression_rj_clone(current, proposed);

  if (proposed->ndatasets == 1) {
    di = 0;
  } else {
    di = rjmcmc_random_choose_int(0, proposed->ndatasets - 1, random);
  }

  /* We want to select either partition 0, or 2 .. npartitions as partition
     1 has no bearing on the result */
  iy = rjmcmc_random_choose_int(0, proposed->npartitions - 2, random);
  if (iy > 0) {
    iy ++;
  }

  xl = position_map1d_position_of_index(proposed->p, iy);
  xr = position_map1d_next_position(proposed->p, xl);

  order = proposed->models[di].order[iy];
  data = datasets[di];

  n = dataset1d_range(datasets[di], xl, xr, &xi, &xj);
  if (n <= 0) {
    return 0;
  }

  prob = 1.0/proposed->models[di].ppratio[iy];

  if (update_partition(proposed,
		       data,
		       di,
		       iy,
		       xi,
		       xj,
		       xl,
		       xr,
		       random,
		       normal,
		       &curve_prob) < 0) {
    rjmcmc_error("part1d_regression_rj_propose_value: failed to update part\n");
    return 0;
  }

  *value_prob = prob * curve_prob;

  return 1;
}

int 
part1d_regression_rj_propose_lambda(const part1d_regression_rj_t *current,
				    part1d_regression_rj_t *proposed,
				    const dataset1d_t **datasets,
				    int ndatasets,
				    rjmcmc_uniform_rand_t random,
				    rjmcmc_normal_rand_t normal,
				    double *lambda_prob)
{
  int di;
  double new_s;
  const dataset1d_t *data;

  di = (int)((double)ndatasets * random());
  data = datasets[di];

  part1d_regression_rj_clone(current, proposed);
  new_s = proposed->models[di].lambda + normal() * data->lambdastd;

  if (new_s < data->lambdamin ||
      new_s > data->lambdamax) {
    return 0;
  }

  *lambda_prob = pow(proposed->models[di].lambda/new_s,
		     data->npoints);
  proposed->models[di].lambda = new_s;
  return -1;
}

int 
part1d_regression_rj_propose_move(const part1d_regression_rj_t *current,
				  part1d_regression_rj_t *proposed,
				  const dataset1d_t **datasets,
				  int ndatasets,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *move_prob)
{
  int iy;

  int iyop;
  int iynp;

  double old_x;
  double new_x;

  double xl, xr;

  int di;
  int xi;
  int xj;
  int n;
  
  double prob;
  double curve_prob;

  const dataset1d_t *data;

  /* The 2 end points can't be moved, if this is all there is then give up */
  if (current->npartitions <= 2) {
    return 0;
  }

  part1d_regression_rj_clone(current, proposed);

  iy = rjmcmc_random_choose_int(2, proposed->npartitions - 1, random);

  old_x = position_map1d_position_of_index(proposed->p, iy);
  new_x = old_x + normal() * proposed->pd;

  if (new_x <= proposed->xmin ||
      new_x >= proposed->xmax) {
    return 0;
  }

  iyop = position_map1d_predecessor_of_index(proposed->p, iy);
  if (iyop < 0) {
    rjmcmc_error(
	    "part1d_regression_rj_propose_move: "
	    "failed to find old precedessor of point\n");
    return 0;
  }

  /*
   * Store the prior/proposal ratio for the 2 regions that will definitely
   * be affected by the move (a third will be affected if the move is large).
   */
  prob = 1.0;
  for (di = 0; di < proposed->ndatasets; di ++) {
    prob /= proposed->models[di].ppratio[iy];
    prob /= proposed->models[di].ppratio[iyop];
  }

  if (position_map1d_move(proposed->p,
			old_x,
			new_x) < 0) {
    rjmcmc_error("part1d_regression_rj_propose_move: "
		 "failed to move point\n");
    return 0;
  }

  if (!partitions_valid(proposed,
			datasets,
			ndatasets)) {
    return 0;
  }

  iynp = position_map1d_predecessor_of_index(proposed->p, iy);
  if (iynp < 0) {
    rjmcmc_error("part1d_regression_rj_propose_move: "
		 "failed to find new predecessor predecessor\n");
    return 0;
  }

  /*
   * We need to update the partitions which have been affected by the
   * move, this will be 2 partitions if it is a small move (the boundary
   * moved hasn't crossed another boundary) or 3 partitions for a large
   * move.
   */

  /* First the moved partition */
  xl = position_map1d_position_of_index(proposed->p, iy);
  xr = position_map1d_next_position(proposed->p, xl);

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    n = dataset1d_range(datasets[di], xl, xr, &xi, &xj);
    if (n <= 1) {
      return 0;
    }
    
    if (update_partition(proposed,
  			 data,
  			 di,
  			 iy,
  			 xi,
  			 xj,
			 xl,
			 xr,
  			 random,
  			 normal,
  			 &curve_prob) < 0) {
      rjmcmc_error("part1d_regression_rj_propose_move: failed to update part\n");
      return 0;
    }

    prob *= curve_prob;
  }

  /* Second, the partition that originally contained the moved partition */
  xl = position_map1d_position_of_index(proposed->p, iyop);
  xr = position_map1d_next_position(proposed->p, xl);

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    n = dataset1d_range(datasets[di], xl, xr, &xi, &xj);
    if (n <= 1) {
      return 0;
    }
    
    if (update_partition(proposed,
  			 data,
  			 di,
  			 iyop,
  			 xi,
  			 xj,
			 xl,
			 xr,
  			 random,
  			 normal,
  			 &curve_prob) < 0) {
      rjmcmc_error("part1d_regression_rj_propose_move: failed to update part\n");
      return 0;
    }

    prob *= curve_prob;
  }

  /*
   * Finally the partition the moved partition moved into if this was a large 
   * move
   */
  if (iyop != iynp) {

    /*
     * Add in this prior/proposed ratio
     */
    for (di = 0; di < proposed->ndatasets; di ++) {
      prob /= proposed->models[di].ppratio[iynp];
    }

    xl = position_map1d_position_of_index(proposed->p, iynp);
    xr = position_map1d_next_position(proposed->p, xl);
    
    for (di = 0; di < ndatasets; di ++) {
      data = datasets[di];
      
      n = dataset1d_range(datasets[di], xl, xr, &xi, &xj);
      if (n <= 1) {
	return 0;
      }
      
      if (update_partition(proposed,
			   data,
			   di,
			   iynp,
			   xi,
			   xj,
			   xl,
			   xr,
			   random,
			   normal,
			   &curve_prob) < 0) {
	rjmcmc_error("part1d_regression_rj_propose_move: failed to update part\n");
	return 0;
      }

      prob *= curve_prob;
    }
  }

  *move_prob = prob;

  return 1;
}

/*
 * Internal methods
 */

static model_t *models_create(int max_part,
			      int max_order,
			      int ndatasets)
{
  model_t *m;
  int di;

  m = (model_t*)malloc(sizeof(model_t) * ndatasets);
  if (m == NULL) {
    return NULL;
  }
  
  for (di = 0; di < ndatasets; di ++) {
    m[di].a = rjmcmc_create_array_2d(max_part, max_order + 1);
    if (m[di].a == NULL) {
      return NULL;
    }

    m[di].order = rjmcmc_create_int_array_1d(max_part);
    if (m[di].order == NULL) {
      return NULL;
    }

    m[di].pk = rjmcmc_create_array_2d(max_part, max_order + 1);
    if (m[di].pk == NULL) {
      return NULL;
    }

    m[di].kcdf = rjmcmc_create_array_2d(max_part, max_order + 1);
    if (m[di].kcdf == NULL) {
      return NULL;
    }
    
    m[di].prior_product = rjmcmc_create_array_2d(max_part, max_order + 1);
    if (m[di].prior_product == NULL) {
      return NULL;
    }

    m[di].ppratio = rjmcmc_create_array_1d(max_part);
    if (m[di].ppratio == NULL) {
      return NULL;
    }
  }


  return m;
}

static void models_destroy(int max_part,
			   int max_order,
			   int ndatasets,
			   model_t *m)
{
  int di;

  if (m != NULL) {
    
    for (di = 0; di < ndatasets; di ++) {
      rjmcmc_destroy_array_2d(max_part, m[di].a);
      rjmcmc_destroy_int_array_1d(m[di].order);

      rjmcmc_destroy_array_2d(max_part, m[di].pk);
      rjmcmc_destroy_array_2d(max_part, m[di].kcdf);
      rjmcmc_destroy_array_2d(max_part, m[di].prior_product);
      rjmcmc_destroy_array_1d(m[di].ppratio);
    }

    free(m);
  }
}

static void models_clone(int max_partitions,
			 int max_order,
			 int ndatasets,
			 const model_t *src,
			 model_t *dst)
{
  int di;
  int pi;
  int ci;
  int oi;

  RJMCMC_NULLCHECKVOID(src, "models_clone: src null\n");
  RJMCMC_NULLCHECKVOID(dst, "models_clone: dst null\n");

  for (di = 0; di < ndatasets; di ++) {
    for (pi = 0; pi < max_partitions; pi ++) {
      for (ci = 0; ci <= max_order; ci ++) {
	dst[di].a[pi][ci] = src[di].a[pi][ci];
      }
      dst[di].order[pi] = src[di].order[pi];
      
      for (oi = 0; oi <= max_order; oi ++) {
	dst[di].pk[pi][oi] = src[di].pk[pi][oi];
	dst[di].kcdf[pi][oi] = src[di].kcdf[pi][oi];
	dst[di].prior_product[pi][oi] = src[di].prior_product[pi][oi];
      }

      dst[di].ppratio[pi] = src[di].ppratio[pi];
    }

    dst[di].lambda = src[di].lambda;
  }
       

}

static void models_delete(int max_partitions,
			  int max_order,
			  int ndatasets,
			  int del_iy,
			  int npart,
			  model_t *m)
{
  int di;
  int pi;
  int ci;
  int oi;

  for (di = 0; di < ndatasets; di ++) {
    for (pi = del_iy + 1; pi < npart; pi ++) {
      for (ci = 0; ci <= max_order; ci ++) {
	m[di].a[pi - 1][ci] = m[di].a[pi][ci];
      }
      m[di].order[pi - 1] = m[di].order[pi];

      for (oi = 0; oi <= max_order; oi ++) {
	m[di].pk[pi - 1][oi] = m[di].pk[pi][oi];
	m[di].kcdf[pi - 1][oi] = m[di].kcdf[pi][oi];
	m[di].prior_product[pi - 1][oi] = m[di].prior_product[pi][oi];
      }

      m[di].ppratio[pi - 1] = m[di].ppratio[pi];
    }
  }
}

static int resample_partition(part1d_regression_rj_t *proposed,
			      const dataset1d_t *data,
			      int di,
			      int pi,
			      int xi,
			      int xj,
			      double xl,
			      double xr,
			      rjmcmc_uniform_rand_t random,
			      rjmcmc_uniform_rand_t normal,
			      double *prob)
{
  int order;
  double curve_prob;
  int i;

  /*
   * Choose the new order
   */
  order = rjmcmc_random_choose_interval(proposed->models[di].kcdf[pi],
					proposed->max_order + 1,
					random);
  /*
   * Compute the best fit given the order
   */
  if (curvefit_compute(data,
		       xi,
		       xj,
		       order,
		       proposed->cf) < 0) {
    rjmcmc_error(
	    "update_partition: failed to compute curvefit (%d %d %d)\n",
	    xi,
	    xj,
	    order);
    return -1;
  }


  /*
   * Sample a random curve using a multivariate normal distribution 
   * from the best fit curve.
   */
  if (curvefit_sample(proposed->cf,
		      normal,
		      proposed->models[di].a[pi],
		      order + 1,
		      &curve_prob) < 0) {
    rjmcmc_error("update_partition: failed to sample curve\n");
    return -1;
  }

  proposed->models[di].order[pi] = order;

  /*
   * The value here is the ratio of the prior over the proposal.
   */
  proposed->models[di].ppratio[pi] = 
    proposed->models[di].prior_product[pi][order]/
    (proposed->models[di].pk[pi][order] * curve_prob);
  if (proposed->models[di].ppratio[pi] == 0.0) {
    rjmcmc_error("ppratio underflow: %g %g %g\n",
	    proposed->models[di].prior_product[pi][order],
	    proposed->models[di].pk[pi][order],
	    curve_prob);
  }

  *prob = proposed->models[di].ppratio[pi];

  return 0;
}
			     
static int update_partition(part1d_regression_rj_t *proposed,
			    const dataset1d_t *data,
			    int di,
			    int pi,
			    int xi,
			    int xj,
			    double xl,
			    double xr,
			    rjmcmc_uniform_rand_t random,
			    rjmcmc_normal_rand_t normal,
			    double *prob)
{
  int oi;
  int i;

  if (curvefit_evaluate_pk(proposed->cf,
			   data,
			   xi,
			   xj,
			   proposed->max_order,
			   NULL,
			   proposed->auto_z,
			   proposed->mean_misfit,
			   proposed->detCm,
			   proposed->autoprior,
			   proposed->S,
			   proposed->models[di].pk[pi],
			   proposed->models[di].kcdf[pi],
			   proposed->mean,
			   proposed->sigma) < 0) {
    rjmcmc_error("update_partition: failed to determine pk\n");
    return -1;
  }

  /*
   * Now update the prior product for each order
   */
  for (oi = 0; oi <= proposed->max_order; oi ++) {
    proposed->models[di].prior_product[pi][oi] = 1.0;

    for (i = 0; i <= oi; i ++) {
      proposed->models[di].prior_product[pi][oi] *= 
	2.0 * DEFAULT_AUTO_Z * proposed->sigma[oi][i];
    }
  }
  
  if (resample_partition(proposed,
			 data,
			 di,
			 pi,
			 xi,
			 xj,
			 xl,
			 xr,
			 random,
			 normal,
			 prob) < 0) {
    rjmcmc_error("update_partition: failed to resample curve\n");
    return -1;
  }
  
  return 0;
}

int
part1d_regression_rj_partitions(const part1d_regression_rj_t *current)
{
  return current->npartitions - 1;
}

double
part1d_regression_rj_lambda(const part1d_regression_rj_t *current,
			    int di)
{
  return current->models[di].lambda;
}

int
part1d_regression_rj_order(const part1d_regression_rj_t *current,
			   int di)
{
  return current->models[di].order[0];
}

int 
part1d_regression_rj_order_sum(const part1d_regression_rj_t *current,
			       int di)
{
  int i;
  int s;
  
  s = 0;
  for (i = 0; i < current->npartitions; i ++) {
    s += current->models[di].order[i] + 1;
  }

  return s;
}
			       
double
part1d_regression_rj_partition_position(const part1d_regression_rj_t *current,
					int pi)
{
  return position_map1d_position_of_index(current->p, pi);
}

struct _partitions_valid_data {
  const part1d_regression_rj_t *current;
  const dataset1d_t **datasets;
  int ndatasets;
  int invalid_count;
  int max_order;
};

static int partitions_valid_cb(void *user_arg,
			       double xmin,
			       double xmax,
			       int iy,
			       int riy)
{
  struct _partitions_valid_data *d = (struct _partitions_valid_data*)user_arg;
  int i;
  const dataset1d_t *data = d->datasets[0];
  int c;

  c = 0;
  for (i = 0; i < data->npoints; i ++) {
    if (data->points[i].x >= xmin && 
	data->points[i].x <= xmax) {
      c ++;
    }
  }
   
  if (c <= d->max_order) {
    d->invalid_count ++;
  }

  /* if (d->current->models[0].ppratio[iy] <= 0.0) { */
  /*   printf("ppratio invalid: %g\n", d->current->models[0].ppratio[iy]); */
  /* } */
  
  return 0;
}

static int partitions_valid(part1d_regression_rj_t *p,
			    const dataset1d_t **datasets,
			    int ndatasets)
{
  struct _partitions_valid_data d;

  d.current = p;
  d.datasets = datasets;
  d.ndatasets = ndatasets;
  d.invalid_count = 0;
  d.max_order = p->max_order;

  if (position_map1d_traverse_intervals(p->p,
					partitions_valid_cb,
					&d) < 0) {
    rjmcmc_error("partitions_valid: failed to traverse intervals\n");
    return 0;
  }

  return d.invalid_count == 0;
}
