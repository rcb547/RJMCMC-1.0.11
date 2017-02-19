
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "rjmcmc/part1d_zero.h"

#include "rjmcmc/position_map1d.h"
#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_defines.h"

struct _model {
  
  double *a;
  double lambda;

  double *mean;
  double *var;

};

typedef struct _model model_t;
  
struct _part1d_zero {

  /*
   * Constant parameters
   */
  int min_partitions;
  int max_partitions;
  int ndatasets;

  double xmin;
  double xmax;
  double pd;

  /*
   * Varying parameter
   */
  int npartitions;

  /*
   * Coordinates of the partition boundaries
   */
  position_map1d_t *p;

  /*
   * Models for each partition
   */
  model_t *models;

};

static model_t *models_create(int max_part,
			      int ndatasets);
static void models_destroy(int max_part,
			   int ndatasets,
			   model_t *m);
static void models_clone(int max_partitions,
			 int ndatasets,
			 const model_t *src,
			 model_t *dst);

static void models_delete(int max_partitions,
			  int ndatasets,
			  int del_iy,
			  int npart,
			  model_t *m);

static double value_at(const part1d_zero_t *current,
		       int di,
		       double x);

part1d_zero_t *
part1d_zero_create(int min_partitions,
		   int max_partitions,
		   int ndatasets,
		   double xmin, 
		   double xmax,
		   double pd)
{
  part1d_zero_t *r;
  
  r = (part1d_zero_t*)malloc(sizeof(part1d_zero_t));
  if (r == NULL) {
    return NULL;
  }

  r->min_partitions = min_partitions;
  if (r->min_partitions < 2) {
    r->min_partitions = 2;
  }

  r->max_partitions = max_partitions;
  r->ndatasets = ndatasets;

  r->xmin = xmin;
  r->xmax = xmax;
  r->pd = pd;

  r->npartitions = 0;

  r->p = position_map1d_create(max_partitions, xmin, xmax);
  if (r->p == NULL) {
    return NULL;
  }

  r->models = models_create(max_partitions,
			    ndatasets);
  if (r->models == NULL) {
    return NULL;
  }

  return r;
}

void
part1d_zero_destroy(part1d_zero_t *p)
{
  if (p != NULL) {

    position_map1d_destroy(p->p);

    models_destroy(p->max_partitions,
		   p->ndatasets,
		   p->models);

    free(p);

  }
}

void
part1d_zero_clone(const part1d_zero_t *src,
			part1d_zero_t *dst)
{
  RJMCMC_NULLCHECKVOID(src, "part1d_zero_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "part1d_zero_clone: null dst\n");

  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "part1d_zero_clone: size mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->ndatasets != dst->ndatasets,
			    "part1d_zero_clone: count mismatch\n");

  position_map1d_clone(src->p, dst->p);
  models_clone(src->max_partitions,
	       src->ndatasets,
	       src->models, dst->models);
  dst->npartitions = src->npartitions;
}

struct initialize_cb_data {
  part1d_zero_t *p;
  const dataset1d_t **datasets;
  int ndatasets;
  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;
};

static int
initialize_cb(void *user_arg,
	      double xmin,
	      double xmax,
	      int iy,
	      int riy)
{
  struct initialize_cb_data *d = (struct initialize_cb_data *)user_arg;
  int di;

  const dataset1d_t *data;

  int n;
  int xi;
  int xj;


  for (di = 0; di < d->ndatasets; di ++) {
    data = d->datasets[di];

    n = dataset1d_range(data, xmin, xmax, &xi, &xj);
    if (n < 2) {
      d->p->models[di].mean[iy] = 0.0;
      d->p->models[di].var[iy] = 0.0;
    } else {
      dataset1d_mean_variance(data, xi, xj, 
			      &(d->p->models[di].mean[iy]),
			      &(d->p->models[di].var[iy]));

     
    }
    
    d->p->models[di].a[iy] = 
      rjmcmc_random_choose_double(data->ymin, 
				  data->ymax,
				  d->random);

  }
  
  return 0;
}

int 
part1d_zero_initialize(part1d_zero_t *p,
		       const dataset1d_t **datasets,
		       int ndatasets,
		       rjmcmc_uniform_rand_t random,
		       rjmcmc_normal_rand_t normal)
{
  int npart;
  int pi;
  int di;

  double x;
  int i;

  struct initialize_cb_data cb_data;
  const dataset1d_t *data;

  npart = 2;
  
  for (pi = 2; pi < npart; pi ++) {

    x = rjmcmc_random_choose_double(p->xmin,
				    p->xmax,
				    random);
    position_map1d_insert(p->p, 
			x,
			pi);
  }

  p->npartitions = npart;

  cb_data.p = p;
  cb_data.datasets = datasets;
  cb_data.ndatasets = ndatasets;
  cb_data.random = random;
  cb_data.normal = normal;

  position_map1d_traverse_intervals(p->p,
				  initialize_cb,
				  (void*)(&cb_data));

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    if (data->lambdastd > 0.0) {
      p->models[di].lambda = 
	rjmcmc_random_choose_double(data->lambdamin,
				    data->lambdamax,
				    random);
    }
  }
				  
  return 0;
}

double 
part1d_zero_value_at(const part1d_zero_t *current,
		     double x)
{
  return value_at(current, 0, x);
}

static double
value_at(const part1d_zero_t *current,
	 int di,
	 double x)
{
  int pi;

  pi = position_map1d_predecessor_of_point(current->p, x);

  if (pi == 1) {
    /* If the predecessor is index 1 that means we're asking for the
     * value at the right most point and what we actually want is 
     * the predecessors value.
     */
    pi = position_map1d_predecessor_of_index(current->p, pi);
  }

  if (pi < 0 || pi == 1) {
    rjmcmc_error("value_at(part1d_zero): invalid index %d %f\n", pi, x);
    return -DBL_MAX;
  }
  
  return current->models[di].a[pi];
}

int
part1d_zero_evaluate(const part1d_zero_t *current,
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
part1d_zero_misfit(part1d_zero_t *p,
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
      
      dy = y - dataset->points[i].y;
      n = dataset->points[i].n;
      dsum += (dy*dy)/(2.0*n*n*l2);
    }

    sum += dsum;
  }

  return sum;

}

int
part1d_zero_propose_birth(const part1d_zero_t *current,
				part1d_zero_t *proposed,
				const dataset1d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal,
				double *birth_prob)
{
  double new_x;
  int new_iy;

  double prob;
  int di;

  double new_x_right;
  int prev_iy;
  double prev_x;

  int xi;
  int xj;
  int n;

  double phi;

  const dataset1d_t *data;
  
  if (current->npartitions == current->max_partitions) {
    /* rjmcmc_error(*/
    /* 	    "part1d_zero_propose_birth: " */
    /* 	    "%d %d\n", */
    /* 	    current->npartitions, */
    /* 	    current->max_partitions); */
    return 0;
  }
  
  part1d_zero_clone(current, proposed);

  new_x = random() * (proposed->xmax - proposed->xmin) + proposed->xmin;
  new_iy = proposed->npartitions;

  if (position_map1d_insert(proposed->p, new_x, new_iy) < 0) {
    rjmcmc_error(
	    "part1d_zero_propose_birth: "
	    "failed to add new point\n");
    return 0;
  }

  proposed->npartitions ++;

  new_x_right = position_map1d_next_position(proposed->p, new_x);
  if (new_x_right < new_x) {
    rjmcmc_error("part1d_zero_propose_birth: "
		 "failed to find right extent of new point\n");
    return 0;
  }

  prev_iy = position_map1d_predecessor_of_index(proposed->p, new_iy);
  if (prev_iy < 0) {
    rjmcmc_error("part1d_zero_propose_birth: "
		 "failed to find predecessor\n");
    return 0;
  }
    
  prev_x = position_map1d_position_of_index(proposed->p, prev_iy);

  prob = 1.0;

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];
    
    /*
     * First store the probability of the current model based on the
     * store mean/variance (these are about to be overwritten.
     */
    if (proposed->models[di].var[prev_iy] > 0.0) {
      prob *= 
	rjmcmc_gaussian_probability(proposed->models[di].a[prev_iy] - 
				    proposed->models[di].mean[prev_iy],
				    sqrt(proposed->models[di].var[prev_iy]));
    } else {
      prob *= 1.0/(data->ymax - data->ymin);
    }

    /*
     * Determine the mean/variance in the 2 new regions
     */

    n = dataset1d_range(data, prev_x, new_x, &xi, &xj);
    if (n >= 2) {
      dataset1d_mean_variance(data, xi, xj, 
			      &(proposed->models[di].mean[prev_iy]),
			      &(proposed->models[di].var[prev_iy]));
    } else {

      /* If we're birthing an empty partition, abort */
      return 0;
      proposed->models[di].mean[prev_iy] = 0.0;
      proposed->models[di].var[prev_iy] = 0.0;
    }

    n = dataset1d_range(data, new_x, new_x_right, &xi, &xj);
    if (n >= 2) {
      dataset1d_mean_variance(data, xi, xj, 
			      &(proposed->models[di].mean[new_iy]),
			      &(proposed->models[di].var[new_iy]));
    } else {

      /* If we're birthing an empty partition, abort */
      return 0;
      proposed->models[di].mean[new_iy] = 0.0;
      proposed->models[di].var[new_iy] = 0.0;
    }

    /*
     * Sample two new values from the mean/variance in the 2 new 
     * regions
     */

    /* New parition left */
    if (proposed->models[di].var[prev_iy] > 0.0) {
      phi = sqrt(proposed->models[di].var[prev_iy]) *
	normal();
      
      proposed->models[di].a[prev_iy] = 
	proposed->models[di].mean[prev_iy] + phi;
      
      prob /= 
	rjmcmc_gaussian_probability(phi, 
				    sqrt(proposed->models[di].var[prev_iy]));
    } else {
      proposed->models[di].a[prev_iy] = 
	rjmcmc_random_choose_double(data->ymin, data->ymax, random);
      
      prob *= (data->ymax - data->ymin);
    }

    /* New partition right */
    if (proposed->models[di].var[new_iy] > 0.0) {
      phi = sqrt(proposed->models[di].var[new_iy]) *
	normal();
      
      proposed->models[di].a[new_iy] = 
	proposed->models[di].mean[new_iy] + phi;
      
      prob /= 
	rjmcmc_gaussian_probability(phi, 
				    sqrt(proposed->models[di].var[new_iy]));
    } else {
      proposed->models[di].a[new_iy] = 
	rjmcmc_random_choose_double(data->ymin, data->ymax, random);
      
      prob *= (data->ymax - data->ymin);
    }      
    
  }

  *birth_prob = prob;

  return 1;
}

int 
part1d_zero_propose_death(const part1d_zero_t *current,
			  part1d_zero_t *proposed,
			  const dataset1d_t **datasets,
			  int ndatasets,
			  rjmcmc_uniform_rand_t random,
			  rjmcmc_normal_rand_t normal,
			  double *death_prob)
{
  int del_iy;
  double deleted_pos;
  int old_iy;
  int new_iy;
  int next_iy;
  int di;
  
  double left_x;
  double right_x;
  int xi;
  int xj;
  int n;
  
  double phi;
  double prob;

  const dataset1d_t *data;

  part1d_zero_clone(current, proposed);

  if (proposed->npartitions <= 2 || proposed->npartitions <= proposed->min_partitions) {
    /* Can't remove any more points */
    /* rjmcmc_error( */
    /* 	    "part1d_zero_propose_death:" */
    /* 	    "too few partitions %d %d\n", */
    /* 	    proposed->npartitions, */
    /* 	    proposed->min_partitions); */

    return 0;
  }
  
  /* Remove one value at random, note that we can't remove the two endpoints
   * which are always located and iy indices of 0 and 1, hence the 2's in 
   * the random choice of index to remove
   */
  del_iy = rjmcmc_random_choose_int(2, proposed->npartitions - 1, random);
  deleted_pos = position_map1d_position_of_index(proposed->p, del_iy);

  old_iy = position_map1d_predecessor_of_index(proposed->p, del_iy);
  if (old_iy < 0) {
    rjmcmc_error("part1d_zero_proposed_death: "
	    "failed to find predecessor of deleted point\n");
    return 0;
  }

  /* 
   * Record the previous mean/variance value
   */
  prob = 1.0;
  for (di = 0; di < proposed->ndatasets; di ++) {
    data = datasets[di];

    /*
     * Deleted partition (c)
     */
    if (proposed->models[di].var[del_iy] > 0.0) {

      prob *= 
	rjmcmc_gaussian_probability(proposed->models[di].a[del_iy] - 
				    proposed->models[di].mean[del_iy],
				    sqrt(proposed->models[di].var[del_iy]));

    } else {
      prob /= (data->ymax - data->ymin);
    }

    /*
     * Predecessor of deleted partition (b)
     */
    if (proposed->models[di].var[old_iy] > 0.0) {

      prob *=
	rjmcmc_gaussian_probability(proposed->models[di].a[old_iy] - 
				    proposed->models[di].mean[old_iy],
				    sqrt(proposed->models[di].var[old_iy]));
    } else {
      prob /= (data->ymax - data->ymin);
    }
  }

  if (position_map1d_delete(proposed->p, deleted_pos, del_iy) < 0) {
    rjmcmc_error("part1d_zero_propose_death: failed to delete position\n");
    return 0;
  }

  models_delete(proposed->max_partitions,
		proposed->ndatasets,
		del_iy,
		proposed->npartitions,
		proposed->models);
  proposed->npartitions --;

  new_iy = position_map1d_predecessor_of_point(proposed->p, deleted_pos);
  if (new_iy < 0) {
    rjmcmc_error("part1d_zero_propose_death: failed to find predecessor\n");
    return 0;
  }

  next_iy = position_map1d_successor_of_point(proposed->p, deleted_pos);
  if (new_iy < 0) {
    rjmcmc_error("part1d_zero_propose_death: failed to find successor\n");
    return 0;
  }

  /*
   * Now update the value in the previous partition (which now encompasses the deleted
   * partition).
   */

  left_x = position_map1d_position_of_index(proposed->p, new_iy);
  right_x = position_map1d_position_of_index(proposed->p, next_iy);

  for (di = 0; di < ndatasets; di ++) {
    
    data = datasets[di];

    n = dataset1d_range(data, left_x, right_x, &xi, &xj);
    if (n >= 2) {
      dataset1d_mean_variance(data, xi, xj, 
			      &(proposed->models[di].mean[new_iy]),
			      &(proposed->models[di].var[new_iy]));

      phi = sqrt(proposed->models[di].var[new_iy]) *
	normal();
      
      proposed->models[di].a[new_iy] = 
	proposed->models[di].mean[new_iy] + phi;
      
      prob /= 
	rjmcmc_gaussian_probability(phi, 
				    sqrt(proposed->models[di].var[new_iy]));

    } else {
      
      proposed->models[di].mean[new_iy] = 0.0;
      proposed->models[di].var[new_iy] = 0.0;

      proposed->models[di].a[new_iy] = 
	rjmcmc_random_choose_double(data->ymin,
				    data->ymax,
				    random);

      prob *= (data->ymax - data->ymin);
    }
  }

  *death_prob = prob;

  return 1;
}

int 
part1d_zero_propose_value(const part1d_zero_t *current, 
			  part1d_zero_t *proposed,
			  const dataset1d_t **datasets,
			  int ndatasets,
			  rjmcmc_uniform_rand_t random,
			  rjmcmc_normal_rand_t normal,
			  double *value_prob)
{
  int di;
  int iy;

  const dataset1d_t *data;

  double phi;
  double phiprime;
  double std;

  part1d_zero_clone(current, proposed);

  if (ndatasets == 1) {
    di = 0;
  } else {
    di = rjmcmc_random_choose_int(0,
				  ndatasets - 1,
				  random);
  }

  /*
   * We want to choose an index in the set (0, 2, 3, ..., npartitions - 1)
   * as partition index 1 has no bearing on the result.
   */
  if (proposed->npartitions == 2) {
    iy = 0;
  } else {
    iy = rjmcmc_random_choose_int(0, 
				  proposed->npartitions - 2,
				  random);
  }

  if (iy > 0) {
    iy ++;
  }

  data = datasets[di];

  if (proposed->models[di].var[iy] > 0.0) {

    std = sqrt(proposed->models[di].var[iy]);

    phi = proposed->models[di].a[iy] - 
      proposed->models[di].mean[iy];

    phiprime = normal() * std;

    proposed->models[di].a[iy] = 
      proposed->models[di].mean[iy] + phiprime;

    *value_prob = 
      rjmcmc_gaussian_probability(phi, std)/
      rjmcmc_gaussian_probability(phiprime, std);
  } else {

    proposed->models[di].a[iy] = 
      rjmcmc_random_choose_double(data->ymin,
				  data->ymax,
				  random);

    *value_prob = 1.0;
  }

  return 1;
}


int 
part1d_zero_propose_move(const part1d_zero_t *current,
				  part1d_zero_t *proposed,
				  const dataset1d_t **datasets,
				  int ndatasets,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *move_prob)
{
  int a_iy;
  int b_iy;
  int c_iy;

  double old_x;
  double new_x;

  double prob;

  int di;
  const dataset1d_t *data;

  double a_left_x;
  double a_right_x;
  double b_left_x;
  double b_right_x;
  double c_left_x;
  double c_right_x;

  int n;
  int xi;
  int xj;

  int tiy;

  double phi;

  /* The 2 end points can't be moved, if this is all there is then give up */
  if (current->npartitions <= 2) {
    return 0;
  }

  part1d_zero_clone(current, proposed);

  b_iy = rjmcmc_random_choose_int(2, 
				     proposed->npartitions - 1,
				     random);

  old_x = position_map1d_position_of_index(proposed->p, b_iy);
  new_x = old_x + normal() * proposed->pd;

  if (new_x <= proposed->xmin ||
      new_x >= proposed->xmax) {
    return 0;
  }

  if (position_map1d_move(proposed->p,
			old_x,
			new_x) < 0) {
    rjmcmc_error(
	    "part1d_zero_propose_move: "
	    "failed to move point\n");
    return 0;
  }

  c_iy = position_map1d_predecessor_of_index(proposed->p, b_iy);
  if (c_iy < 0) {
    rjmcmc_error("part1d_zero_propose_move: "
	    "failed to find predecessor\n");
    return 0;
  }

  a_iy = position_map1d_predecessor_of_index(current->p, b_iy);
  if (a_iy < 0) {
    rjmcmc_error("part1d_zero_proposed_move: "
	    "failed to find old predecessor.\n");
    return 0;
  }

  /*
   * Record the previous curve probabilities.
   */
  prob = 1.0;
  for (di = 0; di < ndatasets; di ++) {

    data = datasets[di];

    /* a */
    if (current->models[di].var[a_iy] > 0.0) {

      prob *= 
	rjmcmc_gaussian_probability(current->models[di].a[a_iy] - 
				    current->models[di].mean[a_iy],
				    sqrt(current->models[di].var[a_iy]));

    } else {

      prob /= (data->ymax - data->ymin);

    }

    /* b */
    if (current->models[di].var[b_iy] > 0.0) {

      prob *= 
	rjmcmc_gaussian_probability(current->models[di].a[b_iy] - 
				    current->models[di].mean[b_iy],
				    sqrt(current->models[di].var[b_iy]));

    } else {

      prob /= (data->ymax - data->ymin);

    }

    /* c (as required) */
    if (a_iy != c_iy) {

      if (current->models[di].var[c_iy] > 0.0) {
	
	prob *= 
	  rjmcmc_gaussian_probability(current->models[di].a[c_iy] - 
				      current->models[di].mean[c_iy],
				      sqrt(current->models[di].var[c_iy]));
	
      } else {
	
	prob /= (data->ymax - data->ymin);
	
      }
      
    }
  }

  /*
   * Update modified partitions as required 
   */

  a_left_x = position_map1d_position_of_index(proposed->p, a_iy);
  tiy = position_map1d_successor_of_index(proposed->p, a_iy);
  if (tiy < 0) {
    rjmcmc_error(
	    "part1d_zero_propose_move: "
	    "failed to find successor for a\n");
    return 0;
  }
  a_right_x = position_map1d_position_of_index(proposed->p, tiy);

  b_left_x = position_map1d_position_of_index(proposed->p, b_iy);
  tiy = position_map1d_successor_of_index(proposed->p, b_iy);
  if (tiy < 0) {
    rjmcmc_error(
	    "part1d_zero_propose_move: "
	    "failed to find successor for b\n");
    return 0;
  }
  b_right_x = position_map1d_position_of_index(proposed->p, tiy);

  if (c_iy != a_iy) {
    c_left_x = position_map1d_position_of_index(proposed->p, c_iy);
    c_right_x = b_left_x;
  }
  
  for (di = 0; di < ndatasets; di ++) {

    data = datasets[di];

    /* a (a_iy)*/
    n = dataset1d_range(data, a_left_x, a_right_x, &xi, &xj);
    if (n >= 2) {
      dataset1d_mean_variance(data, xi, xj,
			      &(proposed->models[di].mean[a_iy]),
			      &(proposed->models[di].var[a_iy]));

      phi = sqrt(proposed->models[di].var[a_iy]) *
	normal();
      
      proposed->models[di].a[a_iy] = 
	proposed->models[di].mean[a_iy] + phi;
      
      prob /= 
	rjmcmc_gaussian_probability(phi, 
				    sqrt(proposed->models[di].var[a_iy]));
    } else {

      proposed->models[di].mean[a_iy] = 0.0;
      proposed->models[di].var[a_iy] = 0.0;

      proposed->models[di].a[a_iy] =
	rjmcmc_random_choose_double(data->ymin,
				    data->ymax,
				    random);

      prob *= (data->ymax - data->ymin);
    }

    /* b (b_iy) */
    n = dataset1d_range(data, b_left_x, a_right_x, &xi, &xj);
    if (n >= 2) {
      dataset1d_mean_variance(data, xi, xj,
			      &(proposed->models[di].mean[b_iy]),
			      &(proposed->models[di].var[b_iy]));

      phi = sqrt(proposed->models[di].var[b_iy]) *
	normal();
      
      proposed->models[di].a[b_iy] = 
	proposed->models[di].mean[b_iy] + phi;
      
      prob /= 
	rjmcmc_gaussian_probability(phi, 
				    sqrt(proposed->models[di].var[b_iy]));
    } else {

      proposed->models[di].mean[b_iy] = 0.0;
      proposed->models[di].var[b_iy] = 0.0;

      proposed->models[di].a[b_iy] =
	rjmcmc_random_choose_double(data->ymin,
				    data->ymax,
				    random);

      prob *= (data->ymax - data->ymin);
    }

    /* c (c_iy, as required) */
    if (a_iy != c_iy) {
      n = dataset1d_range(data, c_left_x, c_right_x, &xi, &xj);
      if (n >= 2) {
	dataset1d_mean_variance(data, xi, xj,
				&(proposed->models[di].mean[c_iy]),
				&(proposed->models[di].var[c_iy]));
	
	phi = sqrt(proposed->models[di].var[c_iy]) *
	  normal();
	
	proposed->models[di].a[c_iy] = 
	  proposed->models[di].mean[c_iy] + phi;
	
	prob /= 
	  rjmcmc_gaussian_probability(phi, 
				      sqrt(proposed->models[di].var[c_iy]));
      } else {
	
	proposed->models[di].mean[c_iy] = 0.0;
	proposed->models[di].var[c_iy] = 0.0;
	
	proposed->models[di].a[c_iy] =
	  rjmcmc_random_choose_double(data->ymin,
				      data->ymax,
				      random);
	
	prob *= (data->ymax - data->ymin);
      }
    }
  }

  *move_prob = prob;

  return 1;
}

int 
part1d_zero_propose_lambda(const part1d_zero_t *current,
			   part1d_zero_t *proposed,
			   const dataset1d_t **datasets,
			   int ndatasets,
			   rjmcmc_uniform_rand_t random,
			   rjmcmc_normal_rand_t normal,
			   double *lambda_prob)
{
  int di;
  double new_l;
  const dataset1d_t *data;

  if (ndatasets > 1) {
    di = rjmcmc_random_choose_int(0, ndatasets - 1, random);
  } else {
    di = 0;
  }
  data = datasets[di];

  part1d_zero_clone(current, proposed);
  new_l = proposed->models[di].lambda + normal() * data->lambdastd;

  if (new_l < data->lambdamin ||
      new_l > data->lambdamax) {
    return 0;
  }

  *lambda_prob = pow(proposed->models[di].lambda/new_l,
		     data->npoints);
  proposed->models[di].lambda = new_l;
  return 1;
}

/*
 * Internal methods
 */

static model_t *models_create(int max_part,
			      int ndatasets)
{
  model_t *m;
  int di;

  m = (model_t*)malloc(sizeof(model_t) * ndatasets);
  if (m == NULL) {
    return NULL;
  }
  
  for (di = 0; di < ndatasets; di ++) {
    m[di].a = rjmcmc_create_array_1d(max_part);
    if (m[di].a == NULL) {
      return NULL;
    }

    m[di].mean = rjmcmc_create_array_1d(max_part);
    if (m[di].mean == NULL) {
      return NULL;
    }

    m[di].var = rjmcmc_create_array_1d(max_part);
    if (m[di].var == NULL) {
      return NULL;
    }
  }

  return m;
}

static void models_destroy(int max_part,
			   int ndatasets,
			   model_t *m)
{
  int di;

  if (m != NULL) {
    
    for (di = 0; di < ndatasets; di ++) {
      rjmcmc_destroy_array_1d(m[di].a);
      rjmcmc_destroy_array_1d(m[di].mean);
      rjmcmc_destroy_array_1d(m[di].var);
    }

    free(m);
  }
}

static void models_clone(int max_partitions,
			 int ndatasets,
			 const model_t *src,
			 model_t *dst)
{
  int di;
  int pi;

  RJMCMC_NULLCHECKVOID(src, "models_clone: null src\n");
  RJMCMC_NULLCHECKVOID(src, "models_clone: null dst\n");

  for (di = 0; di < ndatasets; di ++) {
    for (pi = 0; pi < max_partitions; pi ++) {
      dst[di].a[pi] = src[di].a[pi];
      dst[di].mean[pi] = src[di].mean[pi];
      dst[di].var[pi] = src[di].var[pi];
    }

    dst[di].lambda = src[di].lambda;
  }
       

}

static void models_delete(int max_partitions,
			  int ndatasets,
			  int del_iy,
			  int npart,
			  model_t *m)
{
  int di;
  int pi;

  for (di = 0; di < ndatasets; di ++) {
    for (pi = del_iy + 1; pi < npart; pi ++) {
      m[di].a[pi - 1] = m[di].a[pi];
      m[di].mean[pi - 1] = m[di].mean[pi];
      m[di].var[pi - 1] = m[di].var[pi];
    }
  }
}

int
part1d_zero_partitions(const part1d_zero_t *current)
{
  return current->npartitions;
}

double
part1d_zero_lambda(const part1d_zero_t *current,
			int di)
{
  return current->models[di].lambda;
}

double
part1d_zero_partition_position(const part1d_zero_t *current,
			       int pi)
{
  return position_map1d_position_of_index(current->p, pi);
}

