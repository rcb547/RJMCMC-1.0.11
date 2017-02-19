
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "rjmcmc/part1d_natural_rj.h"

#include "rjmcmc/position_map1d.h"
#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_defines.h"

struct _model {
  
  double *a;
  double lambda;

};

typedef struct _model model_t;
  
struct _part1d_natural_rj {

  /*
   * Constant parameters
   */
  int min_partitions;
  int max_partitions;
  int ndatasets;

  double xmin;
  double xmax;

  double pv;
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

static double value_at(const part1d_natural_rj_t *current,
		       int di,
		       double x);

static void value_gradient_at(const part1d_natural_rj_t *current,
			      int di,
			      double x,
			      double *v,
			      double *dv);

part1d_natural_rj_t *
part1d_natural_rj_create(int min_partitions,
			 int max_partitions,
			 int ndatasets,
			 double xmin, 
			 double xmax,
			 double pv,
			 double pd)
{
  part1d_natural_rj_t *r;
  
  r = (part1d_natural_rj_t*)malloc(sizeof(part1d_natural_rj_t));
  if (r == NULL) {
    return NULL;
  }

  r->min_partitions = min_partitions;
  r->max_partitions = max_partitions;
  r->ndatasets = ndatasets;

  r->xmin = xmin;
  r->xmax = xmax;
  r->pv = pv;
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
part1d_natural_rj_destroy(part1d_natural_rj_t *p)
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
part1d_natural_rj_clone(const part1d_natural_rj_t *src,
			part1d_natural_rj_t *dst)
{
  RJMCMC_NULLCHECKVOID(src, "part1d_natural_rj_clone: null src\n");
  RJMCMC_NULLCHECKVOID(src, "part1d_natural_rj_clone: null dst\n");

  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "part1d_natural_rj_clone: size mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->ndatasets != dst->ndatasets,
			    "part1d_natural_rj_clone: count mismatch\n");

  position_map1d_clone(src->p, dst->p);
  models_clone(src->max_partitions,
	       src->ndatasets,
	       src->models, dst->models);
  dst->npartitions = src->npartitions;
}

int 
part1d_natural_rj_initialize(part1d_natural_rj_t *p,
			     const dataset1d_t **datasets,
			     int ndatasets,
			     rjmcmc_uniform_rand_t random,
			     rjmcmc_normal_rand_t normal)
{
  int npart;
  int pi;
  int di;

  double x;

  const dataset1d_t *data;

  int i;

  npart = 2;
  
  for (pi = 2; pi < npart; pi ++) {
    x = p->xmin + (p->xmax - p->xmin)*random();
    position_map1d_insert(p->p, 
			x,
			pi);
  }
  
  p->npartitions = npart;
  for (i = 0; i < npart; i ++) {
    for (di = 0; di < ndatasets; di ++) {
      data = datasets[di];
      
      p->models[di].a[i] = data->ymin +
	(data->ymax - data->ymin)*random();
      
      p->models[di].lambda = 
	rjmcmc_random_choose_double(data->lambdamin,
				    data->lambdamax,
				    random);
    }
  }

  return 0;
}

double
part1d_natural_rj_value_at(const part1d_natural_rj_t *current,
			   double x)
{
  return value_at(current, 0, x);
}

static double
value_at(const part1d_natural_rj_t *current,
	 int di,
	 double x)
{
  double v, dv;
  value_gradient_at(current, di, x, &v, &dv);
  return v;
}

static void value_gradient_at(const part1d_natural_rj_t *current,
			      int di,
			      double x,
			      double *v,
			      double *dv)
{
  int iyl;
  int iyr;

  double xl;
  double xr;
  double yl;
  double yr;

  iyl = position_map1d_predecessor_of_point(current->p, x);
  if (iyl < 0) {
    *v = -DBL_MAX;
    *dv = -DBL_MAX;
    return;
  }

  xl = position_map1d_position_of_index(current->p, iyl);

  if (iyl == 1) {
    iyr = iyl;
    xr = xl;

    iyl = position_map1d_predecessor_of_index(current->p, iyl);
    xl = position_map1d_position_of_index(current->p, iyl);

    if (iyl < 0 || iyl == 1) {
      *v = -DBL_MAX;
      *dv = -DBL_MAX;
      return;
    }

  } else {
    iyr = position_map1d_next_index(current->p, xl);

    if (iyl < 0) {
      *v = -DBL_MAX;
      *dv = -DBL_MAX;
      return;
    }

    xr = position_map1d_position_of_index(current->p, iyr);
  }

  yl = current->models[di].a[iyl];
  yr = current->models[di].a[iyr];

  *v = yl + (yr - yl) * (x - xl)/(xr - xl);
  *dv = (yr - yl)/(xr - xl);
}

int
part1d_natural_rj_evaluate(const part1d_natural_rj_t *current,
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

int
part1d_natural_rj_evaluate_gradient(const part1d_natural_rj_t *current,
				    int di,
				    double xmin,
				    double xmax,
				    int nsamples,
				    double *samples)
{
  int i;

  double xi;
  double v;

  for (i = 0; i < nsamples; i ++) {
    xi = xmin + (double)i * (xmax - xmin)/(double)(nsamples - 1);
    value_gradient_at(current, di, xi, &v, &(samples[i]));
  }
   
  return 0;  
}
				    

double
part1d_natural_rj_misfit(part1d_natural_rj_t *p,
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
part1d_natural_rj_propose_birth(const part1d_natural_rj_t *current,
				part1d_natural_rj_t *proposed,
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

  const dataset1d_t *data;
  
  if (current->npartitions == current->max_partitions) {
    /* rjmcmc_error( */
    /* 	    "part1d_natural_rj_propose_birth: " */
    /* 	    "%d %d\n",  */
    /* 	    current->npartitions, */
    /* 	    current->max_partitions); */
    return 0;
  }
  
  part1d_natural_rj_clone(current, proposed);

  new_x = random() * (proposed->xmax - proposed->xmin) + proposed->xmin;
  new_iy = proposed->npartitions;

  if (position_map1d_insert(proposed->p, new_x, new_iy) < 0) {
    rjmcmc_error(
	    "part1d_natural_rj_propose_birth: "
	    "failed to add new point\n");
    return 0;
  }

  proposed->npartitions ++;

  new_x_right = position_map1d_next_position(proposed->p, new_x);
  if (new_x_right < new_x) {
    rjmcmc_error(
	    "part1d_natural_rj_propose_birth: "
	    "failed to find right extent of new point\n");
    return 0;
  }
    
  prob = 1.0;

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    proposed->models[di].a[new_iy] = 
      data->ymin + 
      (data->ymax - data->ymin) * random();
    
  }

  *birth_prob = prob;

  return -1;
}

int 
part1d_natural_rj_propose_death(const part1d_natural_rj_t *current,
				part1d_natural_rj_t *proposed,
				const dataset1d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal,
				double *death_prob)
{
  int del_iy;
  double deleted_pos;
  int new_iy;
  
  double prob;
  
  part1d_natural_rj_clone(current, proposed);

  if (proposed->npartitions <= 2 || proposed->npartitions <= proposed->min_partitions) {
    /* Can't remove any more points */
    /* rjmcmc_error( */
    /* 	    "part1d_natural_rj_propose_death:" */
    /* 	    "too few partitions %d %d\n", */
    /* 	    proposed->npartitions, */
    /* 	    proposed->min_partitions); */

    return 0;
  }
  
  /* Remove one value at random, note that we can't remove the two endpoints
   * which are always located and iy indices of 0 and 1, hence the 2's in 
   * the random choice of index to remove
   */
  del_iy = 2 + (int)(random() * (double)(proposed->npartitions - 2));
  deleted_pos = position_map1d_position_of_index(proposed->p, del_iy);
  
  if (position_map1d_delete(proposed->p, deleted_pos, del_iy) < 0) {
    rjmcmc_error("part1d_natural_rj_propose_death: failed to delete position\n");
    return 0;
  }

  new_iy = position_map1d_predecessor_of_point(proposed->p, deleted_pos);
  if (new_iy < 0) {
    rjmcmc_error("part1d_natural_rj_propose_death: failed to find predecessor\n");
    return 0;
  }

  models_delete(proposed->max_partitions,
		proposed->ndatasets,
		del_iy,
		proposed->npartitions,
		proposed->models);
  proposed->npartitions --;

  prob = 1.0;

  *death_prob = prob;

  return 1;
}

int 
part1d_natural_rj_propose_value(const part1d_natural_rj_t *current, 
				part1d_natural_rj_t *proposed,
				const dataset1d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal,
				double *value_prob)
{
  int di;
  int iy;

  const dataset1d_t *data;

  part1d_natural_rj_clone(current, proposed);

  if (ndatasets == 1) {
    di = 0;
  } else {
    di = rjmcmc_random_choose_int(0,
				  ndatasets - 1,
				  random);
  }

  iy = rjmcmc_random_choose_int(0, 
				proposed->npartitions - 1,
				random);

  data = datasets[di];

  proposed->models[di].a[iy] += proposed->pv * normal();

  if (proposed->models[di].a[iy] > data->ymax ||
      proposed->models[di].a[iy] < data->ymin) {
    return 0;
  }

  return 1;
}

int 
part1d_natural_rj_propose_lambda(const part1d_natural_rj_t *current,
				part1d_natural_rj_t *proposed,
				const dataset1d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal,
				double *lambda_prob)
{
  int di;
  double new_l;
  const dataset1d_t *data;

  di = (int)((double)ndatasets * random());
  data = datasets[di];

  part1d_natural_rj_clone(current, proposed);
  new_l = proposed->models[di].lambda + normal() * data->lambdastd;

  if (new_l < data->lambdamin ||
      new_l > data->lambdamax) {
    return 0;
  }

  *lambda_prob = pow(proposed->models[di].lambda/new_l,
		     data->npoints);
  proposed->models[di].lambda = new_l;
  return -1;
}

int 
part1d_natural_rj_propose_move(const part1d_natural_rj_t *current,
				  part1d_natural_rj_t *proposed,
				  const dataset1d_t **datasets,
				  int ndatasets,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *move_prob)
{
  int move_iy;

  double old_x;
  double new_x;

  /* The 2 end points can't be moved, if this is all there is then give up */
  if (current->npartitions <= 2) {
    return 0;
  }

  part1d_natural_rj_clone(current, proposed);

  move_iy = rjmcmc_random_choose_int(2, 
				     proposed->npartitions - 1,
				     random);

  old_x = position_map1d_position_of_index(proposed->p, move_iy);
  new_x = old_x + normal() * proposed->pd;

  if (new_x <= proposed->xmin ||
      new_x >= proposed->xmax) {
    return 0;
  }

  if (position_map1d_move(proposed->p,
			old_x,
			new_x) < 0) {
    rjmcmc_error("part1d_natural_rj_propose_move: "
		 "failed to move point\n");
    return 0;
  }

  *move_prob = 1.0;

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
  RJMCMC_NULLCHECKVOID(dst, "models_clone: null dst\n");

  for (di = 0; di < ndatasets; di ++) {
    for (pi = 0; pi < max_partitions; pi ++) {
      dst[di].a[pi] = src[di].a[pi];
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
    }
  }
}

int
part1d_natural_rj_partitions(const part1d_natural_rj_t *current)
{
  return current->npartitions - 1;
}

double
part1d_natural_rj_lambda(const part1d_natural_rj_t *current,
			 int di)
{
  return current->models[di].lambda;
}

double
part1d_natural_rj_partition_position(const part1d_natural_rj_t *current,
				     int pi)
{
  return position_map1d_position_of_index(current->p, pi);
}

