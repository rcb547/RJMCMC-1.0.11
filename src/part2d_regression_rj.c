
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/part2d_regression_rj.h>

#include <rjmcmc/position_map2d.h>
#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_defines.h>

struct _model {

  double *a;

  double lambda;

};
typedef struct _model model_t;

struct _part2d_regression_rj {

  /*
   * Constant parameters 
   */
  int min_partitions;
  int max_partitions;
  int ndatasets;
  double xmin;
  double xmax;
  double ymin;
  double ymax;

  double pv;
  double pdx;
  double pdy;
  double pdxy;
  
  /*
   * Varying parameters
   */
  int npartitions;

  /*
   * Coordinates of the partition centres
   */
  position_map2d_t *p;

  /*
   * Models for each dataset
   */
  model_t *models;
};

static model_t *models_create(int max_partitions, int ndatasets);

static void models_destroy(int max_partitions, int ndatasets, model_t *models);

static void models_clone(int max_partitions,
			 int ndatasets,
			 const model_t *src,
			 model_t *dst);

static void models_delete(int max_partitions,
			  int ndatasets,
			  int del_iy,
			  int npart,
			  model_t *models);

static double value_at(const part2d_regression_rj_t *current,
		       int di,
		       double x,
		       double y);

part2d_regression_rj_t *
part2d_regression_rj_create(int min_part,
			    int max_part,
			    int ndatasets,
			    double xmin, 
			    double xmax,
			    double ymin,
			    double ymax,
			    double pv,
			    double pdx,
			    double pdy,
			    double pdxy)
{
  part2d_regression_rj_t *r;

  r = (part2d_regression_rj_t*)malloc(sizeof(part2d_regression_rj_t));
  if (r == NULL) {
    return NULL;
  }

  r->min_partitions = min_part;
  if (r->min_partitions < 1) {
    r->min_partitions = 1;
  }
  r->max_partitions = max_part;
  r->ndatasets = ndatasets;

  r->xmin = xmin;
  r->xmax = xmax;
  r->ymin = ymin;
  r->ymax = ymax;
  
  r->pv = pv;
  r->pdx = pdx;
  r->pdy = pdy;
  r->pdxy = pdxy;

  r->npartitions = 0;

  r->p = position_map2d_create(max_part + 4,
			       r->xmin,
			       r->xmax,
			       r->ymin,
			       r->ymax);
  if (r->p == NULL) {
    return NULL;
  }

  r->models = models_create(max_part, ndatasets);
  if (r->models == NULL) {
    return NULL;
  }

  return r;
}

void
part2d_regression_rj_destroy(part2d_regression_rj_t *p)
{
  if (p != NULL) {
    
    position_map2d_destroy(p->p);

    models_destroy(p->max_partitions,
		   p->ndatasets,
		   p->models);

    free(p);
  }
}

void
part2d_regression_rj_clone(const part2d_regression_rj_t *src,
			   part2d_regression_rj_t *dst)
{
  RJMCMC_NULLCHECKVOID(src, "part2d_regression_rj_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "part2d_regression_rj_clone: null dst\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "part2d_regression_rj_clone: size mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->ndatasets != dst->ndatasets,
			    "part2d_regression_rj_clone: count mismatch\n");

  dst->npartitions = src->npartitions;
  
  position_map2d_clone(src->p, dst->p);

  models_clone(src->max_partitions,
	       src->ndatasets,
	       src->models,
	       dst->models);
}

double
part2d_regression_rj_misfit(part2d_regression_rj_t *p,
			    const dataset2d_t **datasets,
			    int ndatasets)
{
  int di;
  double sum;
  double dsum;

  double z;
  double dz;
  double n;
  double l2;

  int i;

  const dataset2d_t *data;


  sum = 0.0;

  for (di = 0; di < ndatasets; di ++) {

    dsum = 0.0;
    data = datasets[di];

    l2 = 1.0;
    if (data->lambdastd > 0.0) {
      l2 = p->models[di].lambda * p->models[di].lambda;
    }

    for (i = 0; i < data->npoints; i ++) {

      z = value_at(p, di, data->points[i].x, data->points[i].y);
      dz = z - data->points[i].z;
      
      n = data->points[i].n;
      dsum += (dz*dz)/(2.0*n*n*l2);
    }

    sum += dsum;
  }

  return sum;
}


int
part2d_regression_rj_evaluate(part2d_regression_rj_t *current,
			      int di,
			      double xmin,
			      double xmax,
			      int xsamples,
			      double ymin,
			      double ymax,
			      int ysamples,
			      double **z)
{
  int i;
  int j;
  double x;
  double y;

  /* models_dump(current, di); */

  for (i = 0; i < xsamples; i ++) {
    x = xmin + (double)i * (xmax - xmin)/(double)(xsamples - 1);
    for (j = 0; j < ysamples; j ++) {
      y = ymin + (double)j * (ymax - ymin)/(double)(ysamples - 1);
      z[i][j] = value_at(current, di, x, y);
    }
  }
  
  return 0;
}

int 
part2d_regression_rj_initialize(part2d_regression_rj_t *p,
				const dataset2d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal)
{
  int npart;
  int pi;
  int di;

  double x;
  double y;

  bbox2d_t bound;

  const dataset2d_t *data;

  npart = 1;
  
  for (pi = 0; pi < npart; pi ++) {

    x = rjmcmc_random_choose_double(p->xmin, p->xmax, random);
    y = rjmcmc_random_choose_double(p->ymin, p->ymax, random);


    if (position_map2d_insert(p->p, x, y, &bound) < 0) {
      rjmcmc_error("part2d_regression_rj_initialize: failed to insert point\n");
      return -1;
    }

    for (di = 0; di < ndatasets; di ++) {

      data = datasets[di];

      p->models[di].a[pi] = rjmcmc_random_choose_double(data->zmin,
      							data->zmax,
      							random);
    }
  }

  for (di = 0; di < ndatasets; di ++) {
    data = datasets[di];

    if (data->lambdastd > 0.0) {
      p->models[di].lambda = 
	rjmcmc_random_choose_double(data->lambdamin,
				    data->lambdamax,
				    random);
    } else {
      p->models[di].lambda = 0.0;
    }
  }

  p->npartitions = npart;

  return 0;
}

int
part2d_regression_rj_propose_birth(const part2d_regression_rj_t *current,
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *birth_prob)
{
  int di;
  int iy;

  double x;
  double y;

  const dataset2d_t *data;

  double prob;
  double vo;
  double dv;

  bbox2d_t bound;

  part2d_regression_rj_clone(current, proposed);

  if ((proposed->npartitions + 1) == proposed->max_partitions) {
    return 0;
  }

  iy = proposed->npartitions;

  x = rjmcmc_random_choose_double(proposed->xmin, proposed->xmax,
				  random);
  y = rjmcmc_random_choose_double(proposed->ymin, proposed->ymax,
				  random);

  if (position_map2d_insert(proposed->p, x, y, &bound) < 0) {
    rjmcmc_error(
	    "part2d_regression_rj_propose_birth: failed to insert point\n");
    return 0;
  }

  prob = 1.0;
  for (di = 0; di < ndatasets; di ++) {

    data = datasets[di];

    dv = current->pv * normal();
    vo = value_at(current, di, x, y);

    proposed->models[di].a[iy] = vo + dv;
    if (proposed->models[di].a[iy] < data->zmin ||
	proposed->models[di].a[iy] > data->zmax) {
      return 0;
    }

    prob *= rjmcmc_gaussian_probability(dv, current->pv);
  }

  proposed->npartitions ++;
  *birth_prob = prob;
  return -1;
}

int 
part2d_regression_rj_propose_death(const part2d_regression_rj_t *current,
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *death_prob)
{
  int iy;

  double dx;
  double dy;

  int di;
  double dz;
  double prob;
  
  const dataset2d_t *data;
  
  bbox2d_t bound;

  part2d_regression_rj_clone(current, proposed);

  if (proposed->npartitions == proposed->min_partitions ||
      proposed->npartitions <= 1) {
    return 0;
  }
  
  iy = rjmcmc_random_choose_int(0, proposed->npartitions - 1,
				random);
  
  if (position_map2d_position_of_index(proposed->p, iy + 4, &dx, &dy) < 0) {
    rjmcmc_error(
	    "part2d_regression_rj_propose_death: "
	    "failed to find deleted point\n");
    return 0;
  }
	    
  position_map2d_delete(proposed->p, iy + 4, &bound);

  models_delete(proposed->max_partitions,
		proposed->ndatasets,
		iy,
		proposed->npartitions,
		proposed->models);
  proposed->npartitions --;

  prob = 1.0;

  for (di = 0; di < proposed->ndatasets; di ++) {
    data = datasets[di];

    dz = value_at(current, di, dx, dy) -
      value_at(proposed, di, dx, dy);

    prob *= rjmcmc_gaussian_probability(dz, current->pv);
	
  }

  *death_prob = prob;
  return -1;
}

int 
part2d_regression_rj_propose_value(const part2d_regression_rj_t *current, 
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *value_prob)
{
  int di;
  int pi;

  double new_v;

  const dataset2d_t *data;

  part2d_regression_rj_clone(current, proposed);

  if (ndatasets == 1) {
    di = 0;
  } else {
    di = rjmcmc_random_choose_int(0, ndatasets - 1, random);
  }

  data = datasets[di];

  pi = rjmcmc_random_choose_int(0, proposed->npartitions - 1, random);

  new_v = proposed->models[di].a[pi] + (current->pv * normal());
  if (new_v < data->zmin || new_v > data->zmax) {
    return 0;
  }

  proposed->models[di].a[pi] = new_v;
  *value_prob = 1.0;
  return -1;
}

int 
part2d_regression_rj_propose_lambda(const part2d_regression_rj_t *current,
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *lambda_prob)
{
  int di;
  double new_lambda;

  const dataset2d_t *data;

  part2d_regression_rj_clone(current, proposed);

  if (ndatasets == 1) {
    di = 0;
  } else {
    di = rjmcmc_random_choose_int(0, ndatasets - 1, random);
  }

  data = datasets[di];

  RJMCMC_CONDITIONCHECKINT(data->lambdastd <= 0.0,
			   "part2d_regression_rj_propose_lambda: "
			   "invalid lambda standard deviation\n");


  new_lambda = proposed->models[di].lambda + 
    (data->lambdastd * normal());
  
  if (new_lambda < data->lambdamin ||
      new_lambda > data->lambdamax) {
    return 0;
  }

  proposed->models[di].lambda = new_lambda;
  *lambda_prob = (double)data->npoints * 
    (current->models[di].lambda/proposed->models[di].lambda);

  return -1;
}

int 
part2d_regression_rj_propose_move(const part2d_regression_rj_t *current,
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *move_prob)
{
  int iy;

  double mx;
  double my;

  double newx;
  double newy;

  bbox2d_t bound;
  
  part2d_regression_rj_clone(current, proposed);

  iy = rjmcmc_random_choose_int(0, proposed->npartitions - 1, random);
  
  if (position_map2d_position_of_index(proposed->p, iy + 4,
				       &mx, &my) < 0) {
    rjmcmc_error(
	    "part2d_regression_rj_propose_move:"
	    "failed to find move point\n");
    return -1;
  }

  newx = mx + (normal() * proposed->pdx);
  if (newx < proposed->xmin || newx > proposed->xmax) {
    return 0;
  }

  newy = my + (normal() * proposed->pdy);
  if (newy < proposed->ymin || newy > proposed->ymax) {
    return 0;
  }

  if (position_map2d_move(proposed->p, iy + 4, newx, newy, &bound) < 0) {
    rjmcmc_error("part2d_regression_rj_propose_move:"
		 "failed to move point %f %f -> %f %f (%d)\n",
		 mx, my,
		 newx, newy,
		 iy);
    return 0;
  }
	    

  *move_prob = 1.0;
  return -1;
}

int
part2d_regression_rj_partitions(const part2d_regression_rj_t *current)
{
  return current->npartitions;
}

double
part2d_regression_rj_lambda(const part2d_regression_rj_t *current)
{
  return current->models[0].lambda;
}

double
part2d_regression_rj_value_at(const part2d_regression_rj_t *current,
			      double x,
			      double y)
{
  return value_at(current, 0, x, y);
}

int 
part2d_regression_rj_partition_centre(const part2d_regression_rj_t *current,
				      int pi,
				      double *x,
				      double *y)
{
  if (pi >= 0 && pi < current->npartitions) {

    return position_map2d_position_of_index(current->p, pi, x, y);

  }

  return -1;
}



static model_t *models_create(int max_partitions, int ndatasets)
{
  model_t *m;
  int i;

  m = (model_t*)malloc(sizeof(model_t) * ndatasets);
  if (m == NULL) {
    return NULL;
  }

  for (i = 0; i < ndatasets; i ++) {
    m[i].a = rjmcmc_create_array_1d(max_partitions);
    if (m[i].a == NULL) {
      return NULL;
    }
  }

  return m;

}

static void models_destroy(int max_partitions, int ndatasets, model_t *models)
{
  int i;

  if (models != NULL) {

    for (i = 0; i < ndatasets; i ++) {
      rjmcmc_destroy_array_1d(models[i].a);
    }
    free(models);

  }
}

static void models_clone(int max_partitions,
			 int ndatasets,
			 const model_t *src,
			 model_t *dst)
{
  int i;
  int j;

  RJMCMC_NULLCHECKVOID(src, "models_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "models_clone: null dst\n");

  for (i = 0; i < ndatasets; i ++) {
    for (j = 0; j < max_partitions; j ++) {
      dst[i].a[j] = src[i].a[j];
    }

    dst[i].lambda = src[i].lambda;
  }
}

static void models_delete(int max_partitions,
			  int ndatasets,
			  int del_iy,
			  int npart,
			  model_t *models)
{
  int i;
  int j;

  for (i = 0; i < ndatasets; i ++) {
    for (j = del_iy + 1; j < npart; j ++) {
      models[i].a[j - 1] = models[i].a[j];
    }
  }
}

static double value_at(const part2d_regression_rj_t *current,
		       int di,
		       double x,
		       double y)
{
  int iy;
  double nx, ny;

  iy = position_map2d_nearest(current->p, 
			      x, y,
			      &nx, &ny,
			      0);
  if (iy < 0) {
    fprintf(stderr, "value_at: failed to find point: %f %f (%f %f %f %f)\n", x, y, current->xmin, current->xmax, current->ymin, current->ymax);
    return -DBL_MAX;
  }

  return current->models[di].a[iy - 4];
}

