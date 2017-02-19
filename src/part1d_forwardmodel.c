
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rjmcmc/part1d_forwardmodel.h"

#include "rjmcmc/position_map1d.h"
#include "rjmcmc/rjmcmc_util.h"

#include "rjmcmc/forwardmodel_util.h"

#include "rjmcmc/rjmcmc_defines.h"
#include "rjmcmc/rjmcmc_debug.h"

struct _model {
  double *local_parameter;
};

typedef struct _model model_t;
  
struct _part1d_forwardmodel {

  /*
   * Constant parameters
   */
  part1d_fm_type_t type;

  int min_partitions;
  int max_partitions;

  double xmin;
  double xmax;
  double pd;

  int nglobalparameters;
  double *global_parameter;

  /*
   * Varying parameter
   */
  int npartitions;

  /*
   * Coordinates of the partition boundaries
   */
  position_map1d_t *p;

  /*
   * Models (local parameters) for each partition
   */
  int nlocalparameters;
  model_t *models;

  /*
   * Global hierarchical parameters
   */
  int nhierarchicalparameters;
  double *hierarchical_parameters;

  /*
   * Model partition gradients (only for cubic partitions)
   */
  model_t *model_gradients;

  double *local_scratch;
};

static model_t *models_create(int max_part,
			      int nlocalparameters);
static void models_destroy(int max_part,
			   int nlocalparameters,
			   model_t *m);
static void models_clone(int max_partitions,
			 int nlocalparameters,
			 const model_t *src,
			 model_t *dst);

static void models_delete(int max_partitions,
			  int nlocalparameters,
			  int del_iy,
			  int npart,
			  model_t *m);

part1d_forwardmodel_t *
part1d_forwardmodel_create(part1d_fm_type_t type,
			   int min_partitions,
			   int max_partitions,
			   double xmin, 
			   double xmax,
			   double pd,
			   int nglobalparameters,
			   int nlocalparameters,
			   int nhierarchicalparameters)
{
  part1d_forwardmodel_t *r;
  
  r = (part1d_forwardmodel_t*)malloc(sizeof(part1d_forwardmodel_t));
  if (r == NULL) {
    rjmcmc_error("part1d_forwardmodel_create: failed to allocate memory\n");
    return NULL;
  }

  r->type = type;

  if (min_partitions < 2) {
    min_partitions = 2;
  }
  r->min_partitions = min_partitions;
  r->max_partitions = max_partitions;

  r->xmin = xmin;
  r->xmax = xmax;
  r->pd = pd;

  r->npartitions = 0;

  r->p = position_map1d_create(max_partitions, xmin, xmax);
  if (r->p == NULL) {
    rjmcmc_error("part1d_forwardmodel_create: failed to create position map\n");
    return NULL;
  }

  r->nglobalparameters = nglobalparameters;
  r->global_parameter = NULL;
  if (nglobalparameters > 0) {
    r->global_parameter = rjmcmc_create_array_1d(nglobalparameters);
    if (r->global_parameter == NULL) {
      return NULL;
    }
  }

  r->nlocalparameters = nlocalparameters;
  r->models = models_create(max_partitions,
			    nlocalparameters);
  if (r->models == NULL) {
    return NULL;
  }

  r->nhierarchicalparameters = nhierarchicalparameters;
  r->hierarchical_parameters = NULL;
  if (nhierarchicalparameters > 0) {
    r->hierarchical_parameters = rjmcmc_create_array_1d(nhierarchicalparameters);
    if (r->hierarchical_parameters == NULL) {
      return NULL;
    }
  }

  r->model_gradients = NULL;
  if (r->type == PART1D_FM_ZERO_CUBIC) {
    r->model_gradients = models_create(max_partitions,
				       nlocalparameters);
    if (r->model_gradients == NULL) {
      return NULL;
    }
  }

  r->local_scratch = rjmcmc_create_array_1d(nlocalparameters);
  if (r->local_scratch == NULL) {
    return NULL;
  }

  return r;
}

void
part1d_forwardmodel_destroy(part1d_forwardmodel_t *p)
{
  if (p != NULL) {

    position_map1d_destroy(p->p);

    models_destroy(p->max_partitions,
		   p->nlocalparameters,
		   p->models);

    if (p->model_gradients != NULL) {
      models_destroy(p->max_partitions,
		     p->nlocalparameters,
		     p->model_gradients);
    }

    rjmcmc_destroy_array_1d(p->hierarchical_parameters);
    rjmcmc_destroy_array_1d(p->global_parameter);

    rjmcmc_destroy_array_1d(p->local_scratch);

    free(p);

  }
}

void
part1d_forwardmodel_clone(const part1d_forwardmodel_t *src,
			  part1d_forwardmodel_t *dst)
{
  int i;

  RJMCMC_NULLCHECKVOID(src, "part1d_forwardmodel_clone: src null");
  RJMCMC_NULLCHECKVOID(dst, "part1d_forwardmodel_clone: dst null");

  RJMCMC_CONDITIONCHECKVOID(src->type != dst->type,
			    "part1d_forwardmodel_clone: type mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "part1d_forwardmodel_clone: size mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->nglobalparameters != dst->nglobalparameters,
			    "part1d_forwardmodel_clone: global mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->nlocalparameters != dst->nlocalparameters,
			    "part1d_forwardmodel_clone: local mismatch\n");

  position_map1d_clone(src->p, dst->p);
  models_clone(src->max_partitions,
	       src->nlocalparameters,
	       src->models, dst->models);
  dst->npartitions = src->npartitions;

  for (i = 0; i < src->nglobalparameters; i ++) {
    dst->global_parameter[i] = src->global_parameter[i];
  }

  for (i = 0; i < src->nhierarchicalparameters; i ++) {
    dst->hierarchical_parameters[i] = src->hierarchical_parameters[i];
  }
}

int 
part1d_forwardmodel_initialize_fixed(part1d_forwardmodel_t *p,
				     const double *global_parameters,
				     const double *hierarchical_parameters,
				     int npartitions,
				     const double *partition_x,
				     const double **local_parameters)
{
  int i;
  int li;

  p->npartitions = npartitions;

  for (i = 0; i < npartitions; i ++) {
    if (i >= 2) {
      position_map1d_insert(p->p, partition_x[i - 2], i);
    }
    
    for (li = 0; li < p->nlocalparameters; li ++) {
      p->models[i].local_parameter[li] = local_parameters[i][li];
    }
  }

  for (i = 0; i < p->nglobalparameters; i ++) {
    p->global_parameter[i] = global_parameters[i];
  }

  for (i = 0; i < p->nhierarchicalparameters; i ++) {
    p->hierarchical_parameters[i] = hierarchical_parameters[i];
  }
  
  return 0;
}


int 
part1d_forwardmodel_initialize(part1d_forwardmodel_t *p,
			       const forwardmodelparameter_t *global_parameters,
			       int nglobalparameters,
			       const forwardmodelparameter_t *local_parameters,
			       int nlocalparameters,
			       const forwardmodelparameter_t *hierarchical_parameters,
			       int nhierarchicalparameters,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal)
{
  int npart;
  int pi;
  int gi;
  int li;
  int si;

  double x;

  int i;
  RJMCMC_CONDITIONCHECKINT(nglobalparameters != p->nglobalparameters,
			    "part1d_forwardmodel_initialize: global mismatch\n");
  RJMCMC_CONDITIONCHECKINT(nlocalparameters != p->nlocalparameters,
			    "part1d_forwardmodel_initialize: local mismatch\n");
  RJMCMC_CONDITIONCHECKINT(nhierarchicalparameters != p->nhierarchicalparameters,
			    "part1d_forwardmodel_initialize: hiearchical mismatch\n");

  npart = p->min_partitions;
  
  for (gi = 0; gi < nglobalparameters; gi++) {
    p->global_parameter[gi] = 
      rjmcmc_random_choose_double(global_parameters[gi].fmin,
				  global_parameters[gi].fmax,
				  random);
  }

  for (pi = 2; pi < npart; pi ++) {
    x = rjmcmc_random_choose_double(p->xmin,
				    p->xmax,
				    random);
    position_map1d_insert(p->p, 
			x,
			pi);
  }

  p->npartitions = npart;
  for (i = 0; i < npart; i ++) {
    for (li = 0; li < nlocalparameters; li ++) {

      p->models[i].local_parameter[li] = 
	rjmcmc_random_choose_double(local_parameters[li].fmin,
				    local_parameters[li].fmax,
				    random);

    }
  }

  /*
   * Uniform initialization of hierarchical parameters
   */
  for (si = 0; si < nhierarchicalparameters; si ++) {
    p->hierarchical_parameters[si] = 
      rjmcmc_random_choose_double(hierarchical_parameters[si].fmin,
				  hierarchical_parameters[si].fmax,
				  random);
  }

  return 0;
}

int
part1d_forwardmodel_value_at(const part1d_forwardmodel_t *current,
			     double x,
			     int nlocalparameters,
			     double *localparameters)
{
  int iy;
  int riy;
  int i;

  double ix;
  double rix;
  double alpha;

  double a, b, d;

  RJMCMC_NULLCHECKINT(current, "part1d_forwardmodel_value_at: null state\n");
  RJMCMC_CONDITIONCHECKINT(nlocalparameters != current->nlocalparameters,
			   "part1d_forwardmodel_value_at: local mismatch\n");

  if (x < current->xmin || x > current->xmax) {
    rjmcmc_error("part1d_forwardmodel_value_at: out of range\n");
    return -1;
  }

  iy = position_map1d_predecessor_of_point(current->p, x);
  if (iy < 0) {
    rjmcmc_error("part1d_forwardmodel_value_at: failed to find predecessor (%f)\n", x);
    return -1;
  }

  switch (current->type) {
  case PART1D_FM_ZERO:
    if (iy == 1) {
      // Special case if we are at the end of the range
      iy = position_map1d_predecessor_of_index(current->p, iy);
      if (iy < 0) {
	rjmcmc_error("part1d_forwardmodel_value_at: failed to find predecessor of end point\n");
	return -1;
      }
    }

    for (i = 0; i < nlocalparameters; i ++) {
      localparameters[i] = current->models[iy].local_parameter[i];
    }
    break;

  case PART1D_FM_NATURAL:
    
    if (x >= current->xmax) {
      for (i = 0; i < nlocalparameters; i ++) {
	localparameters[i] = current->models[iy].local_parameter[i];
      }
    } else {
      riy = position_map1d_successor_of_point(current->p, x);
      if (riy < 0) {
	rjmcmc_error("part1d_forwardmodel_value_at: failed to find successor (%d %d %f)\n",
		     iy,
		     riy, 
		     x);
	return -1;
      }
      
      ix = position_map1d_position_of_index(current->p, iy);
      rix = position_map1d_position_of_index(current->p, riy);
      
      alpha = (x - ix)/(rix - ix);
      
      for (i = 0; i < nlocalparameters; i ++) {
	localparameters[i] = 
	  (1.0 - alpha) * current->models[iy].local_parameter[i] + 
	  alpha * current->models[riy].local_parameter[i];
      }
    }
    break;

  case PART1D_FM_ZERO_CUBIC:
    
    if (x >= current->xmax) {
      for (i = 0; i < nlocalparameters; i ++) {
	localparameters[i] = current->models[iy].local_parameter[i];
      }
    } else {
      riy = position_map1d_successor_of_point(current->p, x);
      if (riy < 0) {
	rjmcmc_error("part1d_forwardmodel_value_at: failed to find successor (%d %d %f)\n",
		     iy,
		     riy, 
		     x);
	return -1;
      }
      
      ix = position_map1d_position_of_index(current->p, iy);
      rix = position_map1d_position_of_index(current->p, riy);
      
      alpha = (x - ix)/(rix - ix);
      
      for (i = 0; i < nlocalparameters; i ++) {

	d = current->models[iy].local_parameter[i];
	b = 3.0 * (current->models[riy].local_parameter[i] -
		   d);
	a = -2.0 * b/3.0;

	localparameters[i] = 
	  a*alpha*alpha*alpha + b*alpha*alpha + d;
      }
    }
    break;

  default:
    rjmcmc_error("part1d_forwardmodel_value_at: invalid type\n");
    return -1;
  }

  return 0;
}

int 
part1d_forwardmodel_gradient_at(const part1d_forwardmodel_t *current,
				double x,
				int nlocalparameters,
				double *localparameters)
{
  return -1;
}

int
part1d_forwardmodel_propose_birth(const part1d_forwardmodel_t *current,
				  part1d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameters,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *birth_prob)
{
  double new_x;
  int new_iy;

  double prob;

  int li;

  double dv;

  if (current->npartitions == current->max_partitions) {
    /* rjmcmc_error( */
    /* 	    "part1d_forwardmodel_propose_birth: " */
    /* 	    "%d %d\n",  */
    /* 	    current->npartitions, */
    /* 	    current->max_partitions); */
    return 0;
  }
  
  part1d_forwardmodel_clone(current, proposed);

  new_x = rjmcmc_random_choose_double(proposed->xmin, proposed->xmax, random);
  new_iy = proposed->npartitions;

  if (position_map1d_insert(proposed->p, new_x, new_iy) < 0) {
    rjmcmc_error(
	    "part1d_forwardmodel_propose_birth: "
	    "failed to add new point %f %d\n", new_x, new_iy);
    return 0;
  }

  proposed->npartitions ++;

  prob = 1.0;

  if (part1d_forwardmodel_value_at(current, 
				   new_x, 
				   current->nlocalparameters, 
				   proposed->local_scratch) < 0) {
    rjmcmc_error("part1d_forwardmodel_propose_birth: "
		 "failed to find values at new position\n");
    return 0;
  }

  for (li = 0; li < proposed->nlocalparameters; li ++) {
    
    dv = local_parameters[li].fstd_bd * normal();
    proposed->models[new_iy].local_parameter[li] = 
      proposed->local_scratch[li] + dv;

    /* printf("  bc (%d): %f -> %f\n", */
    /* 	   li, */
    /* 	   current->models[old_iy].local_parameter[li], */
    /* 	   proposed->models[new_iy].local_parameter[li]); */

    if (proposed->models[new_iy].local_parameter[li] < local_parameters[li].fmin ||
	proposed->models[new_iy].local_parameter[li] > local_parameters[li].fmax) {
      return 0;
    }

    /* printf("  bp: (%d): %f %f\n", */
    /* 	   li, */
    /* 	   dv, */
    /* 	   local_parameters[li].fstd); */
    prob *= rjmcmc_gaussian_probability(dv, local_parameters[li].fstd_bd);
   
    /* printf(" b: %g %g %g\n", dv, local_parameters[li].fstd, rjmcmc_gaussian_probability(dv, local_parameters[li].fstd));  */

  }

  *birth_prob = prob;

  return 1;
}

int 
part1d_forwardmodel_propose_death(const part1d_forwardmodel_t *current,
				  part1d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameters,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *death_prob)
{
  int del_iy;
  double deleted_pos;

  double prob;

  int li;
  double dv;
  
  part1d_forwardmodel_clone(current, proposed);

  if (proposed->npartitions <= 2 || proposed->npartitions <= proposed->min_partitions) {
    /* Can't remove any more points */

    /* rjmcmc_error( */
    /* 	    "part1d_forwardmodel_propose_death:" */
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
  
  if (position_map1d_delete(proposed->p, deleted_pos, del_iy) < 0) {
    rjmcmc_error("part1d_forwardmodel_propose_death: failed to delete position\n");
    return 0;
  }

  models_delete(proposed->max_partitions,
		proposed->nlocalparameters,
		del_iy,
		proposed->npartitions,
		proposed->models);
  proposed->npartitions --;

  prob = 1.0;

  if (part1d_forwardmodel_value_at(current, 
				   deleted_pos, 
				   current->nlocalparameters, 
				   current->local_scratch) < 0) {
    rjmcmc_error("part1d_forwardmodel_propose_death: "
		 "failed to find values at new position\n");
    return 0;
  }

  if (part1d_forwardmodel_value_at(proposed, 
				   deleted_pos, 
				   proposed->nlocalparameters, 
				   proposed->local_scratch) < 0) {
    rjmcmc_error("part1d_forwardmodel_propose_death: "
		 "failed to find values at new position\n");
    return 0;
  }

  for (li = 0; li < proposed->nlocalparameters; li ++) {
    
    dv = proposed->local_scratch[li] - current->local_scratch[li];

    /* printf(" d: %g %g %g (%g %g)\n", dv, local_parameters[li].fstd, rjmcmc_gaussian_probability(dv, local_parameters[li].fstd), proposed->local_scratch[li], current->local_scratch[li]); */
    prob *= rjmcmc_gaussian_probability(dv, local_parameters[li].fstd_bd);
		
  }

  *death_prob = prob;

  return 1;
}

int 
part1d_forwardmodel_propose_local_value(const part1d_forwardmodel_t *current, 
					part1d_forwardmodel_t *proposed,
					int nglobalparameters,
					const forwardmodelparameter_t *global_parameters,
					int nlocalparameters,
					const forwardmodelparameter_t *local_parameters,
					rjmcmc_uniform_rand_t random,
					rjmcmc_normal_rand_t normal,
					double *value_prob,
					int *li)
{
  int iy;

  part1d_forwardmodel_clone(current, proposed);

  if (nlocalparameters == 1) {
    *li = 0;
  } else {
    *li = rjmcmc_random_choose_int(0,
				   nlocalparameters - 1,
				   random);
  }

  if (current->type == PART1D_FM_ZERO) {
    /* Since indices 0 and 1 are endpoints we don't want to perturb 1 */
    
    iy = rjmcmc_random_choose_int(0, 
				  proposed->npartitions - 2,
				  random);
    if (iy > 0) {
      iy ++; /* 0 or 2..npartitions - 1 */
    }
  } else {

    iy = rjmcmc_random_choose_int(0, 
				  proposed->npartitions - 1,
				  random);
  }

  proposed->models[iy].local_parameter[*li] += 
    local_parameters[*li].fstd_value * normal();

  if (proposed->models[iy].local_parameter[*li] < local_parameters[*li].fmin ||
      proposed->models[iy].local_parameter[*li] > local_parameters[*li].fmax) {
    return 0;
  }

  *value_prob = 0.0;
  return 1;
}

int 
part1d_forwardmodel_propose_global_value(const part1d_forwardmodel_t *current, 
					 part1d_forwardmodel_t *proposed,
					 int nglobalparameters,
					 const forwardmodelparameter_t *global_parameters,
					 int nlocalparameters,
					 const forwardmodelparameter_t *local_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *value_prob)
{
  int gi;

  if (nglobalparameters <= 0) {
    rjmcmc_error("part1d_forwardmodel_propose_global_value: "
		 "no global parameters.\n");
    return 0;
  }

  part1d_forwardmodel_clone(current, proposed);

  if (nglobalparameters == 1) {
    gi = 0;
  } else {
    gi = rjmcmc_random_choose_int(0,
				  nglobalparameters - 1,
				  random);
  }

  proposed->global_parameter[gi] += normal() * global_parameters[gi].fstd_value;

  if (proposed->global_parameter[gi] < global_parameters[gi].fmin ||
      proposed->global_parameter[gi] > global_parameters[gi].fmax) {
    return 0;
  }

  *value_prob = 0.0;
  return 1;
}

int 
part1d_forwardmodel_propose_move(const part1d_forwardmodel_t *current,
				 part1d_forwardmodel_t *proposed,
				 int nglobalparameters,
				 const forwardmodelparameter_t *global_parameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *local_parameters,
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

  part1d_forwardmodel_clone(current, proposed);

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
    rjmcmc_error(
	    "part1d_forwardmodel_propose_move: "
	    "failed to move point\n");
    return 0;
  }

  *move_prob = 1.0;

  return 1;
}

int
part1d_forwardmodel_propose_hierarchical(const part1d_forwardmodel_t *current,
					 part1d_forwardmodel_t *proposed,
					 int nhierarchicalparameters,
					 const forwardmodelparameter_t *hierarchical_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *hierarchical_prob)
{
  double dv;
  int hierarchicalindex;

  part1d_forwardmodel_clone(current, proposed);

  if (nhierarchicalparameters > 1) {
    hierarchicalindex = rjmcmc_random_choose_int(0, 
						 nhierarchicalparameters - 1, 
						 random);
  } else {
    hierarchicalindex = 0;
  }

  dv = normal() * hierarchical_parameters[hierarchicalindex].fstd_value;
  proposed->hierarchical_parameters[hierarchicalindex] += dv;

  if ((proposed->hierarchical_parameters[hierarchicalindex] < 
       hierarchical_parameters[hierarchicalindex].fmin) ||
      (proposed->hierarchical_parameters[hierarchicalindex] > 
       hierarchical_parameters[hierarchicalindex].fmax)) {

    return 0;
  }

  *hierarchical_prob = 0.0;

  return 1;
}

/*
 * Internal methods
 */

static model_t *models_create(int max_part,
			      int nlocalparameters)
{
  model_t *m;
  int iy;

  m = (model_t*)malloc(sizeof(model_t) * max_part);
  if (m == NULL) {
    return NULL;
  }
  
  for (iy = 0; iy < max_part; iy ++) {
    m[iy].local_parameter = rjmcmc_create_array_1d(nlocalparameters);
    if (m[iy].local_parameter == NULL) {
      return NULL;
    }
  }

  return m;
}

static void models_destroy(int max_part,
			   int nlocalparameters,
			   model_t *m)
{
  int iy;

  if (m != NULL) {
    
    for (iy = 0; iy < max_part; iy ++) {
      rjmcmc_destroy_array_1d(m[iy].local_parameter);
    }

    free(m);
  }
}

static void models_clone(int max_partitions,
			 int nlocalparameters,
			 const model_t *src,
			 model_t *dst)
{
  int iy;
  int li;

  RJMCMC_NULLCHECKVOID(src, "models_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "models_clone: null dst\n");

  for (iy = 0; iy < max_partitions; iy ++) {
    for (li = 0; li < nlocalparameters; li ++) {
      dst[iy].local_parameter[li] = src[iy].local_parameter[li];
    }

  }
       

}

static void models_delete(int max_partitions,
			  int nlocalparameters,
			  int del_iy,
			  int npart,
			  model_t *m)
{
  int iy;
  int li;

  for (iy = del_iy + 1; iy < max_partitions; iy ++) {
    for (li = 0; li < nlocalparameters; li ++) {
      m[iy - 1].local_parameter[li] = m[iy].local_parameter[li];
    }
  }
}

int
part1d_forwardmodel_partitions(const part1d_forwardmodel_t *current)
{
  return current->npartitions;
}

double
part1d_forwardmodel_partition_position(const part1d_forwardmodel_t *current,
				       int pi)
{
  return position_map1d_position_of_index(current->p, pi);
}

const double *
part1d_forwardmodel_global_parameters(const part1d_forwardmodel_t *current)
{
  return current->global_parameter;
}

const double *
part1d_forwardmodel_hierarchical_parameters(const part1d_forwardmodel_t *current)
{
  return current->hierarchical_parameters;
}

struct interval_data {
  const part1d_forwardmodel_t *current;

  int xi;
  int xsamples;
  const double *x;

  int nlocalparameters;
  double **local_parameters;
};

static int zero_interval(void *user_arg,
			 double xmin,
			 double xmax,
			 int iy,
			 int riy)
{
  struct interval_data *arg = (struct interval_data *)user_arg;
  int li;

  while ((arg->xi < arg->xsamples) && 
	 (arg->x[arg->xi] < xmax)) {

    for (li = 0; li < arg->nlocalparameters; li ++) {
      arg->local_parameters[li][arg->xi] = 
	arg->current->models[iy].local_parameter[li];
    }

    arg->xi ++;
  }

  if (riy == 1) {
    /* Last interval so fill the rest of the values with the 2nd last value */
    while (arg->xi < arg->xsamples) {

      for (li = 0; li < arg->nlocalparameters; li ++) {
	arg->local_parameters[li][arg->xi] = 
	  arg->current->models[iy].local_parameter[li];
      }

      arg->xi ++;
    }
  }

  return 0;
}

static int natural_interval(void *user_arg,
			    double xmin,
			    double xmax,
			    int iy,
			    int riy)
{
  struct interval_data *arg = (struct interval_data *)user_arg;
  int li;

  double alpha;

  /* printf("ni: %f %f %d %d\n", xmin, xmax, iy, riy); */

  if (riy < 0) {
    rjmcmc_error("natural_interval: invalid riy %d (%d) (%f %f)\n", riy, iy, xmin, xmax);
    return -1;
  }

  while ((arg->xi < arg->xsamples) && 
	 (arg->x[arg->xi] < xmax)) {

    alpha = (arg->x[arg->xi] - xmin)/(xmax - xmin);

    for (li = 0; li < arg->nlocalparameters; li ++) {
      arg->local_parameters[li][arg->xi] = 
	(1.0 - alpha) * arg->current->models[iy].local_parameter[li] +
	alpha * arg->current->models[riy].local_parameter[li];
    }

    arg->xi ++;
  }

  if (riy == 1) {
    /* Last interval so fill the rest of the values with the last value */
    while (arg->xi < arg->xsamples) {

      for (li = 0; li < arg->nlocalparameters; li ++) {
	arg->local_parameters[li][arg->xi] = 
	  arg->current->models[riy].local_parameter[li];
      }

      arg->xi ++;
    }
  }

  return 0;
}

static int zero_cubic_interval(void *user_arg,
			       double xmin,
			       double xmax,
			       int iy,
			       int riy)
{
  struct interval_data *arg = (struct interval_data *)user_arg;
  int li;

  double alpha;

  double a, b, d;

  /* printf("ni: %f %f %d %d\n", xmin, xmax, iy, riy); */

  if (riy < 0) {
    rjmcmc_error("zero_cubic_interval: invalid riy %d (%d) (%f %f)\n", riy, iy, xmin, xmax);
    return -1;
  }

  while ((arg->xi < arg->xsamples) && 
	 (arg->x[arg->xi] < xmax)) {

    alpha = (arg->x[arg->xi] - xmin)/(xmax - xmin);

    for (li = 0; li < arg->nlocalparameters; li ++) {

      d = arg->current->models[iy].local_parameter[li];
      b = 3.0 * (arg->current->models[riy].local_parameter[li] -
		 d);
      a = -2.0 * b/3.0;

      arg->local_parameters[li][arg->xi] = 
	a*alpha*alpha*alpha + b*alpha*alpha + d;
    }

    arg->xi ++;
  }

  if (riy == 1) {
    /* Last interval so fill the rest of the values with the last value */
    while (arg->xi < arg->xsamples) {

      for (li = 0; li < arg->nlocalparameters; li ++) {
	arg->local_parameters[li][arg->xi] = 
	  arg->current->models[riy].local_parameter[li];
      }

      arg->xi ++;
    }
  }

  return 0;
}

int
part1d_forwardmodel_evaluate_local_parameters(const part1d_forwardmodel_t *c,
					      int xsamples,
					      const double *x,
					      double **local_parameters)
{

  struct interval_data data;

  data.current = c;
  data.xi = 0;
  data.xsamples = xsamples;
  data.x = x;
  data.nlocalparameters = c->nlocalparameters;
  data.local_parameters = local_parameters;
  
  switch (c->type) {

  case PART1D_FM_ZERO:
    return position_map1d_traverse_intervals(c->p,
					   zero_interval,
					   &data);

  case PART1D_FM_NATURAL:
    return position_map1d_traverse_intervals(c->p,
					   natural_interval,
					   &data);

  case PART1D_FM_ZERO_CUBIC:
    return position_map1d_traverse_intervals(c->p,
					   zero_cubic_interval,
					   &data);
    

  default:
    rjmcmc_error("part1d_forwardmodel_evaluate_local_parameters: "
		 "invalid type\n");
    return -1;
  }
}

int 
part1d_forwardmodel_partition_fill_list(const part1d_forwardmodel_t *current,
					double *positions,
					int *npartitions)
{
  return position_map1d_fill_list(current->p, positions, npartitions);
}

int
part1d_forwardmodel_hierarchical_fill_list(const part1d_forwardmodel_t *current,
				    double *hierarchical,
				    int *nhierarchical)
{
  int i;

  for (i = 0; i < current->nhierarchicalparameters; i ++) {
    hierarchical[i] = current->hierarchical_parameters[i];
  }
  *nhierarchical = current->nhierarchicalparameters;

  return 0;
}

