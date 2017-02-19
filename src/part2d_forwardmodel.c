
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rjmcmc/part2d_forwardmodel.h"

#include "rjmcmc/position_map2d.h"
#include "rjmcmc/rjmcmc_util.h"

#include "rjmcmc/rjmcmc_defines.h"

struct _model {
  double *local_parameter;
};

typedef struct _model model_t;
  
struct _part2d_forwardmodel {

  /*
   * Constant parameters
   */
  part2d_fm_type_t type;

  int min_partitions;
  int max_partitions;
  int include_corners;

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;

  int nglobalparameters;
  double *global_parameter;

  /*
   * Global hierarchical parameters
   */
  int nhierarchicalparameters;
  double *hierarchical_parameter;

  /*
   * Varying parameter
   */
  int npartitions;

  /*
   * Coordinates of the partition boundaries
   */
  position_map2d_t *p;

  /*
   * Models (local parameters) for each partition
   */
  int nlocalparameters;
  model_t *models;

  /*
   * Model partition gradients (only for cubic partitions)
   */
  model_t *model_gradients;

  /*
   * Temporary storage for one set of localparameters
   */
  double *t_lp;

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

part2d_forwardmodel_t *
part2d_forwardmodel_create(part2d_fm_type_t type,
			   int min_partitions,
			   int max_partitions,
			   double xmin, 
			   double xmax,
			   double ymin,
			   double ymax,
			   double pdx,
			   double pdy,
			   int nglobalparameters,
			   int nlocalparameters,
			   int nhierarchicalparameters,
			   int include_corners)
{
  part2d_forwardmodel_t *r;
  
  r = (part2d_forwardmodel_t*)malloc(sizeof(part2d_forwardmodel_t));
  if (r == NULL) {
    return NULL;
  }

  r->type = type;

  if (min_partitions < 1) {
    min_partitions = 1;
  }

  if (min_partitions > 1) {
    r->min_partitions = min_partitions;
  } else {
    r->min_partitions = 1;
  }
  
  r->max_partitions = max_partitions;

  r->xmin = xmin;
  r->xmax = xmax;
  r->ymin = ymin;
  r->ymax = ymax;
  r->pdx = pdx;
  r->pdy = pdy;

  r->npartitions = 0;

  r->p = position_map2d_create(max_partitions,
			       r->xmin,
			       r->xmax,
			       r->ymin,
			       r->ymax);
  if (r->p == NULL) {
    return NULL;
  }

  r->nglobalparameters = nglobalparameters;
  r->global_parameter = NULL;
  if (r->nglobalparameters > 0) {
    r->global_parameter = rjmcmc_create_array_1d(nglobalparameters);
    if (r->global_parameter == NULL) {
      return NULL;
    }
  }

  r->nhierarchicalparameters = nhierarchicalparameters;
  r->hierarchical_parameter = NULL;
  if (r->nhierarchicalparameters > 0) {
    r->hierarchical_parameter = rjmcmc_create_array_1d(nhierarchicalparameters);
    if (r->hierarchical_parameter == NULL) {
      return NULL;
    }
  }

  r->nlocalparameters = nlocalparameters;
  r->models = models_create(max_partitions,
			    nlocalparameters);
  if (r->models == NULL) {
    return NULL;
  }

  r->model_gradients = NULL;
  if (r->type == PART2D_FM_MONOTONE_CUBIC) {
    r->model_gradients = models_create(max_partitions,
				       nlocalparameters);
    if (r->model_gradients == NULL) {
      return NULL;
    }
  }

  r->include_corners = include_corners;

  r->t_lp = rjmcmc_create_array_1d(nlocalparameters);
  if (r->t_lp == NULL) {
    return NULL;
  }

  return r;
}

void
part2d_forwardmodel_destroy(part2d_forwardmodel_t *p)
{
  if (p != NULL) {

    position_map2d_destroy(p->p);

    models_destroy(p->max_partitions,
		   p->nlocalparameters,
		   p->models);

    if (p->model_gradients != NULL) {
      models_destroy(p->max_partitions,
		     p->nlocalparameters,
		     p->model_gradients);
    }

    rjmcmc_destroy_array_1d(p->global_parameter);
    rjmcmc_destroy_array_1d(p->hierarchical_parameter);

    rjmcmc_destroy_array_1d(p->t_lp);
    free(p);

  }
}


static int write_int(FILE *fp, int i);
static int write_double(FILE *fp, double d);

static int read_int(FILE *fp, int *i);
static int read_double(FILE *fp, double *d);

int
part2d_forwardmodel_save(const part2d_forwardmodel_t *p,
			 const char *filename)
{
  FILE *fp;
  int i;
  int j;
  double x;
  double y;

  fp = fopen(filename, "w");
  RJMCMC_NULLCHECKINT(fp, "part2d_forwardmodel_save: failed to create file\n");

  RJMCMC_INTCHECKINT(write_int(fp, (int)(p->type)), "part2d_forwardmodel_save: failed to write type.\n");
  
  RJMCMC_INTCHECKINT(write_int(fp, p->min_partitions), 
		     "part2d_forwardmodel_save: failed to write min partitions.\n");
  RJMCMC_INTCHECKINT(write_int(fp, p->max_partitions), 
		     "part2d_forwardmodel_save: failed to write max partitions.\n");
  RJMCMC_INTCHECKINT(write_int(fp, p->include_corners), 
		     "part2d_forwardmodel_save: failed to write include corners.\n");

  RJMCMC_INTCHECKINT(write_double(fp, p->xmin),
		     "part2d_forwardmodel_save: failed to write xmin\n");
  RJMCMC_INTCHECKINT(write_double(fp, p->xmax),
		     "part2d_forwardmodel_save: failed to write xmax\n");
  RJMCMC_INTCHECKINT(write_double(fp, p->ymin),
		     "part2d_forwardmodel_save: failed to write ymin\n");
  RJMCMC_INTCHECKINT(write_double(fp, p->ymax),
		     "part2d_forwardmodel_save: failed to write ymax\n");

  RJMCMC_INTCHECKINT(write_double(fp, p->pdx),
		     "part2d_forwardmodel_save: failed to write pdx\n");
  RJMCMC_INTCHECKINT(write_double(fp, p->pdy),
		     "part2d_forwardmodel_save: failed to write pdy\n");

  RJMCMC_INTCHECKINT(write_int(fp, p->nglobalparameters),
		     "part2d_forwardmodel_save: failed to write nglobalparameters\n");
  RJMCMC_INTCHECKINT(write_int(fp, p->nlocalparameters),
		     "part2d_forwardmodel_save: failed to write nlocalparameters\n");
  RJMCMC_INTCHECKINT(write_int(fp, p->nhierarchicalparameters),
		     "part2d_forwardmodel_save: failed to write nhierarchicalparameters\n");

  /*
   * Global parameters
   */
  for (i = 0; i < p->nglobalparameters; i ++) {
    RJMCMC_INTCHECKINT(write_double(fp, p->global_parameter[i]),
		       "part2d_forwardmodel_save: failed to write global parameter\n");
  }

  /*
   * Hierarchical parameters 
   */
  for (i = 0; i < p->nhierarchicalparameters; i ++) {
    RJMCMC_INTCHECKINT(write_double(fp, p->hierarchical_parameter[i]),
		       "part2d_forwardmodel_save: failed to write hierarchical parameter\n");
  }

  /*
   * Point locations
   */
  RJMCMC_INTCHECKINT(write_int(fp, p->npartitions),
		     "part2d_forwardmodel_save: failed to write npartitions\n");
  for (i = 0; i < p->npartitions; i ++) {
    RJMCMC_INTCHECKINT(position_map2d_position_of_index(p->p, i + 4, &x, &y), 
		       "part2d_forwardmodel_save: failed to get position\n");

    RJMCMC_INTCHECKINT(write_double(fp, x),
		       "part2d_forwardmodel_save: failed to write position x\n");
    RJMCMC_INTCHECKINT(write_double(fp, y),
		       "part2d_forwardmodel_save: failed to write position y\n");
  }

  /*
   * Local parameters
   */
  for (j = 0; j < p->npartitions; j ++) {
    for (i = 0; i < p->nlocalparameters; i ++) {
      RJMCMC_INTCHECKINT(write_double(fp, p->models[j].local_parameter[i]),
			 "part2d_forwardmodel_save: failed to write local parameter\n");
    }
  }
  
  /*
   * Local gradients (Monotone cubic is not implemented yet).
   */
  if (p->type == PART2D_FM_MONOTONE_CUBIC) {
    for (j = 0; j < p->npartitions; j ++) {
      for (i = 0; i < p->nlocalparameters; i ++) {
	RJMCMC_INTCHECKINT(write_double(fp, p->model_gradients[j].local_parameter[i]),
			   "part2d_forwardmodel_save: failed to write local parameter gradient\n");
      }
    }
  }
    
  fclose(fp);

  return 0;
}

part2d_forwardmodel_t *
part2d_forwardmodel_load(const char *filename)
{
  FILE *fp;
  int type;
  int min_partitions;
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;
  int nglobalparameters;
  int nlocalparameters;
  int nhierarchicalparameters;
  int include_corners;

  bbox2d_t bound;

  int i;
  int j;
  
  double x;
  double y;

  part2d_forwardmodel_t *r;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return NULL;
  }

  RJMCMC_INTCHECKPTR(read_int(fp, &type), "part2d_forwardmodel_load: failed to read type.\n");
  
  RJMCMC_INTCHECKPTR(read_int(fp, &min_partitions), 
		     "part2d_forwardmodel_load: failed to read min partitions.\n");
  RJMCMC_INTCHECKPTR(read_int(fp, &max_partitions), 
		     "part2d_forwardmodel_load: failed to read max partitions.\n");
  RJMCMC_INTCHECKPTR(read_int(fp, &include_corners), 
		     "part2d_forwardmodel_load: failed to read include corners.\n");

  RJMCMC_INTCHECKPTR(read_double(fp, &xmin),
		     "part2d_forwardmodel_load: failed to read xmin\n");
  RJMCMC_INTCHECKPTR(read_double(fp, &xmax),
		     "part2d_forwardmodel_load: failed to read xmax\n");
  RJMCMC_INTCHECKPTR(read_double(fp, &ymin),
		     "part2d_forwardmodel_load: failed to read ymin\n");
  RJMCMC_INTCHECKPTR(read_double(fp, &ymax),
		     "part2d_forwardmodel_load: failed to read ymax\n");

  RJMCMC_INTCHECKPTR(read_double(fp, &pdx),
		     "part2d_forwardmodel_load: failed to read pdx\n");
  RJMCMC_INTCHECKPTR(read_double(fp, &pdy),
		     "part2d_forwardmodel_load: failed to read pdy\n");

  RJMCMC_INTCHECKPTR(read_int(fp, &nglobalparameters),
		     "part2d_forwardmodel_load: failed to read nglobalparameters\n");
  RJMCMC_INTCHECKPTR(read_int(fp, &nlocalparameters),
		     "part2d_forwardmodel_load: failed to read nlocalparameters\n");
  RJMCMC_INTCHECKPTR(read_int(fp, &nhierarchicalparameters),
		     "part2d_forwardmodel_load: failed to read nhierarchicalparameters\n");


  r = part2d_forwardmodel_create(type,
				 min_partitions,
				 max_partitions,
				 xmin,
				 xmax,
				 ymin,
				 ymax,
				 pdx,
				 pdy,
				 nglobalparameters,
				 nlocalparameters,
				 nhierarchicalparameters,
				 include_corners);
  if (r == NULL) {
    return NULL;
  }

  /*
   * Global parameters
   */
  for (i = 0; i < r->nglobalparameters; i ++) {
    RJMCMC_INTCHECKPTR(read_double(fp, &(r->global_parameter[i])),
		       "part2d_forwardmodel_load: failed to read global parameter\n");
  }

  /*
   * Hierarchical parameters 
   */
  for (i = 0; i < r->nhierarchicalparameters; i ++) {
    RJMCMC_INTCHECKPTR(read_double(fp, &(r->hierarchical_parameter[i])),
		       "part2d_forwardmodel_load: failed to read hierarchical parameter\n");
  }

  /*
   * Point locations
   */
  RJMCMC_INTCHECKPTR(read_int(fp, &(r->npartitions)),
		     "part2d_forwardmodel_load: failed to read npartitions\n");
  for (i = 0; i < r->npartitions; i ++) {
    
    RJMCMC_INTCHECKPTR(read_double(fp, &x),
		       "part2d_forwardmodel_load: failed to read position x\n");
    RJMCMC_INTCHECKPTR(read_double(fp, &y),
		       "part2d_forwardmodel_load: failed to read position y\n");

    RJMCMC_INTCHECKPTR(position_map2d_insert(r->p, x, y, &bound),
		       "part2d_forwardmodel_load: failed to add point\n");
      
  }
  

  /*
   * Local parameters
   */
  for (j = 0; j < r->npartitions; j ++) {
    for (i = 0; i < r->nlocalparameters; i ++) {
      RJMCMC_INTCHECKPTR(read_double(fp, &(r->models[j].local_parameter[i])),
			 "part2d_forwardmodel_load: failed to read local parameter\n");
    }
  }
  
  /*
   * Local gradients (Monotone cubic is not implemented yet).
   */
  if (r->type == PART2D_FM_MONOTONE_CUBIC) {
    for (j = 0; j < r->npartitions; j ++) {
      for (i = 0; i < r->nlocalparameters; i ++) {
	RJMCMC_INTCHECKPTR(read_double(fp, &(r->model_gradients[j].local_parameter[i])),
			   "part2d_forwardmodel_load: failed to read local parameter gradient\n");
      }
    }
  }
    
  fclose(fp);

  return r;
}


void
part2d_forwardmodel_clone(const part2d_forwardmodel_t *src,
			  part2d_forwardmodel_t *dst)
{
  int i;

  RJMCMC_NULLCHECKVOID(src, "part2d_forwardmodel_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "part2d_forwardmodel_clone: null dst\n");

  RJMCMC_CONDITIONCHECKVOID(src->type != dst->type,
			    "part2d_forwardmodel_clone: type mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "part2d_forwardmodel_clone: size mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->nglobalparameters != dst->nglobalparameters,
			    "part2d_forwardmodel_clone: global mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->nlocalparameters != dst->nlocalparameters,
			    "part2d_forwardmodel_clone: local mismatch\n");
  RJMCMC_CONDITIONCHECKVOID(src->nhierarchicalparameters != 
			    dst->nhierarchicalparameters,
			    "part2d_forwardmodel_clone: "
			    "hierarchicall mismatch\n");

  position_map2d_clone(src->p, dst->p);
  models_clone(src->max_partitions,
	       src->nlocalparameters,
	       src->models, dst->models);
  dst->npartitions = src->npartitions;

  for (i = 0; i < src->nglobalparameters; i ++) {
    dst->global_parameter[i] = src->global_parameter[i];
  }

  for (i = 0; i < src->nhierarchicalparameters; i ++) {
    dst->hierarchical_parameter[i] = src->hierarchical_parameter[i];
  }

}

int 
part2d_forwardmodel_initialize(part2d_forwardmodel_t *p,
			       const forwardmodelparameter_t *global_parameters,
			       int nglobalparameters,
			       const forwardmodelparameter_t *local_parameters,
			       int nlocalparameters,
			       const forwardmodelparameter_t *hierarchical_parameters,
			       int nhierarchicalparameters,
			       int initial_partitions,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal)
{
  int npart;
  int pi;
  int gi;
  int hi;
  int li;

  double x;
  double y;

  int i;

  bbox2d_t bound;

  RJMCMC_CONDITIONCHECKINT(nglobalparameters != p->nglobalparameters,
			   "part2d_forwardmodel_initialize: "
			   "global mismatch\n");
  RJMCMC_CONDITIONCHECKINT(nlocalparameters != p->nlocalparameters,
			   "part2d_forwardmodel_initialize: "
			   "local mismatch\n");
  RJMCMC_CONDITIONCHECKINT(nhierarchicalparameters != 
			   p->nhierarchicalparameters,
			   "part2d_forwardmodel_initialize: "
			   "hierarchical mismatch\n");

  npart = initial_partitions;
  if (npart < p->min_partitions) {
    npart = p->min_partitions;
  }
  if (npart > p->max_partitions) {
    npart = p->max_partitions;
  }
  
  for (gi = 0; gi < nglobalparameters; gi++) {
    p->global_parameter[gi] =
      rjmcmc_random_choose_double(global_parameters[gi].fmin,
  				  global_parameters[gi].fmax,
  				  random);
  }

  for (hi = 0; hi < nhierarchicalparameters; hi++) {
    p->hierarchical_parameter[hi] =
      rjmcmc_random_choose_double(hierarchical_parameters[hi].fmin,
  				  hierarchical_parameters[hi].fmax,
  				  random);
  }

  for (pi = 0; pi < npart; pi ++) {
    x = rjmcmc_random_choose_double(p->xmin,
    				    p->xmax,
    				    random);
    y = rjmcmc_random_choose_double(p->ymin,
    				    p->ymax,
    				    random);

    position_map2d_insert(p->p, 
			  x,
			  y,
			  &bound);
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

  return 0;
}

int
part2d_forwardmodel_value_at(const part2d_forwardmodel_t *current,
			     double x,
			     double y,
			     int nlocalparameters,
			     double *localparameters)
{
  int iy;
  int i;
  double nx;
  double ny;

  RJMCMC_NULLCHECKINT(current, "part2d_forwardmodel_value_at: null object\n");
  RJMCMC_CONDITIONCHECKINT(nlocalparameters != current->nlocalparameters,
			   "part2d_forwardmodel_value_at: local mismatch\n");

  if (x < current->xmin || x > current->xmax) {
    printf("x range: %f %f %f\n", x, current->xmin, current->xmax);
    return -1;
  }

  if (y < current->ymin || y > current->ymax) {
    printf("y range: %f %f %f\n", y, current->ymin, current->ymax);
    return -1;
  }

  iy = position_map2d_nearest(current->p, x, y, &nx, &ny, 0);
  if (iy < 0) {
    printf("failed to find point\n");
    return -1;
  }
  iy -= 4;

  for (i = 0; i < nlocalparameters; i ++) {
    localparameters[i] = current->models[iy].local_parameter[i];
  }

  return 0;
}

int 
part2d_forwardmodel_gradient_at(const part2d_forwardmodel_t *current,
				double x,
				double y,
				int nlocalparameters,
				double *localparameters)
{
  return -1;
}

int
part2d_forwardmodel_propose_birth(const part2d_forwardmodel_t *current,
				  part2d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameters,
				  rjmcmc_uniform_rand_t uniform,
				  rjmcmc_normal_rand_t normal,
				  bbox2d_t *bound,
				  double *birth_prob)
{
  double new_x;
  double new_y;

  double prob;

  int old_iy;
  double old_x;
  double old_y;

  int li;

  double dv;

  part2d_forwardmodel_clone(current, proposed);

  if ((proposed->npartitions + 1) >= proposed->max_partitions) {
    /* rjmcmc_error( */
    /* 	    "part2d_forwardmodel_propose_birth: " */
    /* 	    "%d %d\n",  */
    /* 	    current->npartitions, */
    /* 	    current->max_partitions); */
    return 0;
  }
  

  new_x = rjmcmc_random_choose_double(proposed->xmin, proposed->xmax, uniform);
  new_y = rjmcmc_random_choose_double(proposed->ymin, proposed->ymax, uniform);

  /* Find the enclosing partition in the old model for the new point */
  old_iy = position_map2d_nearest(current->p, new_x, new_y, &old_x, &old_y, 0);
  if (old_iy < 0) {
    rjmcmc_error(
	    "part2d_forwardmodel_propose_birth: failed to find old partition\n");
    return 0;
  }
  old_iy -= 4;

  prob = 1.0;


  for (li = 0; li < proposed->nlocalparameters; li ++) {
    
    if (local_parameters[li].fstd_bd > 0.0) {

      /*
       * Birth with Gaussian proposal
       */

      dv = local_parameters[li].fstd_bd * normal();
      proposed->t_lp[li] = current->models[old_iy].local_parameter[li] + dv;
      
      if ((proposed->t_lp[li] < local_parameters[li].fmin) ||
	  (proposed->t_lp[li] > local_parameters[li].fmax)) {
	return 0;
      }
      
      prob *= rjmcmc_gaussian_probability(dv, local_parameters[li].fstd_bd);

    } else {
      
      /*
       * Birth from prior
       */

      proposed->t_lp[li] = uniform() * (local_parameters[li].fmax - local_parameters[li].fmin) + 
	local_parameters[li].fmin;

      prob *= 1.0/(local_parameters[li].fmax - local_parameters[li].fmin);

    }
  }

  if (part2d_forwardmodel_addpoint(proposed, new_x, new_y,
				   proposed->nlocalparameters,
				   proposed->t_lp,
				   bound) < 0) {
    return 0;
  }

  *birth_prob = prob;
  return -1;
}

int 
part2d_forwardmodel_propose_death(const part2d_forwardmodel_t *current,
				  part2d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameters,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  bbox2d_t *bound,
				  double *death_prob)
{
  int del_iy;
  double del_x;
  double del_y;

  int new_iy;
  double new_x;
  double new_y;

  double prob;

  int li;
  double dv;

  part2d_forwardmodel_clone(current, proposed);

  if (proposed->npartitions <= proposed->min_partitions) {
    return 0;
  }
  
  del_iy = rjmcmc_random_choose_int(0, proposed->npartitions - 1, random);
  if (position_map2d_position_of_index(proposed->p, del_iy + 4, &del_x, &del_y) < 0) {
    return 0;
  }

  if (part2d_forwardmodel_delpoint(proposed, del_iy, bound) < 0) {
    return 0;
  }
  
  new_iy = position_map2d_nearest(proposed->p, del_x, del_y, &new_x, &new_y, 0);
  if (new_iy < 0) {
    rjmcmc_error("part2d_forwardmodel_propose_death: failed to find predecessor %d\n", proposed->npartitions);
    return 0;
  }

  new_iy -= 4;

  prob = 1.0;

  for (li = 0; li < proposed->nlocalparameters; li ++) {
    
    dv = 
      proposed->models[new_iy].local_parameter[li] - 
      current->models[del_iy].local_parameter[li];

    if (local_parameters[li].fstd_bd > 0.0) {
      prob *= rjmcmc_gaussian_probability(dv, local_parameters[li].fstd_bd);
    } else {
      prob *= 1.0/(local_parameters[li].fmax - local_parameters[li].fmin);
    }
		
  }

  /* if (position_map2d_validate(proposed->p) < 0) { */
  /*   fprintf(stderr, "part2d_forwardmodel_propose_death: proposed invalidate\n"); */
  /*   exit(-1); */
  /* } */

  *death_prob = prob;
  return -1;
}

int 
part2d_forwardmodel_propose_local_value(const part2d_forwardmodel_t *current, 
					part2d_forwardmodel_t *proposed,
					int nglobalparameters,
					const forwardmodelparameter_t *global_parameters,
					int nlocalparameters,
					const forwardmodelparameter_t *local_parameters,
					rjmcmc_uniform_rand_t random,
					rjmcmc_normal_rand_t normal,
					bbox2d_t *bound,
					double *value_prob)
{
  int li;
  int iy;

  part2d_forwardmodel_clone(current, proposed);

  if (nlocalparameters == 1) {
    li = 0;
  } else {
    li = rjmcmc_random_choose_int(0,
				  nlocalparameters - 1,
				  random);
  }

  iy = rjmcmc_random_choose_int(0, 
				proposed->npartitions - 1,
				random);

  proposed->models[iy].local_parameter[li] += local_parameters[li].fstd_value * normal();

  if (proposed->models[iy].local_parameter[li] < local_parameters[li].fmin ||
      proposed->models[iy].local_parameter[li] > local_parameters[li].fmax) {
    return 0;
  }

  *value_prob = 1.0;
  return 1;
}


int 
part2d_forwardmodel_propose_global_value(const part2d_forwardmodel_t *current, 
					 part2d_forwardmodel_t *proposed,
					 int nglobalparameters,
					 const forwardmodelparameter_t *global_parameters,
					 int nlocalparameters,
					 const forwardmodelparameter_t *local_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *value_prob)
{
  int gi;

  part2d_forwardmodel_clone(current, proposed);

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

  *value_prob = 1.0;
  return 1;
}

int 
part2d_forwardmodel_propose_move(const part2d_forwardmodel_t *current,
				 part2d_forwardmodel_t *proposed,
				 int nglobalparameters,
				 const forwardmodelparameter_t *global_parameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *local_parameters,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 bbox2d_t *bound,
				 double *move_prob)
{
  int move_iy;

  double old_x;
  double old_y;

  double new_x;
  double new_y;
  
  part2d_forwardmodel_clone(current, proposed);

  move_iy = rjmcmc_random_choose_int(0, 
				     proposed->npartitions - 1,
				     random);

  if (position_map2d_position_of_index(proposed->p, move_iy + 4, &old_x, &old_y) < 0) {
    rjmcmc_error(
	    "part2d_forwardmodel_propose_move: "
	    "failed to find position of selected index\n");
    return 0;
  }

  new_x = old_x + normal() * proposed->pdx;
  if (new_x <= proposed->xmin ||
      new_x >= proposed->xmax) {
    return 0;
  }

  new_y = old_y + normal() * proposed->pdy;
  if (new_y <= proposed->ymin ||
      new_y >= proposed->ymax) {
    return 0;
  }

  if (part2d_forwardmodel_movepoint(proposed,
				    move_iy,
				    new_x,
				    new_y,
				    bound) < 0){
    return 0;
  }

  *move_prob = 1.0;

  return 1;
}

int
part2d_forwardmodel_propose_hierarchical(const part2d_forwardmodel_t *current,
					 part2d_forwardmodel_t *proposed,
					 int nhierarchicalparameters,
					 const forwardmodelparameter_t *hierarchical_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *hierarchical_prob)
{
  double dv;
  int hierarchicalindex;

  part2d_forwardmodel_clone(current, proposed);

  if (nhierarchicalparameters > 1) {
    hierarchicalindex = rjmcmc_random_choose_int(0, 
						 nhierarchicalparameters - 1, 
						 random);
  } else {
    hierarchicalindex = 0;
  }

  dv = normal() * hierarchical_parameters[hierarchicalindex].fstd_value;
  proposed->hierarchical_parameter[hierarchicalindex] += dv;

  /* printf("ph: %f %f %f %f %f\n", */
  /* 	 hierarchical_parameters[hierarchicalindex].fmin, */
  /* 	 hierarchical_parameters[hierarchicalindex].fmax, */
  /* 	 hierarchical_parameters[hierarchicalindex].fstd_value, */
  /* 	 dv, */
  /* 	 proposed->hierarchical_parameter[hierarchicalindex]); */
  if ((proposed->hierarchical_parameter[hierarchicalindex] < 
       hierarchical_parameters[hierarchicalindex].fmin) ||
      (proposed->hierarchical_parameter[hierarchicalindex] > 
       hierarchical_parameters[hierarchicalindex].fmax)) {

    return 0;
  }

  *hierarchical_prob = 1.0;
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
part2d_forwardmodel_partitions(const part2d_forwardmodel_t *current)
{
  return current->npartitions;
}

const double *
part2d_forwardmodel_global_parameters(const part2d_forwardmodel_t *current)
{
  return current->global_parameter;
}

const double *
part2d_forwardmodel_hierarchical_parameters(const part2d_forwardmodel_t *current)
{
  return current->hierarchical_parameter;
}

int
part2d_forwardmodel_evaluate_local_parameters(const part2d_forwardmodel_t *c,
					      int xsamples,
					      const double *x,
					      int ysamples,
					      const double *y,
					      int nlocalparameters,
					      double ***local_parameters)

{
  int ix;
  int iy;

  int mi;
  int li;

  double nx;
  double ny;
  
  switch (c->type) {

  case PART2D_FM_ZERO:
    for (iy = 0; iy < ysamples; iy ++) {
      for (ix = 0; ix < xsamples; ix ++) {

	mi = position_map2d_nearest(c->p, x[ix], y[iy], &nx, &ny, 0);
	if (mi < 0) {
	  return -1;
	}

	mi -= 4;

	for (li = 0; li < nlocalparameters; li ++) {
	  local_parameters[li][ix][iy] = c->models[mi].local_parameter[li];
	}
      }
    }
    break;

  default:
    rjmcmc_error("part2d_forwardmodel_evaluate_local_parameters: "
		 "invalid type\n");
    return -1;
  }

  return 0;
}

int
part2d_forwardmodel_partition_centre(const part2d_forwardmodel_t *c,
				     int pi,
				     double *x,
				     double *y)
{
  return position_map2d_position_of_index(c->p, pi, x, y);
}

part2d_fm_type_t
part2d_forwardmodel_type(const part2d_forwardmodel_t *p)
{
  return p->type;
}

int 
part2d_forwardmodel_min_partitions(const part2d_forwardmodel_t *p)
{
  return p->min_partitions;
}

int 
part2d_forwardmodel_max_partitions(const part2d_forwardmodel_t *p)
{
  return p->max_partitions;
}

double
part2d_forwardmodel_pdx(const part2d_forwardmodel_t *p)
{
  return p->pdx;
}

double
part2d_forwardmodel_pdy(const part2d_forwardmodel_t *p)
{
  return p->pdy;
}

static int write_int(FILE *fp, int i)
{
  if (fwrite(&i, sizeof(int), 1, fp) != 1) {
    return -1;
  }
  return 0;
}

static int write_double(FILE *fp, double d)
{
  if (fwrite(&d, sizeof(double), 1, fp) != 1) {
    return -1;
  }
  return 0;
}

static int read_int(FILE *fp, int *i)
{
  if (fread(i, sizeof(int), 1, fp) != 1) {
    return -1;
  }
  return 0;
}

static int read_double(FILE *fp, double *d)
{
  if (fread(d, sizeof(double), 1, fp) != 1) {
    return -1;
  }
  return 0;
}

int
part2d_forwardmodel_addpoint(part2d_forwardmodel_t *c,
			     double x,
			     double y,
			     int nlocalparameters,
			     const double *parameters,
			     bbox2d_t *bound)
{
  int new_iy;
  int i;

  new_iy = c->npartitions;

  if (position_map2d_insert(c->p, x, y, bound) < 0) {
    rjmcmc_error(
	    "part2d_forwardmodel_propose_birth: "
	    "failed to add new point\n");
    return -1;
  }

  c->npartitions ++;

  for (i = 0; i < nlocalparameters; i ++) {
    c->models[new_iy].local_parameter[i] = parameters[i];
  }

  return 0;
}

int
part2d_forwardmodel_delpoint(part2d_forwardmodel_t *c,
			     int pi,
			     bbox2d_t *bound)
{
  if (position_map2d_delete(c->p, pi + 4, bound) < 0) {
    rjmcmc_error("part2d_forwardmodel_delpoint: failed to delete position\n");
    return -1;
  }

  models_delete(c->max_partitions,
		c->nlocalparameters,
		pi,
		c->npartitions,
		c->models);
  c->npartitions --;

  return 0;
}

int
part2d_forwardmodel_movepoint(part2d_forwardmodel_t *c,
			      int pi,
			      double new_x,
			      double new_y,
			      bbox2d_t *bound)
{
  if (position_map2d_move(c->p,
			  pi + 4,
			  new_x,
			  new_y,
			  bound) < 0) {
    rjmcmc_error(
	    "part2d_forwardmodel_movepoint: "
	    "failed to move point\n");

    return -1;
  }

  return 0;
}

