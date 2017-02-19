#ifndef part2d_forwardmodel_h
#define part2d_forwardmodel_h

#include <rjmcmc/forwardmodelparameter.h>
#include <rjmcmc/bbox2d.h>
#include <rjmcmc/rjmcmc_random.h>

typedef struct _part2d_forwardmodel part2d_forwardmodel_t;

typedef enum {
  PART2D_FM_ZERO,
  PART2D_FM_NATURAL,
  PART2D_FM_MONOTONE_CUBIC
} part2d_fm_type_t;

part2d_forwardmodel_t *
part2d_forwardmodel_create(part2d_fm_type_t type,
			   int min_part,
			   int max_part,
			   double xmin, 
			   double xmax,
			   double ymin,
			   double ymax,
			   double pdx,
			   double pdy,
			   int nglobalparameters,
			   int nlocalparameters,
			   int nhierarchicalparameters,
			   int includecorners);

void
part2d_forwardmodel_destroy(part2d_forwardmodel_t *p);

int
part2d_forwardmodel_save(const part2d_forwardmodel_t *p,
			 const char *filename);

part2d_forwardmodel_t *
part2d_forwardmodel_load(const char *filename);

void
part2d_forwardmodel_clone(const part2d_forwardmodel_t *src,
			  part2d_forwardmodel_t *dst);

int
part2d_forwardmodel_value_at(const part2d_forwardmodel_t *current,
			     double x,
			     double y,
			     int nlocalparameters,
			     double *localparameters);

int 
part2d_forwardmodel_gradient_at(const part2d_forwardmodel_t *current,
				double x,
				double y,
				int nlocalparameters,
				double *localparameters);

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
			       rjmcmc_normal_rand_t normal);

int
part2d_forwardmodel_propose_birth(const part2d_forwardmodel_t *current,
				  part2d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameter,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  bbox2d_t *bound,
				  double *birth_prob);

int 
part2d_forwardmodel_propose_death(const part2d_forwardmodel_t *current,
				  part2d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameter,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  bbox2d_t *bound,
				  double *death_prob);

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
					double *value_prob);

int 
part2d_forwardmodel_propose_local_value_natural(const part2d_forwardmodel_t *current, 
						part2d_forwardmodel_t *proposed,
						int nglobalparameters,
						const forwardmodelparameter_t *global_parameters,
						int nlocalparameters,
						const forwardmodelparameter_t *local_parameters,
						rjmcmc_uniform_rand_t random,
						rjmcmc_normal_rand_t normal,
						bbox2d_t *bound,
						double *value_prob);

int 
part2d_forwardmodel_propose_global_value(const part2d_forwardmodel_t *current, 
					 part2d_forwardmodel_t *proposed,
					 int nglobalparameters,
					 const forwardmodelparameter_t *global_parameters,
					 int nlocalparameters,
					 const forwardmodelparameter_t *local_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *value_prob);

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
				 double *move_prob);

int
part2d_forwardmodel_propose_hierarchical(const part2d_forwardmodel_t *current,
					 part2d_forwardmodel_t *proposed,
					 int nhierarchicalparameters,
					 const forwardmodelparameter_t *hierarchical_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *hierarchical_prob);

int
part2d_forwardmodel_partitions(const part2d_forwardmodel_t *current);

int
part2d_forwardmodel_partition_centre(const part2d_forwardmodel_t *c,
				     int pi,
				     double *x,
				     double *y);

const double *
part2d_forwardmodel_global_parameters(const part2d_forwardmodel_t *current);

const double *
part2d_forwardmodel_hierarchical_parameters(const part2d_forwardmodel_t *current);

int
part2d_forwardmodel_evaluate_local_parameters(const part2d_forwardmodel_t *c,
					      int xsamples,
					      const double *x,
					      int ysamples,
					      const double *y,
					      int nlocalparameters,
					      double ***local_parameters);

int
part2d_forwardmodel_partition_centre(const part2d_forwardmodel_t *c,
				     int pi,
				     double *x,
				     double *y);

part2d_fm_type_t
part2d_forwardmodel_type(const part2d_forwardmodel_t *p);

int 
part2d_forwardmodel_min_partitions(const part2d_forwardmodel_t *p);

int 
part2d_forwardmodel_max_partitions(const part2d_forwardmodel_t *p);

double
part2d_forwardmodel_pdx(const part2d_forwardmodel_t *p);

double
part2d_forwardmodel_pdy(const part2d_forwardmodel_t *p);
	  

/*
 * Internal functions (exposed only for testing)
 */
int
part2d_forwardmodel_addpoint(part2d_forwardmodel_t *c,
			     double x,
			     double y,
			     int nlocalparameters,
			     const double *parameters,
			     bbox2d_t *bound);

int
part2d_forwardmodel_delpoint(part2d_forwardmodel_t *c,
			     int pi,
			     bbox2d_t *bound);

int
part2d_forwardmodel_movepoint(part2d_forwardmodel_t *c,
			      int pi,
			      double new_x,
			      double new_y,
			      bbox2d_t *bound);

#endif /* part2d_forwardmodel_h */
