#ifndef part1d_forwardmodel_h
#define part1d_forwardmodel_h

#include <rjmcmc/forwardmodelparameter.h>
#include <rjmcmc/rjmcmc_random.h>

typedef struct _part1d_forwardmodel part1d_forwardmodel_t;

typedef enum {
  PART1D_FM_ZERO,
  PART1D_FM_NATURAL,
  PART1D_FM_ZERO_CUBIC
} part1d_fm_type_t;

part1d_forwardmodel_t *
part1d_forwardmodel_create(part1d_fm_type_t type,
			   int min_part,
			   int max_part,
			   double xmin, 
			   double xmax,
			   double pd,
			   int nglobalparameters,
			   int nlocalparameters,
			   int nhierarchicalparameters);

void
part1d_forwardmodel_destroy(part1d_forwardmodel_t *p);

void
part1d_forwardmodel_clone(const part1d_forwardmodel_t *src,
			  part1d_forwardmodel_t *dst);

int
part1d_forwardmodel_value_at(const part1d_forwardmodel_t *current,
			     double x,
			     int nlocalparameters,
			     double *localparameters);

int 
part1d_forwardmodel_gradient_at(const part1d_forwardmodel_t *current,
				double x,
				int nlocalparameters,
				double *localparameters);

int 
part1d_forwardmodel_initialize_fixed(part1d_forwardmodel_t *p,
				     const double *global_parameters,
				     const double *hierarchical_parameters,
				     int npartitions,
				     const double *partition_x,
				     const double **local_parameters);

int 
part1d_forwardmodel_initialize(part1d_forwardmodel_t *p,
			       const forwardmodelparameter_t *global_parameters,
			       int nglobalparameters,
			       const forwardmodelparameter_t *local_parameters,
			       int nlocalparameters,
			       const forwardmodelparameter_t *hierarchical_parameters,
			       int nhierarchicalparameters,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal);

int
part1d_forwardmodel_propose_birth(const part1d_forwardmodel_t *current,
				  part1d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameter,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *birth_prob);

int 
part1d_forwardmodel_propose_death(const part1d_forwardmodel_t *current,
				  part1d_forwardmodel_t *proposed,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameter,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *death_prob);

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
					int *li);

int 
part1d_forwardmodel_propose_global_value(const part1d_forwardmodel_t *current, 
					 part1d_forwardmodel_t *proposed,
					 int nglobalparameters,
					 const forwardmodelparameter_t *global_parameters,
					 int nlocalparameters,
					 const forwardmodelparameter_t *local_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *value_prob);

int 
part1d_forwardmodel_propose_move(const part1d_forwardmodel_t *current,
				 part1d_forwardmodel_t *proposed,
				 int nglobalparameters,
				 const forwardmodelparameter_t *global_parameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *local_parameters,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 double *move_prob);

int
part1d_forwardmodel_propose_hierarchical(const part1d_forwardmodel_t *current,
					 part1d_forwardmodel_t *proposed,
					 int nhierarchicalparameters,
					 const forwardmodelparameter_t *hierarchical_parameters,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 double *hierarchical_prob);


int
part1d_forwardmodel_partitions(const part1d_forwardmodel_t *current);

double
part1d_forwardmodel_partition_position(const part1d_forwardmodel_t *current,
				       int pi);

int 
part1d_forwardmodel_partition_fill_list(const part1d_forwardmodel_t *current,
					double *positions,
					int *npartitions);

const double *
part1d_forwardmodel_global_parameters(const part1d_forwardmodel_t *current);

const double *
part1d_forwardmodel_hierarchical_parameters(const part1d_forwardmodel_t *current);

int
part1d_forwardmodel_evaluate_local_parameters(const part1d_forwardmodel_t *c,
					      int xsamples,
					      const double *x,
					      double **local_parameters);
				      
int
part1d_forwardmodel_hierarchical_fill_list(const part1d_forwardmodel_t *current,
				    double *hierarchical,
				    int *nhierarchical);

			  

#endif /* part1d_forwardmodel_h */
