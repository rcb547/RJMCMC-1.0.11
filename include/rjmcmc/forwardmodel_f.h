#ifndef forwardmodel_f_h
#define forwardmodel_f_h

#include <rjmcmc/forwardmodel.h>

resultsetfm_t *
single_forwardmodel_f(int burnin,
		      int total,
		      rjmcmc_uniform_rand_t *random,
		      rjmcmc_normal_rand_t *normal,
		      int nparameters,
		      const forwardmodelparameter_t *parameters,
		      single_fm_likelihood_t *lp,
		      void **user_arg,
		      int samples,
		      double credible_interval,
		      int results);

resultsetfm_t *
single_forwardmodel_hierarchical_f(int burnin,
				   int total,
				   rjmcmc_uniform_rand_t *random,
				   rjmcmc_normal_rand_t *normal,
				   int nparameters,
				   const forwardmodelparameter_t *parameters,
				   int nhierarchicalparameters,
				   const forwardmodelparameter_t *hierarchicalparameters,
				   single_fm_likelihood_hierarchical_t *lp,
				   void **user_arg,
				   int samples,
				   double credible_interval,
				   int results);

typedef double (*part1d_fm_likelihood_f_t)(void *userarg,
					   int npartitions,
					   const double *partition_boundaries,
					   int nglobalparameters,
					   const double *global_parameters,
					   part1d_fm_likelihood_state_t *state,
					   part1d_fm_value_at_t value_at,
					   part1d_fm_value_at_t gradient_at);

typedef double (*part1d_fm_hierarchical_likelihood_f_t)(void *userarg,
							int npartitions,
							const double *partition_boundaries,
							int nglobalparameters,
							const double *global_parameters,
							int hierarchical,
							int nhierarchicalparameters,
							const double *hierarchical_parameters,
							part1d_fm_hierarchical_likelihood_state_t *state,
							part1d_fm_value_at_t value_at,
							part1d_fm_value_at_t gradient_at,
							double *logdetce);

resultset1dfm_t *
part1d_forwardmodel_f(int burnin,
		      int total,
		      int minpart,
		      int maxpart,
		      double minx,
		      double maxx,
		      int xsamples,
		      int ysamples,
		      double credible_interval,
		      double pd,
		      rjmcmc_uniform_rand_t *random,
		      rjmcmc_normal_rand_t *normal,
		      int nglobalparameters,
		      const forwardmodelparameter_t *global_parameters,
		      int nlocalparameters,
		      const forwardmodelparameter_t *local_parameters,
		      part1d_fm_likelihood_f_t *lp,
		      void *user_arg,
		      int results);

resultset1dfm_t *
part1d_forwardmodel_hierarchical_f(int burnin,
				   int total,
				   int minpart,
				   int maxpart,
				   double minx,
				   double maxx,
				   int xsamples,
				   int ysamples,
				   double credible_interval,
				   double pd,
				   rjmcmc_uniform_rand_t *random,
				   rjmcmc_normal_rand_t *normal,
				   int nglobalparameters,
				   const forwardmodelparameter_t *
				   global_parameters,
				   int nlocalparameters,
				   const forwardmodelparameter_t *
				   local_parameters,
				   int nhierarchicalparameters,
				   const forwardmodelparameter_t *
				   hierarchical_parameters,
				   part1d_fm_hierarchical_likelihood_f_t *lp,
				   void *user_arg,
				   int results);

resultset1dfm_t *
part1d_forwardmodel_natural_f(int burnin,
			      int total,
			      int minpart,
			      int maxpart,
			      double minx,
			      double maxx,
			      int xsamples,
			      int ysamples,
			      double credible_interval,
			      double pd,
			      rjmcmc_uniform_rand_t *random,
			      rjmcmc_normal_rand_t *normal,
			      int nglobalparameters,
			      const forwardmodelparameter_t *global_parameters,
			      int nlocalparameters,
			      const forwardmodelparameter_t *local_parameters,
			      part1d_fm_likelihood_f_t *lp,
			      void *user_arg,
			      int results);

resultset1dfm_t *
part1d_forwardmodel_natural_hierarchical_f(int burnin,
					   int total,
					   int minpart,
					   int maxpart,
					   double minx,
					   double maxx,
					   int xsamples,
					   int ysamples,
					   double credible_interval,
					   double pd,
					   rjmcmc_uniform_rand_t *random,
					   rjmcmc_normal_rand_t *normal,
					   int nglobalparameters,
					   const forwardmodelparameter_t *
					   global_parameters,
					   int nlocalparameters,
					   const forwardmodelparameter_t *
					   local_parameters,
					   int nhierarchicalparameters,
					   const forwardmodelparameter_t *
					   hierarchical_parameters,
					   part1d_fm_hierarchical_likelihood_f_t *lp,
					   void *user_arg,
					   int results);

void *
part1d_forwardmodel_natural_hierarchical_init_f(int burnin,
						int total,
						int minpart,
						int maxpart,
						double minx,
						double maxx,
						int xsamples,
						int ysamples,
						double credible_interval,
						double pd,
						rjmcmc_uniform_rand_t *random,
						rjmcmc_normal_rand_t *normal,
						int nglobalparameters,
						const forwardmodelparameter_t *
						global_parameters,
						int nlocalparameters,
						const forwardmodelparameter_t *
						local_parameters,
						int nhierarchicalparameters,
						const forwardmodelparameter_t *
						hierarchical_parameters,
						part1d_fm_hierarchical_likelihood_f_t *lp,
						void *user_arg,
						int results);

int
part1d_forwardmodel_natural_hierarchical_step_f(void *state);

resultset1dfm_t *
part1d_forwardmodel_natural_hierarchical_finish_f(void *state);


/*
 * 2D Forward model wrappers 
 */
typedef double (*part2d_fm_likelihood_f_t)(void *userarg,
					   int nglobalparameters,
					   const double *global_parameters,
					   part2d_fm_likelihood_state_t *state,
					   part2d_fm_value_at_t value_at,
					   part2d_fm_value_at_t gradient_at,
					   const bbox2d_t *bound);

typedef double (*part2d_fm_hierarchical_likelihood_f_t)(
  void *userarg,
  int nglobalparameters,
  const double *global_parameters,
  int hierarchical,
  int nhierarchicalvalues,
  const double *hierarchical_values,
  part2d_fm_likelihood_state_t *state,
  part2d_fm_value_at_t value_at,
  part2d_fm_value_at_t gradient_at,
  const bbox2d_t *bound,
  double *logdetce);

resultset2dfm_t *
part2d_forwardmodel_f(int burnin,
		      int total,
		      int thin,
		      int minpart,
		      int maxpart,
		      int initpart,
		      double minx,
		      double maxx,
		      double miny,
		      double maxy,
		      int xsamples,
		      int ysamples,
		      int zsamples,
		      double credible_interval,
		      double pdx,
		      double pdy,
		      rjmcmc_uniform_rand_t *random,
		      rjmcmc_normal_rand_t *normal,
		      int nglobalparameters,
		      const forwardmodelparameter_t *global_parameters,
		      int nlocalparameters,
		      const forwardmodelparameter_t *local_parameters,
		      part2d_fm_likelihood_f_t *lp,
		      void *user_arg,
		      int results);

resultset2dfm_t *
part2d_forwardmodel_hierarchical_f(int burnin,
				   int total,
				   int thin, 
				   int minpart,
				   int maxpart,
				   int initpart,
				   double minx,
				   double maxx,
				   double miny,
				   double maxy,
				   int xsamples,
				   int ysamples,
				   int zsamples,
				   double credible_interval,
				   double pdx,
				   double pdy,
				   rjmcmc_uniform_rand_t *random,
				   rjmcmc_normal_rand_t *normal,
				   int nglobalparameters,
				   const forwardmodelparameter_t *
				   global_parameters,
				   int nlocalparameters,
				   const forwardmodelparameter_t *
				   local_parameters,
				   int nhierarchicalparameters,
				   const forwardmodelparameter_t *
				   hierarchical_parameters,
				   part2d_fm_hierarchical_likelihood_f_t *lp,
				   void *user_arg,
				   int results);

#if defined(HAVE_MPI_H)

resultsetfm_t *
mpi_single_forwardmodel_f(int burnin,
			  int total, 
			  rjmcmc_uniform_rand_t *random,
			  rjmcmc_normal_rand_t *normal,
			  int nparameters,
			  const forwardmodelparameter_t *parameters,
			  single_fm_likelihood_t *lp,
			  void *user_arg,
			  int samples,
			  double credible_interval,
			  int results,
			  int mpisize,
			  int mpirank,
			  int root,
			  int comm);

resultsetfm_t *
mpi_single_forwardmodel_hierarchical_f(int burnin,
				       int total, 
				       rjmcmc_uniform_rand_t *random,
				       rjmcmc_normal_rand_t *normal,
				       int nparameters,
				       const forwardmodelparameter_t *parameters,
				       int nhierarchicalparameters,
				       const forwardmodelparameter_t *hierarchicalparameters,
				       single_fm_likelihood_hierarchical_t *lp,
				       void *user_arg,
				       int samples,
				       double credible_interval,
				       int results,
				       int mpisize,
				       int mpirank,
				       int root,
				       int comm);

resultset1dfm_t *
mpi_part1d_forwardmodel_f(int burnin,
			  int total,
			  int minpart,
			  int maxpart,
			  double minx,
			  double maxx,
			  int xsamples,
			  int ysamples,
			  double credible_interval,
			  double pd,
			  rjmcmc_uniform_rand_t *random,
			  rjmcmc_normal_rand_t *normal,
			  int nglobalparameters,
			  const forwardmodelparameter_t *global_parameters,
			  int nlocalparameters,
			  const forwardmodelparameter_t *local_parameters,
			  part1d_fm_likelihood_f_t *lp,
			  void *user_arg,
			  int results,
			  int mpisize,
			  int mpirank,
			  int root,
			  int comm);

resultset1dfm_t *
mpi_part1d_forwardmodel_natural_f(int burnin,
				  int total,
				  int minpart,
				  int maxpart,
				  double xmin,
				  double xmax,
				  int xsamples,
				  int ysamples,
				  double credible_interval,
				  double pd,
				  rjmcmc_uniform_rand_t *random,
				  rjmcmc_normal_rand_t *normal,
				  int nglobalparameters,
				  const forwardmodelparameter_t *global_parameters,
				  int nlocalparameters,
				  const forwardmodelparameter_t *local_parameters,
				  part1d_fm_likelihood_f_t *lp,
				  void *user_arg,
				  int results,
				  int mpisize,
				  int mpirank,
				  int root,
				  int comm);


resultset1dfm_t *
mpi_part1d_forwardmodel_hierarchical_f(int burnin,
				       int total,
				       int minpart,
				       int maxpart,
				       double minx,
				       double maxx,
				       int xsamples,
				       int ysamples,
				       double credible_interval,
				       double pd,
				       rjmcmc_uniform_rand_t *random,
				       rjmcmc_normal_rand_t *normal,
				       int nglobalparameters,
				       const forwardmodelparameter_t *
				       global_parameters,
				       int nlocalparameters,
				       const forwardmodelparameter_t *
				       local_parameters,
				       int nhierarchicalparameters,
				       const forwardmodelparameter_t *
				       hierarchical_parameters,
				       part1d_fm_hierarchical_likelihood_f_t *lp,
				       void *user_arg,
				       int results,
				       int mpisize,
				       int mpirank,
				       int root,
				       int comm);

resultset1dfm_t *
mpi_part1d_forwardmodel_natural_hierarchical_f(int burnin,
					       int total,
					       int minpart,
					       int maxpart,
					       double xmin,
					       double xmax,
					       int xsamples,
					       int ysamples,
					       double credible_interval,
					       double pd,
					       rjmcmc_uniform_rand_t *random,
					       rjmcmc_normal_rand_t *normal,
					       int nglobalparameters,
					       const forwardmodelparameter_t *
					       global_parameters,
					       int nlocalparameters,
					       const forwardmodelparameter_t *
					       local_parameters,
					       int nhierarchicalparameters,
					       const forwardmodelparameter_t *
					       hierarchical_parameters,
					       part1d_fm_hierarchical_likelihood_f_t *lp,
					       void *user_arg,
					       int results,
					       int mpisize,
					       int mpirank,
					       int root,
					       int comm);


resultset2dfm_t *
mpi_part2d_forwardmodel_f(int burnin,
			  int total,
			  int thin,
			  int minpart,
			  int maxpart,
			  int initpart,
			  double minx,
			  double maxx,
			  double miny,
			  double maxy,
			  int xsamples,
			  int ysamples,
			  int zsamples,
			  double credible_interval,
			  double pdx,
			  double pdy,
			  rjmcmc_uniform_rand_t *random,
			  rjmcmc_normal_rand_t *normal,
			  int nglobalparameters,
			  const forwardmodelparameter_t *global_parameters,
			  int nlocalparameters,
			  const forwardmodelparameter_t *local_parameters,
			  part2d_fm_likelihood_f_t *lp,
			  void *user_arg,
			  int results,
			  int mpisize,
			  int mpirank,
			  int root,
			  int comm);

resultset2dfm_t *
mpi_part2d_forwardmodel_hierarchical_f(int burnin,
				       int total,
				       int thin,
				       int minpart,
				       int maxpart,
				       int initpart,
				       double minx,
				       double maxx,
				       double miny,
				       double maxy,
				       int xsamples,
				       int ysamples,
				       int zsamples,
				       double credible_interval,
				       double pdx,
				       double pdy,
				       rjmcmc_uniform_rand_t *random,
				       rjmcmc_normal_rand_t *normal,
				       int nglobalparameters,
				       const forwardmodelparameter_t *
				       global_parameters,
				       int nlocalparameters,
				       const forwardmodelparameter_t *
				       local_parameters,
				       int nhierarchicalparameters,
				       const forwardmodelparameter_t *
				       hierarchical_parameters,
				       part2d_fm_hierarchical_likelihood_f_t *lp,
				       void *user_arg,
				       int results,
				       int mpisize,
				       int mpirank,
				       int root,
				       int comm);

resultset2dfm_t *
mpi_part2d_forwardmodel_hierarchical_restartable_f(const char *infile_template,
						   const char *outfile_template,
						   int burnin,
						   int total,
						   int thin,
						   int minpart,
						   int maxpart,
						   int initpart,
						   double minx,
						   double maxx,
						   double miny,
						   double maxy,
						   int xsamples,
						   int ysamples,
						   int zsamples,
						   double credible_interval,
						   double pdx,
						   double pdy,
						   rjmcmc_uniform_rand_t *random,
						   rjmcmc_normal_rand_t *normal,
						   int nglobalparameters,
						   const forwardmodelparameter_t *
						   global_parameters,
						   int nlocalparameters,
						   const forwardmodelparameter_t *
						   local_parameters,
						   int nhierarchicalparameters,
						   const forwardmodelparameter_t *
						   hierarchical_parameters,
						   part2d_fm_hierarchical_likelihood_f_t *lp,
						   void *user_arg,
						   int results,
						   int mpisize,
						   int mpirank,
						   int root,
						   int comm);

#endif /* HAVE_MPI_H */

#endif /* forwardmodel_f_h */
