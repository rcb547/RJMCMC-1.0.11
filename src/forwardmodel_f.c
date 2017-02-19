
#include <stdlib.h>

#include <rjmcmc/rjmcmc.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/forwardmodel_f.h>

#include <rjmcmc/forwardmodel_mpi.h>

/*
 * These are wrappers to enable fortran programs to use the rjmcmc
 * routines. These deal with most unix environments but there will
 * probably be some cases where pointer arithmetic between C and 
 * fortran compilers won't match.
 *
 * Also note that all functions here most be in all lower case to 
 * enable linking with fortran (gfortran).
 */

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
		      int results)
{
  return single_forwardmodel(burnin,
			     total,
			     *random,
			     *normal,
			     nparameters,
			     parameters,
			     *lp,
			     user_arg,
			     samples,
			     credible_interval,
			     results);
}

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
				   int results)
{
  return single_forwardmodel_hierarchical(burnin,
					  total,
					  *random,
					  *normal,
					  nparameters,
					  parameters,
					  nhierarchicalparameters,
					  hierarchicalparameters,
					  *lp,
					  user_arg,
					  samples,
					  credible_interval,
					  results);
}

struct part1d_fm_value_at_wrap {
  part1d_fm_likelihood_state_t *state;
  part1d_fm_value_at_t value_at;
  part1d_fm_value_at_t gradient_at;
};

static int
part1d_fm_value_at_wrap(void *state, 
			double x,
			int maxsize,
			double *values)
{

  struct part1d_fm_value_at_wrap *wrap = 
    (struct part1d_fm_value_at_wrap*)state;
  int i;
  const double *_values;

  _values = wrap->value_at(wrap->state, x);
  if (_values == NULL) {
    for (i = 0; i < maxsize; i ++) {
      values[i] = 0.0;
    }
  } else {
    for (i = 0; i < maxsize; i ++) {
      values[i] = _values[i];
    }
  }

  return maxsize;
}

static int
part1d_fm_gradient_at_wrap(void *state, 
			   double x,
			   int maxsize,
			   double *gradients)
{
  struct part1d_fm_value_at_wrap *wrap = 
    (struct part1d_fm_value_at_wrap*)state;
  int i;
  const double *_gradients;

  _gradients = wrap->gradient_at(wrap->state, x);

  for (i = 0; i < maxsize; i ++) {
    gradients[i] = _gradients[i];
  }

  return maxsize;
}


struct part1d_forwardmodel_wrap {
  void *user_arg;
  part1d_fm_likelihood_f_t lp;
};

static double 
part1d_fm_likelihood_f_wrap(void *user_arg,
			    int npartitions,
			    const double *partition_boundaries,
			    int nglobalparameters,
			    const double *global_parameters,
			    part1d_fm_likelihood_state_t *state,
			    part1d_fm_value_at_t value_at,
			    part1d_fm_value_at_t gradient_at)
{
  struct part1d_forwardmodel_wrap *wrap;

  struct part1d_fm_value_at_wrap value_wrap;

  value_wrap.state = state;
  value_wrap.value_at = value_at;
  value_wrap.gradient_at = gradient_at;

  wrap = (struct part1d_forwardmodel_wrap*)user_arg;

  return (wrap->lp)(wrap->user_arg,
		    npartitions,
		    partition_boundaries,
		    nglobalparameters,
		    global_parameters,
		    (part1d_fm_likelihood_state_t*)&value_wrap,
		    (part1d_fm_value_at_t)&part1d_fm_value_at_wrap,
		    (part1d_fm_value_at_t)&part1d_fm_gradient_at_wrap);
}

struct part1d_hierarchical_forwardmodel_wrap {
  void *user_arg;
  part1d_fm_hierarchical_likelihood_f_t lp;
  void *step_state;
};

static double 
part1d_fm_hierarchical_likelihood_f_wrap(void *user_arg,
					 int npartitions,
					 const double *partition_boundaries,
					 int nglobalparameters,
					 const double *global_parameters,
					 int hierarchical,
					 int nhierarchicalparameters,
					 const double *hierarchical_parameters,
					 part1d_fm_hierarchical_likelihood_state_t *state,
					 part1d_fm_hierarchical_value_at_t value_at,
					 part1d_fm_hierarchical_value_at_t gradient_at,
					 double *logdetce)
{
  struct part1d_hierarchical_forwardmodel_wrap *wrap;

  struct part1d_fm_value_at_wrap value_wrap;

  value_wrap.state = (part1d_fm_likelihood_state_t *)state;
  value_wrap.value_at = (part1d_fm_value_at_t)value_at;
  value_wrap.gradient_at = (part1d_fm_value_at_t)gradient_at;

  wrap = (struct part1d_hierarchical_forwardmodel_wrap*)user_arg;

  return (wrap->lp)(wrap->user_arg,
		    npartitions,
		    partition_boundaries,
		    nglobalparameters,
		    global_parameters,
		    hierarchical,
		    nhierarchicalparameters,
		    hierarchical_parameters,
		    (part1d_fm_hierarchical_likelihood_state_t*)&value_wrap,
		    (part1d_fm_value_at_t)&part1d_fm_value_at_wrap,
		    (part1d_fm_value_at_t)&part1d_fm_gradient_at_wrap,
		    logdetce);
}

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
		      int results)
{
  struct part1d_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return part1d_forwardmodel(burnin,
			     total,
			     minpart,
			     maxpart,
			     minx,
			     maxx,
			     xsamples,
			     ysamples,
			     credible_interval,
			     pd,
			     *random,
			     *normal,
			     nglobalparameters,
			     global_parameters,
			     nlocalparameters,
			     local_parameters,
			     part1d_fm_likelihood_f_wrap,
			     &wrap,
			     results);
}

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
				   int results)
{
  struct part1d_hierarchical_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return part1d_forwardmodel_hierarchical(burnin,
					  total,
					  minpart,
					  maxpart,
					  minx,
					  maxx,
					  xsamples,
					  ysamples,
					  credible_interval,
					  pd,
					  *random,
					  *normal,
					  nglobalparameters,
					  global_parameters,
					  nlocalparameters,
					  local_parameters,
					  nhierarchicalparameters,
					  hierarchical_parameters,
					  part1d_fm_hierarchical_likelihood_f_wrap,
					  &wrap,
					  results);
}

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
			      int results)
{
  struct part1d_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return part1d_forwardmodel_natural(burnin,
				     total,
				     minpart,
				     maxpart,
				     minx,
				     maxx,
				     xsamples,
				     ysamples,
				     credible_interval,
				     pd,
				     *random,
				     *normal,
				     nglobalparameters,
				     global_parameters,
				     nlocalparameters,
				     local_parameters,
				     part1d_fm_likelihood_f_wrap,
				     &wrap,
				     results);
}

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
					   int results)
{
  struct part1d_hierarchical_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return part1d_forwardmodel_natural_hierarchical(burnin,
						  total,
						  minpart,
						  maxpart,
						  minx,
						  maxx,
						  xsamples,
						  ysamples,
						  credible_interval,
						  pd,
						  *random,
						  *normal,
						  nglobalparameters,
						  global_parameters,
						  nlocalparameters,
						  local_parameters,
						  nhierarchicalparameters,
						  hierarchical_parameters,
						  part1d_fm_hierarchical_likelihood_f_wrap,
						  &wrap,
						  results);
}

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
						int results)
{
  struct part1d_hierarchical_forwardmodel_wrap *wrap;

  wrap = (struct part1d_hierarchical_forwardmodel_wrap *)
    malloc(sizeof(struct part1d_hierarchical_forwardmodel_wrap));

  wrap->user_arg = user_arg;
  wrap->lp = *lp;

  wrap->step_state = 
    part1d_forwardmodel_natural_hierarchical_init(burnin,
						  total,
						  minpart,
						  maxpart,
						  minx,
						  maxx,
						  xsamples,
						  ysamples,
						  credible_interval,
						  pd,
						  *random,
						  *normal,
						  nglobalparameters,
						  global_parameters,
						  nlocalparameters,
						  local_parameters,
						  nhierarchicalparameters,
						  hierarchical_parameters,
						  part1d_fm_hierarchical_likelihood_f_wrap,
						  wrap,
						  results);
  if (wrap->step_state == NULL) {
    free(wrap);
    return NULL;
  }

  return wrap;
}

int
part1d_forwardmodel_natural_hierarchical_step_f(void *state)
{
  struct part1d_hierarchical_forwardmodel_wrap *s = 
    (struct part1d_hierarchical_forwardmodel_wrap *)state;
  return part1d_forwardmodel_natural_hierarchical_step(s->step_state);

}

resultset1dfm_t *
part1d_forwardmodel_natural_hierarchical_finish_f(void *state)
{
  resultset1dfm_t *r;
  struct part1d_hierarchical_forwardmodel_wrap *s = 
    (struct part1d_hierarchical_forwardmodel_wrap *)state;

  r = part1d_forwardmodel_natural_hierarchical_finish(s->step_state);
  free(s);

  return r;
}


/*
 * 2D Forward model wrappers
 */

struct part2d_fm_value_at_wrap {
  part2d_fm_likelihood_state_t *state;
  part2d_fm_value_at_t value_at;
  part2d_fm_value_at_t gradient_at;
};

struct part2d_forwardmodel_wrap {
  void *user_arg;
  part2d_fm_likelihood_f_t lp;
};

struct part2d_forwardmodel_hierarchical_wrap {
  void *user_arg;
  part2d_fm_hierarchical_likelihood_f_t lp;
};

static int
part2d_fm_value_at_wrap(void *state, 
			double x,
			double y,
			int maxsize,
			double *values)
{
  struct part2d_fm_value_at_wrap *wrap = 
    (struct part2d_fm_value_at_wrap*)state;
  int i;
  const double *_values;

  _values = wrap->value_at(wrap->state, x, y);
  if (_values == NULL) {
    for (i = 0; i < maxsize; i ++) {
      values[i] = 0.0;
    }
  } else {
    for (i = 0; i < maxsize; i ++) {
      values[i] = _values[i];
    }
  }

  return maxsize;
}

static int
part2d_fm_gradient_at_wrap(void *state, 
			   double x,
			   double y,
			   int maxsize,
			   double *gradients)
{
  struct part2d_fm_value_at_wrap *wrap = 
    (struct part2d_fm_value_at_wrap*)state;
  int i;
  const double *_gradients;

  _gradients = wrap->gradient_at(wrap->state, x, y);

  for (i = 0; i < maxsize; i ++) {
    gradients[i] = _gradients[i];
  }

  return maxsize;
}

static double 
part2d_fm_likelihood_f_wrap(void *user_arg,
			    int nglobalparameters,
			    const double *global_parameters,
			    part2d_fm_likelihood_state_t *state,
			    part2d_fm_value_at_t value_at,
			    part2d_fm_value_at_t gradient_at,
			    const bbox2d_t *bound)
{
  struct part2d_forwardmodel_wrap *wrap;

  struct part2d_fm_value_at_wrap value_wrap;

  value_wrap.state = state;
  value_wrap.value_at = value_at;
  value_wrap.gradient_at = gradient_at;

  wrap = (struct part2d_forwardmodel_wrap*)user_arg;

  return (wrap->lp)(wrap->user_arg,
		    nglobalparameters,
		    global_parameters,
		    (part2d_fm_likelihood_state_t*)&value_wrap,
		    (part2d_fm_value_at_t)&part2d_fm_value_at_wrap,
		    (part2d_fm_value_at_t)&part2d_fm_gradient_at_wrap,
		    bound);
}

static double 
part2d_fm_hierarchical_likelihood_f_wrap(void *user_arg,
					 int nglobalparameters,
					 const double *global_parameters,
					 int hierarchical,
					 int nhierarchicalvalues,
					 const double *hierarchical_values,
					 part2d_fm_likelihood_state_t *state,
					 part2d_fm_value_at_t value_at,
					 part2d_fm_value_at_t gradient_at,
					 const bbox2d_t *bound,
					 double *logdetce)
{
  struct part2d_forwardmodel_hierarchical_wrap *wrap;

  struct part2d_fm_value_at_wrap value_wrap;

  value_wrap.state = state;
  value_wrap.value_at = value_at;
  value_wrap.gradient_at = gradient_at;

  wrap = (struct part2d_forwardmodel_hierarchical_wrap*)user_arg;

  return (wrap->lp)(wrap->user_arg,
		    nglobalparameters,
		    global_parameters,
		    hierarchical,
		    nhierarchicalvalues,
		    hierarchical_values,
		    (part2d_fm_likelihood_state_t*)&value_wrap,
		    (part2d_fm_value_at_t)&part2d_fm_value_at_wrap,
		    (part2d_fm_value_at_t)&part2d_fm_gradient_at_wrap,
		    bound,
		    logdetce);
}

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
		      int results)
{
  struct part2d_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return part2d_forwardmodel(burnin,
			     total,
			     thin, 
			     minpart,
			     maxpart,
			     initpart,
			     minx,
			     maxx,
			     miny,
			     maxy,
			     xsamples,
			     ysamples,
			     zsamples,
			     credible_interval,
			     pdx,
			     pdy,
			     *random,
			     *normal,
			     nglobalparameters,
			     global_parameters,
			     nlocalparameters,
			     local_parameters,
			     part2d_fm_likelihood_f_wrap,
			     &wrap,
			     results);
}

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
				   int results)
{
  struct part2d_forwardmodel_hierarchical_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return part2d_forwardmodel_hierarchical(burnin,
					  total,
					  thin,
					  minpart,
					  maxpart,
					  initpart,
					  minx,
					  maxx,
					  miny,
					  maxy,
					  xsamples,
					  ysamples,
					  zsamples,
					  credible_interval,
					  pdx,
					  pdy,
					  *random,
					  *normal,
					  nglobalparameters,
					  global_parameters,
					  nlocalparameters,
					  local_parameters,
					  nhierarchicalparameters,
					  hierarchical_parameters,
					  part2d_fm_hierarchical_likelihood_f_wrap,
					  &wrap,
					  results);
}

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
			  int comm)
{
  return MPI_single_forwardmodel(burnin,
				 total,
				 *random,
				 *normal,
				 nparameters,
				 parameters,
				 *lp,
				 user_arg,
				 samples,
				 credible_interval,
				 results,
				 mpisize,
				 mpirank,
				 root,
				 MPI_COMM_WORLD);
}

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
				       int comm)
{
  return MPI_single_forwardmodel_hierarchical(burnin,
					      total,
					      *random,
					      *normal,
					      nparameters,
					      parameters,
					      nhierarchicalparameters,
					      hierarchicalparameters,
					      *lp,
					      user_arg,
					      samples,
					      credible_interval,
					      results,
					      mpisize,
					      mpirank,
					      root,
					      MPI_COMM_WORLD);
}

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
			  int comm)
{
  struct part1d_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return MPI_part1d_forwardmodel(burnin,
				 total,
				 minpart,
				 maxpart,
				 minx,
				 maxx,
				 xsamples,
				 ysamples,
				 credible_interval,
				 pd,
				 *random,
				 *normal,
				 nglobalparameters,
				 global_parameters,
				 nlocalparameters,
				 local_parameters,
				 part1d_fm_likelihood_f_wrap,
				 &wrap,
				 results,
				 mpisize,
				 mpirank,
				 root,
				 MPI_COMM_WORLD);
}

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
				  int comm)
{
  struct part1d_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return MPI_part1d_forwardmodel_natural(burnin,
					 total,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 xsamples,
					 ysamples,
					 credible_interval,
					 pd,
					 *random,
					 *normal,
					 nglobalparameters,
					 global_parameters,
					 nlocalparameters,
					 local_parameters,
					 part1d_fm_likelihood_f_wrap,
					 &wrap,
					 results,
					 mpisize,
					 mpirank,
					 root,
					 MPI_COMM_WORLD);
}


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
				       int comm)
{
  struct part1d_hierarchical_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return MPI_part1d_forwardmodel_hierarchical(burnin,
					      total,
					      minpart,
					      maxpart,
					      minx,
					      maxx,
					      xsamples,
					      ysamples,
					      credible_interval,
					      pd,
					      *random,
					      *normal,
					      nglobalparameters,
					      global_parameters,
					      nlocalparameters,
					      local_parameters,
					      nhierarchicalparameters,
					      hierarchical_parameters,
					      part1d_fm_hierarchical_likelihood_f_wrap,
					      &wrap,
					      results,
					      mpisize,
					      mpirank,
					      root,
					      MPI_COMM_WORLD);
}

resultset1dfm_t *
mpi_part1d_forwardmodel_natural_hierarchical_f(int burnin,
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
					       int comm)
{
  struct part1d_hierarchical_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return MPI_part1d_forwardmodel_natural_hierarchical(burnin,
						      total,
						      minpart,
						      maxpart,
						      minx,
						      maxx,
						      xsamples,
						      ysamples,
						      credible_interval,
						      pd,
						      *random,
						      *normal,
						      nglobalparameters,
						      global_parameters,
						      nlocalparameters,
						      local_parameters,
						      nhierarchicalparameters,
						      hierarchical_parameters,
						      part1d_fm_hierarchical_likelihood_f_wrap,
						      &wrap,
						      results,
						      mpisize,
						      mpirank,
						      root,
						      MPI_COMM_WORLD);
}


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
			  int comm)
{
  struct part2d_forwardmodel_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return MPI_part2d_forwardmodel(burnin,
				 total,
				 thin,
				 minpart,
				 maxpart,
				 initpart,
				 minx,
				 maxx,
				 miny,
				 maxy,
				 xsamples,
				 ysamples,
				 zsamples,
				 credible_interval,
				 pdx,
				 pdy,
				 *random,
				 *normal,
				 nglobalparameters,
				 global_parameters,
				 nlocalparameters,
				 local_parameters,
				 part2d_fm_likelihood_f_wrap,
				 &wrap,
				 results,
				 mpisize,
				 mpirank,
				 root,
				 MPI_COMM_WORLD);
}

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
				       int comm)
{
  struct part2d_forwardmodel_hierarchical_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return MPI_part2d_forwardmodel_hierarchical(burnin,
					      total,
					      thin,
					      minpart,
					      maxpart,
					      initpart,
					      minx,
					      maxx,
					      miny,
					      maxy,
					      xsamples,
					      ysamples,
					      zsamples,
					      credible_interval,
					      pdx,
					      pdy,
					      *random,
					      *normal,
					      nglobalparameters,
					      global_parameters,
					      nlocalparameters,
					      local_parameters,
					      nhierarchicalparameters,
					      hierarchical_parameters,
					      part2d_fm_hierarchical_likelihood_f_wrap,
					      &wrap,
					      results,
					      mpisize,
					      mpirank,
					      root,
					      MPI_COMM_WORLD);
}

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
						   int comm)
{
  struct part2d_forwardmodel_hierarchical_wrap wrap;

  wrap.user_arg = user_arg;
  wrap.lp = *lp;

  return MPI_part2d_forwardmodel_hierarchical_restartable(infile_template,
							  outfile_template,
							  burnin,
							  total,
							  thin,
							  minpart,
							  maxpart,
							  initpart,
							  minx,
							  maxx,
							  miny,
							  maxy,
							  xsamples,
							  ysamples,
							  zsamples,
							  credible_interval,
							  pdx,
							  pdy,
							  *random,
							  *normal,
							  nglobalparameters,
							  global_parameters,
							  nlocalparameters,
							  local_parameters,
							  nhierarchicalparameters,
							  hierarchical_parameters,
							  part2d_fm_hierarchical_likelihood_f_wrap,
							  &wrap,
							  results,
							  mpisize,
							  mpirank,
							  root,
							  MPI_COMM_WORLD);
}

#endif /* HAVE_MPI_H */
