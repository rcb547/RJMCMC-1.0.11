
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/engine.h>
#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/forwardmodel_mpi.h>
#include <rjmcmc/part1d_forwardmodel.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_debug.h>
 
static double part1d_init(void *arg);
static int part1d_select(void *arg);
static void *part1d_perturb(void *arg, int proc);
static double part1d_misfit(void *arg, void *state);
static int part1d_accept(void *arg, 
			 double current,
			 double proposed);
static int part1d_sample(void *arg,
			 int i);


struct part1dfm {

  resultset1dfm_t *results;

  part1d_forwardmodel_t *current;
  double current_like;

  part1d_forwardmodel_t *proposed;
  double proposed_like;

  int minpart;
  int maxpart;

  double xmin;
  double xmax;
  
  int nprocesses;

  int out;
  int accepted;
  int process;

  double birth_prob;
  double death_prob;
  double move_prob;
  double value_prob;
  double global_value_prob;

  int nglobalparameters;
  const forwardmodelparameter_t *global_parameters;

  int nlocalparameters;
  const forwardmodelparameter_t *local_parameters;

  double proddelta_l;

  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;

  part1d_fm_likelihood_t lp;
  void *user_arg;

  int xsamples;

  double *mf_global_parameters;

  double *mf_values;
  double *mf_gradients;

  double *mf_partitions;

  double *x;
  double **y;

  int index;
  rjmcmc_engine_cb_t cb;

};
  
struct _part1d_fm_likelihood_state {
  part1d_forwardmodel_t *model;

  int nlocalparameters;
  double *values;
  double *gradients;
};

static const double *
part1d_fm_likelihood_value_callback(part1d_fm_likelihood_state_t *state,
				    double x)
{
  if (part1d_forwardmodel_value_at(state->model,
				   x,
				   state->nlocalparameters,
				   state->values) < 0) {
    return NULL;
  } else {
    return state->values;
  }
}

static const double *
part1d_fm_likelihood_gradient_callback(part1d_fm_likelihood_state_t *state,
				       double x)
{
  if (part1d_forwardmodel_gradient_at(state->model,
				      x,
				      state->nlocalparameters,
				      state->gradients) < 0) {
    return NULL;
  } else {
    return state->gradients;
  }
}

resultset1dfm_t *
part1d_forwardmodel(int burnin,
		    int total,
		    int minpart,
		    int maxpart,
		    double xmin,
		    double xmax,
		    int xsamples,
		    int ysamples,
		    double credible_interval,
		    double pd,
		    rjmcmc_uniform_rand_t random,
		    rjmcmc_normal_rand_t normal,
		    int nglobalparameters,
		    const forwardmodelparameter_t *global_parameters,
		    int nlocalparameters,
		    const forwardmodelparameter_t *local_parameters,
		    part1d_fm_likelihood_t lp,
		    void *user_arg,
		    int results)
{
  struct part1dfm s;
  int i;
  int xs;

  if (nlocalparameters <= 0) {
    rjmcmc_error("part1d_forwardmodel: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part1d_forwardmodel: "
		 "a forward model function must be provided\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = 4;
  } else {
    s.nprocesses = 5;
  }

  if (minpart < 2) {
    rjmcmc_error(
	    "part1d_forwardmodel: "
	    "minimum partitions must be 2 or greater\n");
    return NULL;
  }

  if (maxpart <= minpart) {
    rjmcmc_error(
	    "part1d_forwardmodel: "
	    "maximum number of partitions must be greater than the minimum\n");
    return NULL;
  }
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.xmin = xmin;
  s.xmax = xmax;
  
  s.results = resultset1dfm_create(burnin,
				   total,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   0,
				   xsamples,
				   ysamples,
				   maxpart,
				   xmin,
				   xmax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = part1d_forwardmodel_create(PART1D_FM_ZERO,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 pd,
					 nglobalparameters,
					 nlocalparameters,
					 0);
  if (s.current == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s.proposed = part1d_forwardmodel_create(PART1D_FM_ZERO,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  pd,
					  nglobalparameters,
					  nlocalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_global_parameters = NULL;
  if (nglobalparameters > 0) {
    s.mf_global_parameters = rjmcmc_create_array_1d(nglobalparameters);
    if (s.mf_global_parameters == NULL) {
      return NULL;
    }
  }

  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
    return NULL;
  }

  s.mf_partitions = rjmcmc_create_array_1d(maxpart);
  if (s.mf_partitions == NULL) {
    return NULL;
  }
  
  s.mf_gradients = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_gradients == NULL) {
    return NULL;
  }

  s.x = rjmcmc_create_array_1d(xsamples);
  if (s.x == NULL) {
    return NULL;
  }

  s.y = rjmcmc_create_array_2d(nlocalparameters, xsamples);
  if (s.y == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset1dfm_fill_xcoord_vector(s.results, s.x, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    if (local_parameters[i].fmax <= local_parameters[i].fmin) {
      rjmcmc_error("part1d_forwardmodel: invalid local parameter bounds\n");
      return NULL;
    }
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  s.cb.initialize_state = part1d_init;
  s.cb.select_process = part1d_select;
  s.cb.perturb_process = part1d_perturb;
  s.cb.compute_misfit = part1d_misfit;
  s.cb.accept = part1d_accept;
  s.cb.sample = part1d_sample;
  s.cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&s.cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_global_parameters);
  rjmcmc_destroy_array_1d(s.mf_partitions);
  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);
  
  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_2d(s.nlocalparameters, s.y);

  part1d_forwardmodel_destroy(s.current);
  part1d_forwardmodel_destroy(s.proposed);

  resultset1dfm_assemble_results(s.results);

  return s.results;
}

resultset1dfm_t *
part1d_forwardmodel_natural(int burnin,
			    int total,
			    int minpart,
			    int maxpart,
			    double xmin,
			    double xmax,
			    int xsamples,
			    int ysamples,
			    double credible_interval,
			    double pd,
			    rjmcmc_uniform_rand_t random,
			    rjmcmc_normal_rand_t normal,
			    int nglobalparameters,
			    const forwardmodelparameter_t *global_parameters,
			    int nlocalparameters,
			    const forwardmodelparameter_t *local_parameters,
			    part1d_fm_likelihood_t lp,
			    void *user_arg,
			    int results)
{
  struct part1dfm s;
  int i;
  int xs;

  memset(&s, 0, sizeof(struct part1dfm));

  if (nlocalparameters <= 0) {
    rjmcmc_error(
	    "part1d_forwardmodel_natural: "
	    "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = 4;
  } else {
    s.nprocesses = 5;
  }
  
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.xmin = xmin;
  s.xmax = xmax;
  
  s.results = resultset1dfm_create(burnin,
				   total,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   0,
				   xsamples,
				   ysamples,
				   maxpart,
				   xmin,
				   xmax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = part1d_forwardmodel_create(PART1D_FM_NATURAL,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 pd,
					 nglobalparameters,
					 nlocalparameters,
					 0);
  if (s.current == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s.proposed = part1d_forwardmodel_create(PART1D_FM_NATURAL,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  pd,
					  nglobalparameters,
					  nlocalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_global_parameters = NULL;
  if (nglobalparameters > 0) {
    s.mf_global_parameters = rjmcmc_create_array_1d(nglobalparameters);
    if (s.mf_global_parameters == NULL) {
      return NULL;
    }
  }

  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
    return NULL;
  }

  s.mf_partitions = rjmcmc_create_array_1d(maxpart);
  if (s.mf_partitions == NULL) {
    return NULL;
  }
  
  s.mf_gradients = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_gradients == NULL) {
    return NULL;
  }

  s.x = rjmcmc_create_array_1d(xsamples);
  if (s.x == NULL) {
    return NULL;
  }

  s.y = rjmcmc_create_array_2d(nlocalparameters, xsamples);
  if (s.y == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset1dfm_fill_xcoord_vector(s.results, s.x, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  s.cb.initialize_state = part1d_init;
  s.cb.select_process = part1d_select;
  s.cb.perturb_process = part1d_perturb;
  s.cb.compute_misfit = part1d_misfit;
  s.cb.accept = part1d_accept;
  s.cb.sample = part1d_sample;
  s.cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&s.cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_global_parameters);
  rjmcmc_destroy_array_1d(s.mf_partitions);
  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);

  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_2d(s.nlocalparameters, s.y);

  part1d_forwardmodel_destroy(s.current);
  part1d_forwardmodel_destroy(s.proposed);

  resultset1dfm_assemble_results(s.results);

  return s.results;
}

void *
part1d_forwardmodel_natural_init(int burnin,
				 int total,
				 int minpart,
				 int maxpart,
				 double xmin,
				 double xmax,
				 int xsamples,
				 int ysamples,
				 double credible_interval,
				 double pd,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 int nglobalparameters,
				 const forwardmodelparameter_t *global_parameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *local_parameters,
				 part1d_fm_likelihood_t lp,
				 void *user_arg,
				 int results)
{
  struct part1dfm *s;
  int i;
  int xs;

  s = (struct part1dfm*)malloc(sizeof(struct part1dfm));
  if (s == NULL) {
    return NULL;
  }
  memset(s, 0, sizeof(struct part1dfm));

  if (nlocalparameters <= 0) {
    rjmcmc_error(
	    "part1d_forwardmodel_natural: "
	    "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s->nprocesses = 4;
  } else {
    s->nprocesses = 5;
  }
  
  s->minpart = minpart;
  s->maxpart = maxpart;
  s->xmin = xmin;
  s->xmax = xmax;
  
  s->results = resultset1dfm_create(burnin,
				   total,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   0,
				   xsamples,
				   ysamples,
				   maxpart,
				   xmin,
				   xmax,
				   s->nprocesses,
				   credible_interval,
				   results);
  if (s->results == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s->current = part1d_forwardmodel_create(PART1D_FM_NATURAL,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 pd,
					 nglobalparameters,
					 nlocalparameters,
					 0);
  if (s->current == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s->proposed = part1d_forwardmodel_create(PART1D_FM_NATURAL,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  pd,
					  nglobalparameters,
					  nlocalparameters,
					  0);
  if (s->proposed == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s->nglobalparameters = nglobalparameters;
  s->global_parameters = global_parameters;
  
  s->nlocalparameters = nlocalparameters;
  s->local_parameters = local_parameters;
  
  s->random = random;
  s->normal = normal;

  s->xsamples = xsamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s->mf_global_parameters = NULL;
  if (nglobalparameters > 0) {
    s->mf_global_parameters = rjmcmc_create_array_1d(nglobalparameters);
    if (s->mf_global_parameters == NULL) {
      return NULL;
    }
  }

  s->mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s->mf_values == NULL) {
    return NULL;
  }

  s->mf_partitions = rjmcmc_create_array_1d(maxpart);
  if (s->mf_partitions == NULL) {
    return NULL;
  }
  
  s->mf_gradients = rjmcmc_create_array_1d(nlocalparameters);
  if (s->mf_gradients == NULL) {
    return NULL;
  }

  s->x = rjmcmc_create_array_1d(xsamples);
  if (s->x == NULL) {
    return NULL;
  }

  s->y = rjmcmc_create_array_2d(nlocalparameters, xsamples);
  if (s->y == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset1dfm_fill_xcoord_vector(s->results, s->x, &xs);

  s->proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s->proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s->lp = lp;
  s->user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  s->cb.initialize_state = part1d_init;
  s->cb.select_process = part1d_select;
  s->cb.perturb_process = part1d_perturb;
  s->cb.compute_misfit = part1d_misfit;
  s->cb.accept = part1d_accept;
  s->cb.sample = part1d_sample;
  s->cb.arg = (void*)s;

  if (rjmcmc_engine_init(&(s->cb),
			 burnin,
			 total, 
			 1) < 0) {
    return NULL;
  }

  return s;
}

int
part1d_forwardmodel_natural_step(void *_s) 
{
  struct part1dfm *s = (struct part1dfm *)_s;
  return rjmcmc_engine_step(&s->cb);
}

resultset1dfm_t *
part1d_forwardmodel_natural_finish(void *_s)
{
  struct part1dfm *s = (struct part1dfm *)_s;
  resultset1dfm_t *r;

  rjmcmc_destroy_array_1d(s->mf_global_parameters);
  rjmcmc_destroy_array_1d(s->mf_partitions);
  rjmcmc_destroy_array_1d(s->mf_values);
  rjmcmc_destroy_array_1d(s->mf_gradients);

  rjmcmc_destroy_array_1d(s->x);
  rjmcmc_destroy_array_2d(s->nlocalparameters, s->y);

  part1d_forwardmodel_destroy(s->current);
  part1d_forwardmodel_destroy(s->proposed);

  resultset1dfm_assemble_results(s->results);

  r = s->results;
  free(s);

  return r;
}

resultset1dfm_t *
part1d_forwardmodel_zero_cubic(int burnin,
			       int total,
			       int minpart,
			       int maxpart,
			       double xmin,
			       double xmax,
			       int xsamples,
			       int ysamples,
			       double credible_interval,
			       double pd,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal,
			       int nglobalparameters,
			       const forwardmodelparameter_t *global_parameters,
			       int nlocalparameters,
			       const forwardmodelparameter_t *local_parameters,
			       part1d_fm_likelihood_t lp,
			       void *user_arg,
			       int results)
{
  struct part1dfm s;
  int i;
  int xs;

  memset(&s, 0, sizeof(struct part1dfm));

  if (nlocalparameters <= 0) {
    rjmcmc_error(
	    "part1d_forwardmodel_zero_cubic: "
	    "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = 4;
  } else {
    s.nprocesses = 5;
  }
  
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.xmin = xmin;
  s.xmax = xmax;
  
  s.results = resultset1dfm_create(burnin,
				   total,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   0,
				   xsamples,
				   ysamples,
				   maxpart,
				   xmin,
				   xmax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = part1d_forwardmodel_create(PART1D_FM_ZERO_CUBIC,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 pd,
					 nglobalparameters,
					 nlocalparameters,
					 0);
  if (s.current == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s.proposed = part1d_forwardmodel_create(PART1D_FM_ZERO_CUBIC,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  pd,
					  nglobalparameters,
					  nlocalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_global_parameters = NULL;
  if (nglobalparameters > 0) {
    s.mf_global_parameters = rjmcmc_create_array_1d(nglobalparameters);
    if (s.mf_global_parameters == NULL) {
      return NULL;
    }
  }

  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
    return NULL;
  }

  s.mf_partitions = rjmcmc_create_array_1d(maxpart);
  if (s.mf_partitions == NULL) {
    return NULL;
  }
  
  s.mf_gradients = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_gradients == NULL) {
    return NULL;
  }

  s.x = rjmcmc_create_array_1d(xsamples);
  if (s.x == NULL) {
    return NULL;
  }

  s.y = rjmcmc_create_array_2d(nlocalparameters, xsamples);
  if (s.y == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset1dfm_fill_xcoord_vector(s.results, s.x, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  s.cb.initialize_state = part1d_init;
  s.cb.select_process = part1d_select;
  s.cb.perturb_process = part1d_perturb;
  s.cb.compute_misfit = part1d_misfit;
  s.cb.accept = part1d_accept;
  s.cb.sample = part1d_sample;
  s.cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&s.cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_global_parameters);
  rjmcmc_destroy_array_1d(s.mf_partitions);
  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);

  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_2d(s.nlocalparameters, s.y);

  part1d_forwardmodel_destroy(s.current);
  part1d_forwardmodel_destroy(s.proposed);

  resultset1dfm_assemble_results(s.results);

  return s.results;
}

#if defined(HAVE_MPI_H)

resultset1dfm_t *
MPI_part1d_forwardmodel(int burnin,
			int total,
			int minpart,
			int maxpart,
			double xmin,
			double xmax,
			int xsamples,
			int ysamples,
			double credible_interval,
			double pd,
			rjmcmc_uniform_rand_t random,
			rjmcmc_normal_rand_t normal,
			int nglobalparameters,
			const forwardmodelparameter_t *global_parameters,
			int nlocalparameters,
			const forwardmodelparameter_t *local_parameters,
			part1d_fm_likelihood_t lp,
			void *user_arg,
			int results,
			int mpisize,
			int mpirank,
			int root,
			MPI_Comm comm)
{
  struct part1dfm s;
  int i;
  int xs;

  if (nlocalparameters <= 0) {
    rjmcmc_error(
	    "part1d_forwardmodel: "
	    "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = 4;
  } else {
    s.nprocesses = 5;
  }

  if (minpart < 2) {
    rjmcmc_error(
	    "part1d_forwardmodel: "
	    "minimum partitions must be 2 or greater\n");
    return NULL;
  }

  if (maxpart <= minpart) {
    rjmcmc_error(
	    "part1d_forwardmodel: "
	    "maximum number of partitions must be greater than the minimum\n");
    return NULL;
  }
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.xmin = xmin;
  s.xmax = xmax;
  
  s.results = resultset1dfm_create(burnin,
				   total,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   0,
				   xsamples,
				   ysamples,
				   maxpart,
				   xmin,
				   xmax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = part1d_forwardmodel_create(PART1D_FM_ZERO,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 pd,
					 nglobalparameters,
					 nlocalparameters,
					 0);
  if (s.current == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s.proposed = part1d_forwardmodel_create(PART1D_FM_ZERO,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  pd,
					  nglobalparameters,
					  nlocalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_global_parameters = NULL;
  if (nglobalparameters > 0) {
    s.mf_global_parameters = rjmcmc_create_array_1d(nglobalparameters);
    if (s.mf_global_parameters == NULL) {
      return NULL;
    }
  }

  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
    return NULL;
  }

  s.mf_partitions = rjmcmc_create_array_1d(maxpart);
  if (s.mf_partitions == NULL) {
    return NULL;
  }
  
  s.mf_gradients = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_gradients == NULL) {
    return NULL;
  }

  s.x = rjmcmc_create_array_1d(xsamples);
  if (s.x == NULL) {
    return NULL;
  }

  s.y = rjmcmc_create_array_2d(nlocalparameters, xsamples);
  if (s.y == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset1dfm_fill_xcoord_vector(s.results, s.x, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  s.cb.initialize_state = part1d_init;
  s.cb.select_process = part1d_select;
  s.cb.perturb_process = part1d_perturb;
  s.cb.compute_misfit = part1d_misfit;
  s.cb.accept = part1d_accept;
  s.cb.sample = part1d_sample;
  s.cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&s.cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_global_parameters);
  rjmcmc_destroy_array_1d(s.mf_partitions);
  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);
  
  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_2d(s.nlocalparameters, s.y);

  part1d_forwardmodel_destroy(s.current);
  part1d_forwardmodel_destroy(s.proposed);

  MPI_resultset1dfm_assemble_results(s.results, 
				     mpisize,
				     mpirank,
				     root,
				     comm);

  return s.results;
}

resultset1dfm_t *
MPI_part1d_forwardmodel_natural(int burnin,
				int total,
				int minpart,
				int maxpart,
				double xmin,
				double xmax,
				int xsamples,
				int ysamples,
				double credible_interval,
				double pd,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal,
				int nglobalparameters,
				const forwardmodelparameter_t *
				global_parameters,
				int nlocalparameters,
				const forwardmodelparameter_t *
				local_parameters,
				part1d_fm_likelihood_t lp,
				void *user_arg,
				int results,
				int mpisize,
				int mpirank,
				int root,
				MPI_Comm comm)
{
  struct part1dfm s;
  int i;
  int xs;

  if (nlocalparameters <= 0) {
    rjmcmc_error(
	    "part1d_forwardmodel_natural: "
	    "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = 4;
  } else {
    s.nprocesses = 5;
  }
  
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.xmin = xmin;
  s.xmax = xmax;
  
  s.results = resultset1dfm_create(burnin,
				   total,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   0,
				   xsamples,
				   ysamples,
				   maxpart,
				   xmin,
				   xmax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = part1d_forwardmodel_create(PART1D_FM_NATURAL,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 pd,
					 nglobalparameters,
					 nlocalparameters,
					 0);
  if (s.current == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s.proposed = part1d_forwardmodel_create(PART1D_FM_NATURAL,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  pd,
					  nglobalparameters,
					  nlocalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part1d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_global_parameters = NULL;
  if (nglobalparameters > 0) {
    s.mf_global_parameters = rjmcmc_create_array_1d(nglobalparameters);
    if (s.mf_global_parameters == NULL) {
      return NULL;
    }
  }

  
  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
    return NULL;
  }

  s.mf_partitions = rjmcmc_create_array_1d(maxpart);
  if (s.mf_partitions == NULL) {
    return NULL;
  }
  
  s.mf_gradients = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_gradients == NULL) {
    return NULL;
  }

  s.x = rjmcmc_create_array_1d(xsamples);
  if (s.x == NULL) {
    return NULL;
  }

  s.y = rjmcmc_create_array_2d(nlocalparameters, xsamples);
  if (s.y == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset1dfm_fill_xcoord_vector(s.results, s.x, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  s.cb.initialize_state = part1d_init;
  s.cb.select_process = part1d_select;
  s.cb.perturb_process = part1d_perturb;
  s.cb.compute_misfit = part1d_misfit;
  s.cb.accept = part1d_accept;
  s.cb.sample = part1d_sample;
  s.cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&s.cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_global_parameters);
  rjmcmc_destroy_array_1d(s.mf_partitions);
  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);

  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_2d(s.nlocalparameters, s.y);

  part1d_forwardmodel_destroy(s.current);
  part1d_forwardmodel_destroy(s.proposed);

  MPI_resultset1dfm_assemble_results(s.results, 
				     mpisize,
				     mpirank,
				     root,
				     comm);

  return s.results;
}

#endif /* HAVE_MPI_H */

static double part1d_init(void *arg)
{
  struct part1dfm *s = (struct part1dfm *)arg;
  part1d_fm_likelihood_state_t state;

  int npartitions;

  if (part1d_forwardmodel_initialize(s->current,
				     s->global_parameters,
				     s->nglobalparameters,
				     s->local_parameters,
				     s->nlocalparameters,
				     NULL,
				     0,
				     s->random,
				     s->normal) < 0) {
    rjmcmc_error("part1d_init: failed to initialize\n");
    return -1.0;
  }

  state.model = s->current;
  state.nlocalparameters = s->nlocalparameters;
  state.values = s->mf_values;
  state.gradients = s->mf_gradients;

  npartitions = s->maxpart;
  part1d_forwardmodel_partition_fill_list(s->current,
					  s->mf_partitions,
					  &npartitions);

  s->current_like = s->lp(s->user_arg,
			  npartitions,
			  s->mf_partitions,
			  s->nglobalparameters,
			  s->mf_global_parameters,
			  &state,
			  part1d_fm_likelihood_value_callback,
			  part1d_fm_likelihood_gradient_callback);

  return s->current_like;
}

static int part1d_select(void *arg)
{
  struct part1dfm *s = (struct part1dfm *)arg;

  return rjmcmc_random_choose_int(0, s->nprocesses - 1, s->random);
}

static void *part1d_perturb(void *arg, int proc)
{
  struct part1dfm *s = (struct part1dfm *)arg;

  s->process = proc;

  resultset1dfm_propose(s->results, proc);

  s->birth_prob = 0.0;
  s->death_prob = 0.0;
  s->move_prob = 0.0;
  s->value_prob = 0.0;
  switch(proc) {
  case 0: /* Birth */
    s->out = part1d_forwardmodel_propose_birth(s->current,
					       s->proposed,
					       s->nglobalparameters,
					       s->global_parameters,
					       s->nlocalparameters,
					       s->local_parameters,
					       s->random,
					       s->normal,
					       &(s->birth_prob));
    break;

  case 1: /* Death */
    s->out = part1d_forwardmodel_propose_death(s->current,
					       s->proposed,
					       s->nglobalparameters,
					       s->global_parameters,
					       s->nlocalparameters,
					       s->local_parameters,
					       s->random,
					       s->normal,
					       &(s->death_prob));
    break;

  case 2: /* Move */
    s->out = part1d_forwardmodel_propose_move(s->current,
					      s->proposed,
					      s->nglobalparameters,
					      s->global_parameters,
					      s->nlocalparameters,
					      s->local_parameters,
					      s->random,
					      s->normal,
					      &(s->move_prob));
    break;
    
  case 3: /* Value (local) */
    s->out = part1d_forwardmodel_propose_local_value(s->current,
						     s->proposed,
						     s->nglobalparameters,
						     s->global_parameters,
						     s->nlocalparameters,
						     s->local_parameters,
						     s->random,
						     s->normal,
						     &(s->value_prob),
						     &(s->index));

    resultset1dfm_propose_local_value(s->results, s->index);
    
    break;

  case 4: /* Value (global) */
    s->out = part1d_forwardmodel_propose_global_value(s->current,
						      s->proposed,
						      s->nglobalparameters,
						      s->global_parameters,
						      s->nlocalparameters,
						      s->local_parameters,
						      s->random,
						      s->normal,
						      &(s->value_prob));
    break;

  default:
    
    return NULL;

  }

  if (s->out) {
    return s->proposed;
  } else {
    return NULL;
  }
}

static double part1d_misfit(void *arg, void *_state)
{
  struct part1dfm *s = (struct part1dfm *)arg;
  part1d_fm_likelihood_state_t state;
  int npartitions;

  if (s->out) {
    state.model = s->proposed;
    state.nlocalparameters = s->nlocalparameters;
    state.values = s->mf_values;
    state.gradients = s->mf_gradients;
    
    npartitions = s->maxpart;
    part1d_forwardmodel_partition_fill_list(s->proposed,
					    s->mf_partitions,
					    &npartitions);
    
    s->proposed_like = s->lp(s->user_arg,
			     npartitions,
			     s->mf_partitions,
			     s->nglobalparameters,
			     s->mf_global_parameters,
			     &state,
			     part1d_fm_likelihood_value_callback,
			     part1d_fm_likelihood_gradient_callback);
    
    return s->proposed_like;
  } else {
    return FLT_MAX;
  }
}

static int part1d_accept(void *arg, 
			 double current,
			 double proposed)
{
  struct part1dfm *s = (struct part1dfm *)arg;
  double u;

  s->accepted = 0;

  if (s->out == 0) {
    return 0;
  }

  u = log(s->random());

  switch (s->process) {
  case 0: /* Birth */
    if (u < (current - proposed - 
	     log(s->proddelta_l) - log(s->birth_prob))) {
      s->accepted = 1;
    }
    break;

  case 1: /* Death */
    if (u < (current - proposed +
	     log(s->proddelta_l) + log(s->death_prob))) {
      s->accepted = 1;
    }
    break;

  case 2: /* Move */
    if (u < (current - proposed)) {
      s->accepted = 1;
    }      
    break;

  case 3: /* Value */
    if (u < (current - proposed)) {
      s->accepted = 1;
      resultset1dfm_accept_local_value(s->results, s->index);
    }
    break;

  case 4: /* Global Value */
    if (u < (current - proposed)) {
      s->accepted = 1;
    }
    break;

  default:
    rjmcmc_error("part1d_accept: invalid process %d\n", s->process);
    s->accepted = 0;
    break;
  }

  if (s->accepted) {
    part1d_forwardmodel_clone(s->proposed, s->current);
    s->current_like = s->proposed_like;

    resultset1dfm_accept(s->results, s->process);
  }

  return s->accepted;
}



static int part1d_sample(void *arg,
			 int i)
{
  int gi;
  int pi;
  int li;
  int npartitions;
  struct part1dfm *s = (struct part1dfm *)arg;

  const double *gp;

  gp = part1d_forwardmodel_global_parameters(s->current);
  if (gp == NULL && s->nglobalparameters > 0) {
    return -1;
  }

  for (gi = 0; gi < s->nglobalparameters; gi ++) {
    resultset1dfm_sample_global_parameter(s->results, 
					  i, 
					  gi, 
					  gp[gi]);
  }

  if (part1d_forwardmodel_evaluate_local_parameters(s->current,
						    s->xsamples,
						    s->x,
						    s->y) < 0) {
    return -1;
  }

  for (li = 0; li < s->nlocalparameters; li ++) {

    resultset1dfm_sample_local_parameter(s->results,
					 i,
					 li,
					 s->y[li]);
  }

  npartitions = part1d_forwardmodel_partitions(s->current);
  if (npartitions < 0) {
    return -1;
  }

  resultset1dfm_sample_npartitions(s->results,
				   i,
				   npartitions);

  for (pi = 2; pi < npartitions; pi ++) {
    double x = part1d_forwardmodel_partition_position(s->current,
						      pi);
    resultset1dfm_sample_partition_x(s->results,
				     x);
  }

  resultset1dfm_sample_misfit(s->results,
			      i,
			      s->current_like);

  return 0;
}

