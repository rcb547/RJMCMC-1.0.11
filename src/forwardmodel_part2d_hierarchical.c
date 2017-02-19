
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/engine.h>
#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/part2d_forwardmodel.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_debug.h>
 
static double part2d_init(void *arg);
static double part2d_init_restart(void *arg);
static int part2d_select(void *arg);
static void *part2d_perturb(void *arg, int proc);
static double part2d_misfit(void *arg, void *state);
static int part2d_accept(void *arg, 
			 double current,
			 double proposed);
static int part2d_sample(void *arg,
			 int i);

enum {
  BIRTH = 0,
  DEATH,
  MOVE,
  LOCALVALUE,
  HIERARCHICALVALUE,
  GLOBALVALUE,

  NPROCESSES,
  NPROCESSES_NOGLOBAL = GLOBALVALUE
};


struct part2dfm {

  resultset2dfm_t *results;

  part2d_forwardmodel_t *current;
  double current_like;
  double current_logdetce;

  part2d_forwardmodel_t *proposed;
  double proposed_like;
  double proposed_logdetce;

  int minpart;
  int maxpart;
  int initpart;

  bbox2d_t bound;
  bbox2d_t pbound;
  
  int nprocesses;

  int out;
  int accepted;
  int process;

  double birth_prob;
  double death_prob;
  double move_prob;
  double value_prob;
  double global_value_prob;
  double hierarchical_prob;

  int nglobalparameters;
  const forwardmodelparameter_t *global_parameters;

  int nlocalparameters;
  const forwardmodelparameter_t *local_parameters;

  int nhierarchicalparameters;
  const forwardmodelparameter_t *hierarchical_parameters;

  double proddelta_l;

  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;

  part2d_fm_hierarchical_likelihood_t lp;
  void *user_arg;

  int xsamples;
  int ysamples;

  double *mf_values;
  double *mf_gradients;

  double *x;
  double *y;
  double ***z;
  
};
  
struct _part2d_fm_likelihood_state {
  part2d_forwardmodel_t *model;

  int nlocalparameters;
  double *values;
  double *gradients;

  int nhierarchicalparameters;
  double *hierarchical_parameters;
};

static const double *
part2d_fm_likelihood_value_callback(part2d_fm_likelihood_state_t *state,
				    double x,
				    double y)
{
  if (part2d_forwardmodel_value_at(state->model,
				   x,
				   y,
				   state->nlocalparameters,
				   state->values) < 0) {
    return NULL;
  } else {
    return state->values;
  }
}

static const double *
part2d_fm_likelihood_gradient_callback(part2d_fm_likelihood_state_t *state,
				       double x,
				       double y)
{
  if (part2d_forwardmodel_gradient_at(state->model,
				      x,
				      y,
				      state->nlocalparameters,
				      state->gradients) < 0) {
    return NULL;
  } else {
    return state->gradients;
  }
}

resultset2dfm_t *
part2d_forwardmodel_hierarchical(int burnin,
				 int total,
				 int thin,
				 int minpart,
				 int maxpart,
				 int initpart,
				 double xmin,
				 double xmax,
				 double ymin,
				 double ymax,
				 int xsamples,
				 int ysamples,
				 int zsamples,
				 double credible_interval,
				 double pdx,
				 double pdy,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 int nglobalparameters,
				 const forwardmodelparameter_t *
				 global_parameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *
				 local_parameters,
				 int nhierarchicalparameters,
				 const forwardmodelparameter_t *
				 hierarchical_parameters,
				 part2d_fm_hierarchical_likelihood_t lp,
				 void *user_arg,
				 int results)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;
  int i;
  int xs;

  memset(&s, 0, sizeof(s));
  memset(&cb, 0, sizeof(cb));

  if (nlocalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nhierarchicalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "there needs to be at least on hierarchical parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "a forward model function must be provided\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = NPROCESSES_NOGLOBAL;
  } else {
    s.nprocesses = NPROCESSES;
  }
  
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.initpart = initpart;
  s.bound.xmin = xmin;
  s.bound.xmax = xmax;
  s.bound.ymin = ymin;
  s.bound.ymax = ymax;
  
  s.results = resultset2dfm_create(burnin,
				   total,
				   thin,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   nhierarchicalparameters,
				   xsamples,
				   ysamples,
				   zsamples,
				   maxpart,
				   xmin,
				   xmax,
				   ymin,
				   ymax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = part2d_forwardmodel_create(PART2D_FM_ZERO,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 ymin,
					 ymax,
					 pdx,
					 pdy,
					 nglobalparameters,
					 nlocalparameters,
					 nhierarchicalparameters,
					 0);
  if (s.current == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s.proposed = part2d_forwardmodel_create(PART2D_FM_ZERO,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  ymin,
					  ymax,
					  pdx,
					  pdy,
					  nglobalparameters,
					  nlocalparameters,
					  nhierarchicalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;

  s.nhierarchicalparameters = nhierarchicalparameters;
  s.hierarchical_parameters = hierarchical_parameters;
  
  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;
  s.ysamples = ysamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
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

  s.y = rjmcmc_create_array_1d(ysamples);
  if (s.y == NULL) {
    return NULL;
  }

  s.z = rjmcmc_create_array_3d(nlocalparameters, xsamples, ysamples);
  if (s.z == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset2dfm_fill_xcoord_vector(s.results, s.x, &xs);

  xs = ysamples;
  resultset2dfm_fill_ycoord_vector(s.results, s.y, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  cb.initialize_state = part2d_init;
  cb.select_process = part2d_select;
  cb.perturb_process = part2d_perturb;
  cb.compute_misfit = part2d_misfit;
  cb.accept = part2d_accept;
  cb.sample = part2d_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);
  
  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_1d(s.y);
  rjmcmc_destroy_array_3d(s.nlocalparameters, s.xsamples, s.z);

  part2d_forwardmodel_destroy(s.current);
  part2d_forwardmodel_destroy(s.proposed);

  resultset2dfm_assemble_results(s.results);

  return s.results;
}

#if defined(HAVE_MPI_H)

resultset2dfm_t *
MPI_part2d_forwardmodel_hierarchical(int burnin,
				     int total,
				     int thin,
				     int minpart,
				     int maxpart,
				     int initpart,
				     double xmin,
				     double xmax,
				     double ymin,
				     double ymax,
				     int xsamples,
				     int ysamples,
				     int zsamples,
				     double credible_interval,
				     double pdx,
				     double pdy,
				     rjmcmc_uniform_rand_t random,
				     rjmcmc_normal_rand_t normal,
				     int nglobalparameters,
				     const forwardmodelparameter_t *
				     global_parameters,
				     int nlocalparameters,
				     const forwardmodelparameter_t *
				     local_parameters,
				     int nhierarchicalparameters,
				     const forwardmodelparameter_t *
				     hierarchical_parameters,
				     part2d_fm_hierarchical_likelihood_t lp,
				     void *user_arg,
				     int results,
				     int mpisize,
				     int mpirank,
				     int root,
				     MPI_Comm comm)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;
  int i;
  int xs;

  if (nlocalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nhierarchicalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "there needs to be at least on hierarchical parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "a forward model function must be provided\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = NPROCESSES_NOGLOBAL;
  } else {
    s.nprocesses = NPROCESSES;
  }
  
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.initpart = initpart;
  s.bound.xmin = xmin;
  s.bound.xmax = xmax;
  s.bound.ymin = ymin;
  s.bound.ymax = ymax;
  
  s.results = resultset2dfm_create(burnin,
				   total,
				   thin,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   nhierarchicalparameters,
				   xsamples,
				   ysamples,
				   zsamples,
				   maxpart,
				   xmin,
				   xmax,
				   ymin,
				   ymax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = part2d_forwardmodel_create(PART2D_FM_ZERO,
					 minpart,
					 maxpart,
					 xmin,
					 xmax,
					 ymin,
					 ymax,
					 pdx,
					 pdy,
					 nglobalparameters,
					 nlocalparameters,
					 nhierarchicalparameters,
					 0);
  if (s.current == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create current state\n");
    return NULL;
  }

  s.proposed = part2d_forwardmodel_create(PART2D_FM_ZERO,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  ymin,
					  ymax,
					  pdx,
					  pdy,
					  nglobalparameters,
					  nlocalparameters,
					  nhierarchicalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
  s.nhierarchicalparameters = nhierarchicalparameters;
  s.hierarchical_parameters = hierarchical_parameters;

  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;
  s.ysamples = ysamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
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

  s.y = rjmcmc_create_array_1d(ysamples);
  if (s.y == NULL) {
    return NULL;
  }

  s.z = rjmcmc_create_array_3d(nlocalparameters, xsamples, ysamples);
  if (s.z == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset2dfm_fill_xcoord_vector(s.results, s.x, &xs);

  xs = ysamples;
  resultset2dfm_fill_ycoord_vector(s.results, s.y, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks
   */
  cb.initialize_state = part2d_init;
  cb.select_process = part2d_select;
  cb.perturb_process = part2d_perturb;
  cb.compute_misfit = part2d_misfit;
  cb.accept = part2d_accept;
  cb.sample = part2d_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);
  
  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_1d(s.y);
  rjmcmc_destroy_array_3d(s.nlocalparameters, s.xsamples, s.z);

  part2d_forwardmodel_destroy(s.current);
  part2d_forwardmodel_destroy(s.proposed);

  MPI_resultset2dfm_assemble_results(s.results,
				     mpisize,
				     mpirank,
				     root,
				     comm);

  return s.results;
}


resultset2dfm_t *
MPI_part2d_forwardmodel_hierarchical_restartable(const char *model_file_in_template,
						 const char *model_file_out_template,
						 int burnin,
						 int total,
						 int thin,
						 int minpart,
						 int maxpart,
						 int initpart,
						 double xmin,
						 double xmax,
						 double ymin,
						 double ymax,
						 int xsamples,
						 int ysamples,
						 int zsamples,
						 double credible_interval,
						 double pdx,
						 double pdy,
						 rjmcmc_uniform_rand_t random,
						 rjmcmc_normal_rand_t normal,
						 int nglobalparameters,
						 const forwardmodelparameter_t *
						 global_parameters,
						 int nlocalparameters,
						 const forwardmodelparameter_t *
						 local_parameters,
						 int nhierarchicalparameters,
						 const forwardmodelparameter_t *
						 hierarchical_parameters,
						 part2d_fm_hierarchical_likelihood_t lp,
						 void *user_arg,
						 int results,
						 int mpisize,
						 int mpirank,
						 int root,
						 MPI_Comm comm)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;
  int i;
  int xs;

  char filename[256];

  if (nlocalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (nhierarchicalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "there needs to be at least on hierarchical parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel_hierarchical: "
		 "a forward model function must be provided\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    s.nprocesses = NPROCESSES_NOGLOBAL;
  } else {
    s.nprocesses = NPROCESSES;
  }
  
  s.minpart = minpart;
  s.maxpart = maxpart;
  s.initpart = initpart;
  s.bound.xmin = xmin;
  s.bound.xmax = xmax;
  s.bound.ymin = ymin;
  s.bound.ymax = ymax;
  
  s.results = resultset2dfm_create(burnin,
				   total,
				   thin,
				   nglobalparameters,
				   global_parameters,
				   nlocalparameters,
				   local_parameters,
				   nhierarchicalparameters,
				   xsamples,
				   ysamples,
				   zsamples,
				   maxpart,
				   xmin,
				   xmax,
				   ymin,
				   ymax,
				   s.nprocesses,
				   credible_interval,
				   results);
  if (s.results == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create results\n");
    return NULL;
  }

  if (model_file_in_template == NULL ||
      strlen(model_file_in_template) == 0) {
    s.current = part2d_forwardmodel_create(PART2D_FM_ZERO,
					   minpart,
					   maxpart,
					   xmin,
					   xmax,
					   ymin,
					   ymax,
					   pdx,
					   pdy,
					   nglobalparameters,
					   nlocalparameters,
					   nhierarchicalparameters,
					   0);
    if (s.current == NULL) {
      rjmcmc_error("part2d_forwardmodel: failed to create current state\n");
      return NULL;
    }

    cb.initialize_state = part2d_init;

  } else {
    sprintf(filename, model_file_in_template, mpirank);
    s.current = part2d_forwardmodel_load(filename);
    if (s.current == NULL) {
      rjmcmc_error("part2d_forwardmodel: failed to load model from file `%s'\n", 
		   filename);
      return NULL;
    }

    cb.initialize_state = part2d_init_restart;
  }


  s.proposed = part2d_forwardmodel_create(PART2D_FM_ZERO,
					  minpart,
					  maxpart,
					  xmin,
					  xmax,
					  ymin,
					  ymax,
					  pdx,
					  pdy,
					  nglobalparameters,
					  nlocalparameters,
					  nhierarchicalparameters,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
  s.nhierarchicalparameters = nhierarchicalparameters;
  s.hierarchical_parameters = hierarchical_parameters;

  s.random = random;
  s.normal = normal;

  s.xsamples = xsamples;
  s.ysamples = ysamples;

  /*
   * Create temporary arrays for passing to the misfit function
   */
  s.mf_values = rjmcmc_create_array_1d(nlocalparameters);
  if (s.mf_values == NULL) {
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

  s.y = rjmcmc_create_array_1d(ysamples);
  if (s.y == NULL) {
    return NULL;
  }

  s.z = rjmcmc_create_array_3d(nlocalparameters, xsamples, ysamples);
  if (s.z == NULL) {
    return NULL;
  }

  xs = xsamples;
  resultset2dfm_fill_xcoord_vector(s.results, s.x, &xs);

  xs = ysamples;
  resultset2dfm_fill_ycoord_vector(s.results, s.y, &xs);

  s.proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s.proddelta_l *= (local_parameters[i].fmax - local_parameters[i].fmin);
  }

  s.lp = lp;
  s.user_arg = user_arg;

  /*
   * Set the engine callbacks (note init callback set when checking restart above).
   */
  cb.select_process = part2d_select;
  cb.perturb_process = part2d_perturb;
  cb.compute_misfit = part2d_misfit;
  cb.accept = part2d_accept;
  cb.sample = part2d_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.mf_values);
  rjmcmc_destroy_array_1d(s.mf_gradients);
  
  rjmcmc_destroy_array_1d(s.x);
  rjmcmc_destroy_array_1d(s.y);
  rjmcmc_destroy_array_3d(s.nlocalparameters, s.xsamples, s.z);

  part2d_forwardmodel_destroy(s.proposed);

  MPI_resultset2dfm_assemble_results(s.results,
				     mpisize,
				     mpirank,
				     root,
				     comm);

  sprintf(filename, model_file_out_template, mpirank);
  if (part2d_forwardmodel_save(s.current, filename) < 0) {
    fprintf(stderr, "Failed to save final model\n");
    return NULL;
  }

  part2d_forwardmodel_destroy(s.current);

  return s.results;
}

#endif /* HAVE_MPI_H */

static double part2d_init(void *arg)
{
  struct part2dfm *s = (struct part2dfm *)arg;
  part2d_fm_likelihood_state_t state;

  if (part2d_forwardmodel_initialize(s->current,
				     s->global_parameters,
				     s->nglobalparameters,
				     s->local_parameters,
				     s->nlocalparameters,
				     s->hierarchical_parameters,
				     s->nhierarchicalparameters,
				     s->initpart,
				     s->random,
				     s->normal) < 0) {
    rjmcmc_error("part2d_init: failed to initialize\n");
    return -1.0;
  }

  state.model = s->current;
  state.nlocalparameters = s->nlocalparameters;
  state.values = s->mf_values;
  state.gradients = s->mf_gradients;
  s->current_like = s->lp(s->user_arg,
			  s->nglobalparameters,
			  part2d_forwardmodel_global_parameters(s->current),
			  1, /* Require calculation of logdetce */
			  s->nhierarchicalparameters,
			  part2d_forwardmodel_hierarchical_parameters(s->current),
			  &state,
			  part2d_fm_likelihood_value_callback,
			  part2d_fm_likelihood_gradient_callback,
			  &(s->bound),
			  &(s->current_logdetce));

  return s->current_like;
 
}

static double part2d_init_restart(void *arg)
{
  struct part2dfm *s = (struct part2dfm *)arg;
  part2d_fm_likelihood_state_t state;

  state.model = s->current;
  state.nlocalparameters = s->nlocalparameters;
  state.values = s->mf_values;
  state.gradients = s->mf_gradients;
  s->current_like = s->lp(s->user_arg,
			  s->nglobalparameters,
			  part2d_forwardmodel_global_parameters(s->current),
			  1, /* Require calculation of logdetce */
			  s->nhierarchicalparameters,
			  part2d_forwardmodel_hierarchical_parameters(s->current),
			  &state,
			  part2d_fm_likelihood_value_callback,
			  part2d_fm_likelihood_gradient_callback,
			  &(s->bound),
			  &(s->current_logdetce));

  return s->current_like;
}

static int part2d_select(void *arg)
{
  struct part2dfm *s = (struct part2dfm *)arg;

  return rjmcmc_random_choose_int(0, s->nprocesses - 1, s->random);
}

static void *part2d_perturb(void *arg, int proc)
{
  struct part2dfm *s = (struct part2dfm *)arg;

  s->process = proc;

  resultset2dfm_propose(s->results, proc);

  switch(proc) {
  case BIRTH: 
    s->out = part2d_forwardmodel_propose_birth(s->current,
					       s->proposed,
					       s->nglobalparameters,
					       s->global_parameters,
					       s->nlocalparameters,
					       s->local_parameters,
					       s->random,
					       s->normal,
					       &(s->pbound),
					       &(s->birth_prob));
    break;

  case DEATH: 
    s->out = part2d_forwardmodel_propose_death(s->current,
					       s->proposed,
					       s->nglobalparameters,
					       s->global_parameters,
					       s->nlocalparameters,
					       s->local_parameters,
					       s->random,
					       s->normal,
					       &(s->pbound),
					       &(s->death_prob));
    break;

  case MOVE: 
    s->out = part2d_forwardmodel_propose_move(s->current,
					      s->proposed,
					      s->nglobalparameters,
					      s->global_parameters,
					      s->nlocalparameters,
					      s->local_parameters,
					      s->random,
					      s->normal,
					      &(s->pbound),
					      &(s->move_prob));
    break;
    
  case LOCALVALUE: 
    s->out = part2d_forwardmodel_propose_local_value(s->current,
						     s->proposed,
						     s->nglobalparameters,
						     s->global_parameters,
						     s->nlocalparameters,
						     s->local_parameters,
						     s->random,
						     s->normal,
						     &(s->pbound),
						     &(s->value_prob));
    break;

  case GLOBALVALUE:
    s->out = part2d_forwardmodel_propose_global_value(s->current,
						      s->proposed,
						      s->nglobalparameters,
						      s->global_parameters,
						      s->nlocalparameters,
						      s->local_parameters,
						      s->random,
						      s->normal,
						      &(s->value_prob));
    break;

  case HIERARCHICALVALUE:
    s->out = part2d_forwardmodel_propose_hierarchical(s->current,
						      s->proposed,
						      s->nhierarchicalparameters,
						      s->hierarchical_parameters,
						      s->random,
						      s->normal,
						      &(s->hierarchical_prob));
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

static double part2d_misfit(void *arg, void *_state)
{
  struct part2dfm *s = (struct part2dfm *)arg;
  part2d_fm_likelihood_state_t state;
  int hierarchical = 0;

  if (s->out) {

    if (s->process == HIERARCHICALVALUE) {
      hierarchical = 1;
    }

    state.model = s->proposed;
    state.nlocalparameters = s->nlocalparameters;
    state.values = s->mf_values;
    state.gradients = s->mf_gradients;
    s->proposed_like = s->lp(s->user_arg,
			     s->nglobalparameters,
			     part2d_forwardmodel_global_parameters(s->proposed),
			     hierarchical,
			     s->nhierarchicalparameters,
			     part2d_forwardmodel_hierarchical_parameters(s->proposed),
			     &state,
			     part2d_fm_likelihood_value_callback,
			     part2d_fm_likelihood_gradient_callback,
			     &(s->bound),
			     &(s->proposed_logdetce));
  } else {
    s->proposed_like = FLT_MAX;
  }

  return s->proposed_like;
}

static int part2d_accept(void *arg, 
			 double current,
			 double proposed)
{
  struct part2dfm *s = (struct part2dfm *)arg;
  double u;
  
  s->accepted = 0;

  if (s->out == 0) {
    return 0;
  }

  u = log(s->random());

  switch (s->process) {
  case BIRTH: 
    s->accepted = (u < (current - proposed - 
			log(s->proddelta_l) - log(s->birth_prob)));
    break;

  case DEATH:
    s->accepted = (u < (current - proposed + 
			log(s->proddelta_l) + log(s->death_prob)));
    break;

  case MOVE:
    s->accepted = (u < (current - proposed));
    break;

  case LOCALVALUE:
    s->accepted = (u < (current - proposed));
    break;

  case GLOBALVALUE:
    s->accepted = (u < (current - proposed));
    break;

  case HIERARCHICALVALUE:
    s->accepted = (u < (0.5 * (s->current_logdetce - s->proposed_logdetce) + 
			current - proposed));
    break;

  default:
    rjmcmc_error("part2d_accept: invalid process %d\n", s->process);
    s->accepted = 0;
    break;
  }

  if (s->accepted) {
    part2d_forwardmodel_clone(s->proposed, s->current);
    s->current_like = s->proposed_like;
    s->current_logdetce = s->proposed_logdetce;

    resultset2dfm_accept(s->results, s->process);
  }

  return s->accepted;
}



static int part2d_sample(void *arg,
			 int i)
{
  int gi;
  int hi;
  int li;
  int npartitions;
  struct part2dfm *s = (struct part2dfm *)arg;

  const double *gp;
  const double *hp;

  int j;
  double x;
  double y;

  resultset2dfm_sample_misfit(s->results, i, s->current_like);

  gp = part2d_forwardmodel_global_parameters(s->current);
  if (gp == NULL && s->nglobalparameters > 0) {
    return -1;
  }

  for (gi = 0; gi < s->nglobalparameters; gi ++) {
    resultset2dfm_sample_global_parameter(s->results, 
					  i, 
					  gi, 
					  gp[gi]);
  }

  hp = part2d_forwardmodel_hierarchical_parameters(s->current);
  if (hp == NULL) {
    return -1;
  }

  for (hi = 0; hi < s->nhierarchicalparameters; hi ++) {
    resultset2dfm_sample_hierarchical_parameter(s->results,
						i,
						hi,
						hp[hi]);
  }

  if (part2d_forwardmodel_evaluate_local_parameters(s->current,
						    s->xsamples,
						    s->x,
						    s->ysamples,
						    s->y,
						    s->nlocalparameters,
						    s->z) < 0) {
    return -1;
  }

  for (li = 0; li < s->nlocalparameters; li ++) {
    resultset2dfm_sample_local_parameter(s->results,
					 i,
					 li,
					 s->z[li]);
  }

  npartitions = part2d_forwardmodel_partitions(s->current);
  if (npartitions < 0) {
    return -1;
  }

  resultset2dfm_sample_npartitions(s->results,
				   i,
				   npartitions);

  for (j = 0; j < npartitions; j ++) {
    if (part2d_forwardmodel_partition_centre(s->current, j, &x, &y) >= 0) {
      resultset2dfm_sample_centre(s->results,
				x, y);
    } else {
      rjmcmc_error("part2d_sample: Failed to get centre\n");
    }
  }

  return 0;
}

