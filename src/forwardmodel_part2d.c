
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/engine.h>
#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/forwardmodel_mpi.h>
#include <rjmcmc/part2d_forwardmodel.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_defines.h>
#include <rjmcmc/rjmcmc_debug.h>
 
#define RESULTSET_CHECKPOINT_TEMPLATE "%s_%03d_%d.resultset.rjmcmc"
#define MODEL_CHECKPOINT_TEMPLATE "%s_%03d_%d.model.rjmcmc"

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
  GLOBALVALUE,

  NPROCESSES,
  NPROCESSES_NOGLOBAL = GLOBALVALUE
};

struct part2dfm {

  resultset2dfm_t *results;

  part2d_forwardmodel_t *current;
  double current_like;

  part2d_forwardmodel_t *proposed;
  double proposed_like;

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

  int nglobalparameters;
  const forwardmodelparameter_t *global_parameters;

  int nlocalparameters;
  const forwardmodelparameter_t *local_parameters;

  double proddelta_l;

  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;

  part2d_fm_likelihood_t lp;
  void *user_arg;

  int xsamples;
  int ysamples;

  double *mf_global_parameters;

  double *mf_values;
  double *mf_gradients;

  double *x;
  double *y;
  double ***z;

  int checkpoint_interval;
  const char *checkpoint_prefix;

  int mpirank;
};
  
static int common_initialization(struct part2dfm *s,
				 int burnin,
				 int total,
				 int thin,
				 int min_partitions,
				 int max_partitions,
				 int init_partitions,
				 double xmin,
				 double xmax,
				 double ymin,
				 double ymax,
				 double pdx,
				 double pdy,
				 int xsamples,
				 int ysamples,
				 int zsamples,
				 int nprocesses,
				 int nglobalparameters,
				 const forwardmodelparameter_t *globalparameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *localparameters,
				 double credible_interval,
				 int requested_results,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 part2d_fm_likelihood_t lp,
				 void *user_arg);

static int restart_initialization(struct part2dfm *s,
				  int *results,
				  int *burnin,
				  int *total,
				  int *thin,
				  int *xsamples,
				  int *ysamples,
				  int *zsamples,
				  int *nprocesses,
				  int *maxpartitions,
				  double *xmin,
				  double *xmax,
				  double *ymin,
				  double *ymax,
				  double *credibleinterval,
				  int *nhierarchical,
				  int *nglobal,
				  int *nlocal,
				  int mpirank,
				  int checkpoint,
				  const char *checkpoint_prefix,
				  forwardmodelparameter_t **global_parameters,
				  forwardmodelparameter_t **local_parameters,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  part2d_fm_likelihood_t lp,
				  void *user_arg);

static int common_finalization(struct part2dfm *s);


struct _part2d_fm_likelihood_state {
  part2d_forwardmodel_t *model;

  int nlocalparameters;
  double *values;
  double *gradients;
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
part2d_forwardmodel(int burnin,
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
		    const forwardmodelparameter_t *global_parameters,
		    int nlocalparameters,
		    const forwardmodelparameter_t *local_parameters,
		    part2d_fm_likelihood_t lp,
		    void *user_arg,
		    int results)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;
  int nprocesses;

  if (nlocalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel: "
		 "a forward model function must be provided\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    nprocesses = NPROCESSES_NOGLOBAL;
  } else {
    nprocesses = NPROCESSES;
  }

  if (common_initialization(&s,
			    burnin,
			    total,
			    thin,
			    minpart,
			    maxpart,
			    initpart,
			    xmin,
			    xmax,
			    ymin,
			    ymax,
			    pdx,
			    pdy,
			    xsamples,
			    ysamples,
			    zsamples,
			    nprocesses,
			    nglobalparameters,
			    global_parameters,
			    nlocalparameters,
			    local_parameters,
			    credible_interval,
			    results,
			    random,
			    normal,
			    lp,
			    user_arg) < 0) {
    return NULL;
  }

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

  if (common_finalization(&s) < 0) {
    return NULL;
  }

  resultset2dfm_assemble_results(s.results);

  return s.results;
}

resultset2dfm_t *
part2d_forwardmodel_restartable(const char *model_file_in,
				const char *model_file_out,
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
				part2d_fm_likelihood_t lp,
				void *user_arg,
				int results)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;
  int i;
  int xs;

  char filename[256];

  if (nlocalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel_restartable: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel_restartable: "
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
				   0,
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

  if (model_file_in == NULL ||
      strlen(model_file_in) == 0) {
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
					   0,
					   0);
    if (s.current == NULL) {
      rjmcmc_error("part2d_forwardmodel: failed to create current state\n");
      return NULL;
    }

    cb.initialize_state = part2d_init;

  } else {
    s.current = part2d_forwardmodel_load(model_file_in);
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
					  0,
					  0);
  if (s.proposed == NULL) {
    rjmcmc_error("part2d_forwardmodel: failed to create proposed state\n");
    return NULL;
  }
	
  s.nglobalparameters = nglobalparameters;
  s.global_parameters = global_parameters;
  
  s.nlocalparameters = nlocalparameters;
  s.local_parameters = local_parameters;
  
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

  s.checkpoint_interval = 0;
  s.checkpoint_prefix = NULL;
  s.mpirank = 0;
  
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

  resultset2dfm_assemble_results(s.results);

  if (model_file_out != NULL) {
    if (part2d_forwardmodel_save(s.current, model_file_out) < 0) {
      fprintf(stderr, "Failed to save final model\n");
      return NULL;
    }
  }

  part2d_forwardmodel_destroy(s.current);

  return s.results;
}


#if defined(HAVE_MPI_H)
resultset2dfm_t *
MPI_part2d_forwardmodel(int burnin,
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
			const forwardmodelparameter_t *global_parameters,
			int nlocalparameters,
			const forwardmodelparameter_t *local_parameters,
			part2d_fm_likelihood_t lp,
			void *user_arg,
			int results,
			int mpisize,
			int mpirank,
			int root,
			MPI_Comm comm)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;
  int nprocesses;

  if (nlocalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel: "
		 "a forward model function must be provided\n");
    return NULL;
  }

  if (nglobalparameters == 0) {
    nprocesses = NPROCESSES_NOGLOBAL;
  } else {
    nprocesses = NPROCESSES;
  }

  
  if (common_initialization(&s,
			    burnin,
			    total,
			    thin,
			    minpart,
			    maxpart,
			    initpart,
			    xmin,
			    xmax,
			    ymin,
			    ymax,
			    pdx,
			    pdy,
			    xsamples,
			    ysamples,
			    zsamples,
			    nprocesses,
			    nglobalparameters,
			    global_parameters,
			    nlocalparameters,
			    local_parameters,
			    credible_interval,
			    results,
			    random,
			    normal,
			    lp,
			    user_arg) < 0) {
    return NULL;
  }

  s.mpirank = mpirank;

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

  if (common_finalization(&s) < 0) {
    return NULL;
  }

  MPI_resultset2dfm_assemble_results(s.results,
				     mpisize,
				     mpirank,
				     root,
				     comm);

  return s.results;
}

resultset2dfm_t *
MPI_part2d_forwardmodel_checkpoint(int burnin,
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
				   const forwardmodelparameter_t *global_parameters,
				   int nlocalparameters,
				   const forwardmodelparameter_t *local_parameters,
				   part2d_fm_likelihood_t lp,
				   void *user_arg,
				   int results,
				   int checkpoint_interval,
				   const char *checkpoint_prefix,
				   int mpisize,
				   int mpirank,
				   int root,
				   MPI_Comm comm)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;
  int nprocesses;

  char filename[256];
  FILE *fp;

  if (nlocalparameters <= 0) {
    rjmcmc_error("part2d_forwardmodel: "
		 "there needs to be at least one local parameter\n");
    return NULL;
  }

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel: "
		 "a forward model function must be provided\n");
    return NULL;
  }
  
  if (checkpoint_interval <= 0) {
    rjmcmc_error("MPI_part2d_forwradmodel_checkpoint: checkpoint interval must be greater than 0\n");
    return NULL;
  }

  if (checkpoint_prefix == NULL) {
    rjmcmc_error("MPI_part2d_forwardmodel_checkpoint: checkpoint prefix cannot be NULL\n");
    return NULL;
  }

  if (strlen(checkpoint_prefix) > 128) {
    rjmcmc_error("MPI_part2d_forwardmodel_checkpoint: checkpoint prefix too long (limit 128)\n");
    return NULL;
  }

  /*
   * Check if we can create a file at the prefix and fail immediately.
   */
  sprintf(filename, "%s_%03d_%d.rjmcmc", checkpoint_prefix, 0, 0);
  fp = fopen(filename, "w");
  if (fp == NULL) {
    rjmcmc_error("MPI_part2d_forwardmodel_checkpoint: failed to create %s, verify checkpoint prefix.\n", filename);
    return NULL;
  }
  fclose(fp);
  remove(filename);

  if (nglobalparameters == 0) {
    nprocesses = NPROCESSES_NOGLOBAL;
  } else {
    nprocesses = NPROCESSES;
  }
  
  if (common_initialization(&s,
			    burnin,
			    total,
			    thin,
			    minpart,
			    maxpart,
			    initpart,
			    xmin,
			    xmax,
			    ymin,
			    ymax,
			    pdx,
			    pdy,
			    xsamples,
			    ysamples,
			    zsamples,
			    nprocesses,
			    nglobalparameters,
			    global_parameters,
			    nlocalparameters,
			    local_parameters,
			    credible_interval,
			    results,
			    random,
			    normal,
			    lp,
			    user_arg) < 0) {
    return NULL;
  }

  s.checkpoint_interval = checkpoint_interval;
  s.checkpoint_prefix = checkpoint_prefix;

  s.mpirank = mpirank;

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

  if (common_finalization(&s) < 0) {
    return NULL;
  }

  MPI_resultset2dfm_assemble_results(s.results,
				     mpisize,
				     mpirank,
				     root,
				     comm);

  return s.results;
}

resultset2dfm_t *
MPI_part2d_forwardmodel_checkpoint_restart(rjmcmc_uniform_rand_t random,
					   rjmcmc_normal_rand_t normal,
					   part2d_fm_likelihood_t lp,
					   void *user_arg,
					   int checkpoint,
					   int checkpoint_interval,
					   const char *checkpoint_prefix,
					   int mpisize,
					   int mpirank,
					   int root,
					   MPI_Comm comm,
					   int *requested_results,
					   int *burnin,
					   int *total,
					   int *thin,
					   int *xsamples,
					   int *ysamples,
					   int *zsamples,
					   int *nprocesses,
					   int *maxpartitions,
					   double *xmin,
					   double *xmax,
					   double *ymin,
					   double *ymax,
					   double *credibleinterval,
					   int *nhierarchical,
					   int *nglobal,
					   int *nlocal,
					   forwardmodelparameter_t **local_parameters,
					   forwardmodelparameter_t **global_parameters)
{
  struct part2dfm s;
  rjmcmc_engine_cb_t cb;

  if (lp == NULL) {
    rjmcmc_error("part2d_forwardmodel: "
		 "a forward model function must be provided\n");
    return NULL;
  }
  
  if (checkpoint_interval <= 0) {
    rjmcmc_error("MPI_part2d_forwradmodel_checkpoint: checkpoint interval must be greater than 0\n");
    return NULL;
  }

  if (checkpoint_prefix == NULL) {
    rjmcmc_error("MPI_part2d_forwardmodel_checkpoint: checkpoint prefix cannot be NULL\n");
    return NULL;
  }

  if (strlen(checkpoint_prefix) > 128) {
    rjmcmc_error("MPI_part2d_forwardmodel_checkpoint: checkpoint prefix too long (limit 128)\n");
    return NULL;
  }

  if (restart_initialization(&s,
			     requested_results,
			     burnin,
			     total,
			     thin,
			     xsamples,
			     ysamples,
			     zsamples,
			     nprocesses,
			     maxpartitions,
			     xmin,
			     xmax,
			     ymin,
			     ymax,
			     credibleinterval,
			     nhierarchical,
			     nglobal,
			     nlocal,
			     mpirank,
			     checkpoint,
			     checkpoint_prefix,
			     local_parameters,
			     global_parameters,
			     random,
			     normal,
			     lp,
			     user_arg) < 0) {
    return NULL;
  }

  s.checkpoint_interval = checkpoint_interval;
  s.checkpoint_prefix = checkpoint_prefix;

  s.mpirank = mpirank;

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

  if (rjmcmc_engine_restart(&cb,
			    *burnin,
			    *total,
			    1,
			    checkpoint,
			    s.current_like) < 0) {
    return NULL;
  }  

  if (common_finalization(&s) < 0) {
    return NULL;
  }

  MPI_resultset2dfm_assemble_results(s.results,
				     mpisize,
				     mpirank,
				     root,
				     comm);

  return s.results;
}


#endif /* HAVE_MPI_H */

static int common_initialization(struct part2dfm *s,
				 int burnin,
				 int total,
				 int thin,
				 int min_partitions,
				 int max_partitions,
				 int init_partitions,
				 double xmin,
				 double xmax,
				 double ymin,
				 double ymax,
				 double pdx,
				 double pdy,
				 int xsamples,
				 int ysamples,
				 int zsamples,
				 int nprocesses,
				 int nglobalparameters,
				 const forwardmodelparameter_t *globalparameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *localparameters,
				 double credible_interval,
				 int requested_results,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 part2d_fm_likelihood_t lp,
				 void *user_arg)
{
  int xs;
  int ys;
  int i;

  s->results = NULL;

  s->current = NULL;
  s->current_like = 0.0;

  s->proposed = NULL;
  s->proposed_like = 0.0;

  s->minpart = min_partitions;
  s->maxpart = max_partitions;
  s->initpart = init_partitions;

  s->bound.xmin = xmin;
  s->bound.xmax = xmax;
  s->bound.ymin = ymin;
  s->bound.ymax = ymax;
  
  s->pbound.xmin = xmin;
  s->pbound.xmax = xmax;
  s->pbound.ymin = ymin;
  s->pbound.ymax = ymax;
  
  s->nprocesses = nprocesses;

  s->out = 0;
  s->accepted = 0;
  s->process = 0;

  s->birth_prob = 0.0;
  s->death_prob = 0.0;
  s->move_prob = 0.0;
  s->value_prob = 0.0;
  s->global_value_prob = 0.0;

  s->nglobalparameters = nglobalparameters;
  s->global_parameters = globalparameters;

  s->nlocalparameters = nlocalparameters;
  s->local_parameters = localparameters;

  s->proddelta_l = 0.0;

  s->random = random;
  s->normal = normal;

  s->lp = lp;
  s->user_arg = user_arg;

  s->xsamples = xsamples;
  s->ysamples = ysamples;

  s->results = resultset2dfm_create(burnin,
				    total,
				    thin,
				    nglobalparameters,
				    globalparameters,
				    nlocalparameters,
				    localparameters,
				    0, /* nhierarchicalparameters */
				    xsamples,
				    ysamples,
				    zsamples,
				    max_partitions,
				    xmin,
				    xmax,
				    ymin,
				    ymax,
				    nprocesses,
				    credible_interval,
				    requested_results);
  RJMCMC_NULLCHECKINT(s->results, "common_initialization: NULL resultset.\n");

  s->current = part2d_forwardmodel_create(PART2D_FM_ZERO,
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
					  0 /*nhierarchicalparameters*/,
					  0);
  RJMCMC_NULLCHECKINT(s->current, "common_initialization: NULL current model.\n");
  
  s->proposed = part2d_forwardmodel_create(PART2D_FM_ZERO,
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
					   0 /*nhierarchicalparameters*/,
					   0);
  RJMCMC_NULLCHECKINT(s->proposed, "common_initialization: NULL proposed model.\n");
  
  
  s->checkpoint_interval = 0;
  s->checkpoint_prefix = NULL;

  s->mf_global_parameters = NULL;
  if (nglobalparameters > 0) {
    s->mf_global_parameters = rjmcmc_create_array_1d(nglobalparameters);
    RJMCMC_NULLCHECKINT(s->mf_global_parameters, "common_initialization: NULL global array\n");
  }

  /*
   * Create temporary arrays for passing to the misfit function
   */

  s->mf_values = rjmcmc_create_array_1d(nlocalparameters);
  RJMCMC_NULLCHECKINT(s->mf_values, "common_initialization: NULL values array\n");
  
  s->mf_gradients = rjmcmc_create_array_1d(nlocalparameters);
  RJMCMC_NULLCHECKINT(s->mf_gradients, "common_initialization: NULL gradient array\n");

  s->x = rjmcmc_create_array_1d(xsamples);
  RJMCMC_NULLCHECKINT(s->x, "common_initialization: NULL x array\n");

  s->y = rjmcmc_create_array_1d(ysamples);
  RJMCMC_NULLCHECKINT(s->y, "common_initialization: NULL y array\n");

  s->z = rjmcmc_create_array_3d(nlocalparameters, xsamples, ysamples);
  RJMCMC_NULLCHECKINT(s->z, "common_initialization: NULL z array\n");

  xs = xsamples;
  resultset2dfm_fill_xcoord_vector(s->results, s->x, &xs);

  ys = ysamples;
  resultset2dfm_fill_ycoord_vector(s->results, s->y, &ys);

  s->proddelta_l = 1.0;
  for (i = 0; i < nlocalparameters; i ++) {
    s->proddelta_l *= (localparameters[i].fmax - localparameters[i].fmin);
  }

  s->mpirank = 0;

  return 0;
}

static int restart_initialization(struct part2dfm *s,
				  int *requested_results,
				  int *burnin,
				  int *total,
				  int *thin,
				  int *xsamples,
				  int *ysamples,
				  int *zsamples,
				  int *nprocesses,
				  int *maxpartitions,
				  double *xmin,
				  double *xmax,
				  double *ymin,
				  double *ymax,
				  double *credibleinterval,
				  int *nhierarchical,
				  int *nglobal,
				  int *nlocal,
				  int mpirank,
				  int checkpoint,
				  const char *checkpoint_prefix,
				  forwardmodelparameter_t **global_parameters,
				  forwardmodelparameter_t **local_parameters,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  part2d_fm_likelihood_t lp,
				  void *user_arg)
{
  char filename[256];

  int xs;
  int ys;

  int i;

  const double *misfit;

  sprintf(filename, RESULTSET_CHECKPOINT_TEMPLATE, checkpoint_prefix, mpirank, checkpoint);
  s->results = resultset2dfm_load(filename,
				  RESULTSET2DFM_BINARY,
				  requested_results,
				  burnin,
				  total,
				  thin,
				  xsamples,
				  ysamples,
				  zsamples,
				  nprocesses,
				  maxpartitions,
				  xmin,
				  xmax,
				  ymin,
				  ymax,
				  credibleinterval,
				  nhierarchical,
				  nglobal,
				  nlocal,
				  global_parameters,
				  local_parameters);
  RJMCMC_NULLCHECKINT(s->results, "restart_initialization: failed to load resultset.\n");

  sprintf(filename, MODEL_CHECKPOINT_TEMPLATE, checkpoint_prefix, mpirank, checkpoint);
  s->current = part2d_forwardmodel_load(filename);
  RJMCMC_NULLCHECKINT(s->current, "restart_initialization: failed to load model.\n");

  s->proposed = part2d_forwardmodel_create(PART2D_FM_ZERO,
					   part2d_forwardmodel_min_partitions(s->current),
					   *maxpartitions,
					   *xmin,
					   *xmax,
					   *ymin,
					   *ymax,
					   part2d_forwardmodel_pdx(s->current),
					   part2d_forwardmodel_pdy(s->current),
					   *nglobal,
					   *nlocal,
					   0 /*nhierarchicalparameters*/,
					   0 /*include_corners*/);

  misfit = resultset2dfm_get_misfit(s->results);
  RJMCMC_NULLCHECKINT(misfit, "restart_initialization: failed to retrieve misfit data from checkpoint\n");
  s->current_like = misfit[checkpoint - 1];
  if (s->current_like <= 0.0) {
    rjmcmc_error("restart_initialization: 0 likelihood\n");
    return -1;
  }
  s->proposed_like = 0.0;

  s->minpart = 1;
  s->maxpart = *maxpartitions;

  s->bound.xmin = *xmin;
  s->bound.xmax = *xmax;
  s->bound.ymin = *ymin;
  s->bound.ymax = *ymax;
  
  s->pbound.xmin = *xmin;
  s->pbound.xmax = *xmax;
  s->pbound.ymin = *ymin;
  s->pbound.ymax = *ymax;
  
  s->nprocesses = *nprocesses;

  s->out = 0;
  s->accepted = 0;
  s->process = 0;

  s->birth_prob = 0.0;
  s->death_prob = 0.0;
  s->move_prob = 0.0;
  s->value_prob = 0.0;
  s->global_value_prob = 0.0;

  s->nglobalparameters = *nglobal;
  s->global_parameters = *global_parameters;

  s->nlocalparameters = *nlocal;
  s->local_parameters = *local_parameters;

  s->proddelta_l = 0.0;

  s->random = random;
  s->normal = normal;

  s->lp = lp;
  s->user_arg = user_arg;

  s->xsamples = *xsamples;
  s->ysamples = *ysamples;

  s->checkpoint_interval = 0;
  s->checkpoint_prefix = NULL;

  s->mf_global_parameters = NULL;
  if ((*nglobal) > 0) {
    s->mf_global_parameters = rjmcmc_create_array_1d((*nglobal));
    RJMCMC_NULLCHECKINT(s->mf_global_parameters, "common_initialization: NULL global array\n");
  }

  /*
   * Create temporary arrays for passing to the misfit function
   */

  s->mf_values = rjmcmc_create_array_1d((*nlocal));
  RJMCMC_NULLCHECKINT(s->mf_values, "common_initialization: NULL values array\n");
  
  s->mf_gradients = rjmcmc_create_array_1d((*nlocal));
  RJMCMC_NULLCHECKINT(s->mf_gradients, "common_initialization: NULL gradient array\n");

  s->x = rjmcmc_create_array_1d(*xsamples);
  RJMCMC_NULLCHECKINT(s->x, "common_initialization: NULL x array\n");

  s->y = rjmcmc_create_array_1d(*ysamples);
  RJMCMC_NULLCHECKINT(s->y, "common_initialization: NULL y array\n");

  s->z = rjmcmc_create_array_3d((*nlocal), *xsamples, *ysamples);
  RJMCMC_NULLCHECKINT(s->z, "common_initialization: NULL z array\n");

  xs = *xsamples;
  resultset2dfm_fill_xcoord_vector(s->results, s->x, &xs);

  ys = *ysamples;
  resultset2dfm_fill_ycoord_vector(s->results, s->y, &ys);

  s->proddelta_l = 1.0;
  for (i = 0; i < (*nlocal); i ++) {
    s->proddelta_l *= ((*local_parameters)[i].fmax - (*local_parameters)[i].fmin);
  }

  s->mpirank = mpirank;


  return 0;
}




static int common_finalization(struct part2dfm *s)
{
  rjmcmc_destroy_array_1d(s->mf_global_parameters);
  rjmcmc_destroy_array_1d(s->mf_values);
  rjmcmc_destroy_array_1d(s->mf_gradients);
  
  rjmcmc_destroy_array_1d(s->x);
  rjmcmc_destroy_array_1d(s->y);
  rjmcmc_destroy_array_3d(s->nlocalparameters, s->xsamples, s->z);

  part2d_forwardmodel_destroy(s->current);
  part2d_forwardmodel_destroy(s->proposed);

  return 0;
}

static double part2d_init(void *arg)
{
  struct part2dfm *s = (struct part2dfm *)arg;
  part2d_fm_likelihood_state_t state;

  if (part2d_forwardmodel_initialize(s->current,
				     s->global_parameters,
				     s->nglobalparameters,
				     s->local_parameters,
				     s->nlocalparameters,
				     NULL, /*hierarchical_parameters*/
				     0, /*nhierarchicalparameters*/
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
			  s->mf_global_parameters,
			  &state,
			  part2d_fm_likelihood_value_callback,
			  part2d_fm_likelihood_gradient_callback,
			  &(s->bound));

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
			  s->mf_global_parameters,
			  &state,
			  part2d_fm_likelihood_value_callback,
			  part2d_fm_likelihood_gradient_callback,
			  &(s->bound));

  return s->current_like;
}


static int part2d_select(void *arg)
{
  struct part2dfm *s = (struct part2dfm *)arg;
  double u = s->random();

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

  if (s->out) {
    state.model = s->proposed;
    state.nlocalparameters = s->nlocalparameters;
    state.values = s->mf_values;
    state.gradients = s->mf_gradients;
    s->proposed_like = s->lp(s->user_arg,
			     s->nglobalparameters,
			     s->mf_global_parameters,
			     &state,
			     part2d_fm_likelihood_value_callback,
			     part2d_fm_likelihood_gradient_callback,
			     &(s->bound));
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
    /* printf("%10.5f %10.5f: %10.5f %10.5f %10.5f %10.5f\n",  */
    /* 	   u, */
    /* 	   (current - proposed -  */
    /* 	    log(s->proddelta_l) - log(s->birth_prob)), */
    /* 	   current, proposed, log(s->proddelta_l), log(s->birth_prob)); */
    if (u < (current - proposed - 
	     log(s->proddelta_l) - log(s->birth_prob))) {
      s->accepted = 1;
    }
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

  default:
    rjmcmc_error("part2d_accept: invalid process %d\n", s->process);
    s->accepted = 0;
    break;
  }

  if (s->accepted) {
    part2d_forwardmodel_clone(s->proposed, s->current);
    s->current_like = s->proposed_like;

    resultset2dfm_accept(s->results, s->process);
  }

  return s->accepted;
}



static int part2d_sample(void *arg,
			 int i)
{
  int gi;
  int li;
  int npartitions;
  struct part2dfm *s = (struct part2dfm *)arg;

  const double *gp;

  int j;
  double x;
  double y;

  char filename[256];

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

  if ((s->checkpoint_interval > 0) &&
      ((i % s->checkpoint_interval) == 0)) {

    /*
     * Write a checkpoint.
     */
    sprintf(filename, RESULTSET_CHECKPOINT_TEMPLATE, s->checkpoint_prefix, s->mpirank, i);
    if (resultset2dfm_save(s->results, filename, RESULTSET2DFM_BINARY) < 0) {
      rjmcmc_error("part2d_sample: Failed to write results checkpoint.\n");
      return -1;
    }
    sprintf(filename, MODEL_CHECKPOINT_TEMPLATE, s->checkpoint_prefix, s->mpirank, i);
    if (part2d_forwardmodel_save(s->current, filename) < 0) {
      rjmcmc_error("part2d_sample: Failed to write model checkpoint.\n");
      return -1;
    }
  }

  return 0;
}

