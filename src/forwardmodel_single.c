
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/forwardmodel.h>

#include <rjmcmc/engine.h>
#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_debug.h>
#include <rjmcmc/rjmcmc_random.h>

static double single_init(void *arg);
static int single_select(void *arg);
static void *single_perturb(void *arg, int proc);
static double single_misfit(void *arg, void *state);
static int single_accept(void *arg, 
			   double current,
			   double proposed);
static int single_sample(void *arg,
			   int i);

static double single_hierarchical_init(void *arg);
static int single_hierarchical_select(void *arg);
static void *single_hierarchical_perturb(void *arg, int proc);
static double single_hierarchical_misfit(void *arg, void *state);
static int single_hierarchical_accept(void *arg, 
					double current,
					double proposed);
static int single_hierarchical_sample(void *arg,
					int i);


struct fmsingle {

  resultsetfm_t *results;

  double *current;
  double *current_hierarchical;
  double current_logdetce;
  double current_like;

  double *proposed;
  double *proposed_hierarchical;
  double proposed_logdetce;
  double proposed_like;

  int nparameters;
  const forwardmodelparameter_t *parameters;

  int nhierarchicalparameters;
  const forwardmodelparameter_t *hierarchicalparameters;

  int nsamples;

  int out;
  int process;

  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;

  single_fm_likelihood_t lp;
  single_fm_likelihood_hierarchical_t lp_hierarchical;
  
  void *user_arg;
};

resultsetfm_t *
single_forwardmodel(int burnin,
		    int total,
		    rjmcmc_uniform_rand_t random,
		    rjmcmc_normal_rand_t normal,
		    int nparameters,
		    const forwardmodelparameter_t *parameters,
		    single_fm_likelihood_t lp,
		    void *user_arg,
		    int samples,
		    double credible_interval,
		    int results)
{
  struct fmsingle s;
  rjmcmc_engine_cb_t cb;
  int i;

  s.results = resultsetfm_create(burnin,
				 total,
				 nparameters,
				 parameters,
				 0,
				 NULL,
				 samples,
				 credible_interval,
				 results);
  if (s.results == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = rjmcmc_create_array_1d(nparameters);
  if (s.current == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create model\n");
    return NULL;
  }

  s.current_hierarchical = NULL;
  s.current_logdetce = 0.0;
  s.current_like = 0.0;

  s.proposed = rjmcmc_create_array_1d(nparameters);
  if (s.proposed == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create model\n");
    return NULL;
  }

  s.proposed_hierarchical = NULL;
  s.proposed_logdetce = 0.0;
  s.proposed_like = 0.0;

  s.nparameters = nparameters;
  s.parameters = parameters;

  s.nhierarchicalparameters = 0;
  s.hierarchicalparameters = NULL;

  s.nsamples = samples;

  s.out = 0;
  s.process = 0;

  s.random = random;
  s.normal = normal;

  s.lp = lp;
  s.lp_hierarchical = NULL;

  s.user_arg = user_arg;
  
  cb.initialize_state = single_init;
  cb.select_process = single_select;
  cb.perturb_process = single_perturb;
  cb.compute_misfit = single_misfit;
  cb.accept = single_accept;
  cb.sample = single_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.current);
  rjmcmc_destroy_array_1d(s.proposed);

  resultsetfm_assemble_results(s.results);

  return s.results;
}

resultsetfm_t *
single_forwardmodel_hierarchical(int burnin,
				 int total,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 int nparameters,
				 const forwardmodelparameter_t *parameters,
				 int nhierarchicalparameters,
				 const forwardmodelparameter_t *hierarchicalparameters,
				 single_fm_likelihood_hierarchical_t lp,
				 void *user_arg,
				 int samples,
				 double credible_interval,
				 int results)
{
  struct fmsingle s;
  rjmcmc_engine_cb_t cb;

  if (nhierarchicalparameters <= 0) {
    rjmcmc_error("single_forwardmodel_hierarchical: "
		 "invalid number of hierarchical parameters\n");
    return NULL;
  }

  s.results = resultsetfm_create(burnin,
				 total,
				 nparameters,
				 parameters,
				 nhierarchicalparameters,
				 hierarchicalparameters,
				 samples,
				 credible_interval,
				 results);

  if (s.results == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create results\n");
    return NULL;
  }

  s.current = rjmcmc_create_array_1d(nparameters);
  if (s.current == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create model\n");
    return NULL;
  }

  s.current_hierarchical = rjmcmc_create_array_1d(nhierarchicalparameters);
  if (s.current_hierarchical == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create model\n");
    return NULL;
  }

  s.current_logdetce = 0.0;
  s.current_like = 0.0;

  s.proposed = rjmcmc_create_array_1d(nparameters);
  if (s.proposed == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create model\n");
    return NULL;
  }

  s.proposed_hierarchical = rjmcmc_create_array_1d(nhierarchicalparameters);
  if (s.proposed_hierarchical == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create model\n");
    return NULL;
  }

  s.proposed_logdetce = 0.0;
  s.proposed_like = 0.0;

  s.nparameters = nparameters;
  s.parameters = parameters;

  s.nhierarchicalparameters = nhierarchicalparameters;
  s.hierarchicalparameters = hierarchicalparameters;

  s.nsamples = samples;

  s.out = 0;
  s.process = 0;

  s.random = random;
  s.normal = normal;

  s.lp = NULL;
  s.lp_hierarchical = lp;
  s.user_arg = user_arg;
  
  cb.initialize_state = single_hierarchical_init;
  cb.select_process = single_hierarchical_select;
  cb.perturb_process = single_hierarchical_perturb;
  cb.compute_misfit = single_hierarchical_misfit;
  cb.accept = single_hierarchical_accept;
  cb.sample = single_hierarchical_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.current);
  rjmcmc_destroy_array_1d(s.current_hierarchical);

  rjmcmc_destroy_array_1d(s.proposed);
  rjmcmc_destroy_array_1d(s.proposed_hierarchical);

  resultsetfm_assemble_results(s.results);

  return s.results;
}

/*
 * Basic Forward Model Callbacks
 */

static double single_init(void *arg)
{
  int i;

  struct fmsingle *s = (struct fmsingle *)arg;

  for (i = 0; i < s->nparameters; i ++) {

    s->current[i] = s->parameters[i].fmin + 
      (s->random)() * (s->parameters[i].fmax - s->parameters[i].fmin);

  }

  s->current_like = (s->lp)(s->user_arg, s->nparameters, s->current);
  return s->current_like;
}

static int single_select(void *arg)
{
  struct fmsingle *s = (struct fmsingle *)arg;
  int p;

  if (s->nparameters == 1) {
    p = 0; 
  } else {
    p = rjmcmc_random_choose_int(0, s->nparameters - 1, s->random);
  }

  return p;
}

static void *single_perturb(void *arg, int proc)
{
  struct fmsingle *s = (struct fmsingle *)arg;

  s->process = proc;

  memcpy(s->proposed, s->current, sizeof(double) * s->nparameters);

  s->proposed[proc] += s->normal() * s->parameters[proc].fstd_value;
  s->out = 
    (s->proposed[proc] >= s->parameters[proc].fmin) && 
    (s->proposed[proc] <= s->parameters[proc].fmax);

  resultsetfm_propose(s->results, proc);

  return s->proposed;
}

static double single_misfit(void *arg, void *state)
{
  struct fmsingle *s = (struct fmsingle *)arg;

  if (s->out) {
    s->proposed_like = (s->lp)(s->user_arg, s->nparameters, s->proposed);

    return s->proposed_like;
  } else {
    return FLT_MAX;
  }
}

static int single_accept(void *arg, 
			   double current,
			   double proposed)
{
  struct fmsingle *s = (struct fmsingle *)arg;
  double u;
  int accepted;

  if (s->out == 0) {
    return 0;
  }

  u = log(s->random());
  accepted = (u < (current - proposed));

  if (accepted) {
    memcpy(s->current, s->proposed, sizeof(double) * s->nparameters);
    s->current_like = s->proposed_like;

    resultsetfm_accept(s->results, s->process);
  } 

  return accepted;
}

static int single_sample(void *arg,
			   int i)
{
  struct fmsingle *s = (struct fmsingle *)arg;

  resultsetfm_sample(s->results,
		     i,
		     s->current);
  resultsetfm_sample_misfit(s->results,
			    i,
			    s->current_like);

  return 0;
}

/*
 * Hierarchical Callbacks 
 */

static double single_hierarchical_init(void *arg)
{
  int i;

  struct fmsingle *s = (struct fmsingle *)arg;

  for (i = 0; i < s->nparameters; i ++) {
    s->current[i] = 
      rjmcmc_random_choose_double(s->parameters[i].fmin,
				  s->parameters[i].fmax,
				  s->random);
  }

  for (i = 0; i < s->nhierarchicalparameters; i ++) {
    s->current_hierarchical[i] = 
      rjmcmc_random_choose_double(s->hierarchicalparameters[i].fmin,
				  s->hierarchicalparameters[i].fmax,
				  s->random);
  }

  s->current_like = (s->lp_hierarchical)(s->user_arg, 
					 s->nparameters, 
					 s->current,
					 1, /* Require calculation of log(detCe) */
					 s->nhierarchicalparameters,
					 s->current_hierarchical,
					 &(s->current_logdetce));

  return s->current_like;
}

static int single_hierarchical_select(void *arg)
{
  struct fmsingle *s = (struct fmsingle *)arg;

  return rjmcmc_random_choose_int(0, 
				  s->nparameters + 
				  s->nhierarchicalparameters - 1, 
				  s->random);
}

static void *single_hierarchical_perturb(void *arg, int proc)
{
  struct fmsingle *s = (struct fmsingle *)arg;
  int hproc;

  s->process = proc;
  resultsetfm_propose(s->results, s->process);

  memcpy(s->proposed, s->current, sizeof(double) * s->nparameters);
  memcpy(s->proposed_hierarchical, s->current_hierarchical, sizeof(double) * s->nhierarchicalparameters);

  if (proc < s->nparameters) {

    /*
     * Normal parameter peturbation
     */
    s->proposed[proc] += s->normal() * s->parameters[proc].fstd_value;
    s->out = 
      (s->proposed[proc] >= s->parameters[proc].fmin) && 
      (s->proposed[proc] <= s->parameters[proc].fmax);

  } else {

    /*
     * Hierarchical parameter peturbation
     */
    hproc = proc - s->nparameters;

    s->proposed_hierarchical[hproc] += 
      s->normal() * s->hierarchicalparameters[hproc].fstd_value;

    s->out = 
      (s->proposed_hierarchical[hproc] >= s->hierarchicalparameters[hproc].fmin) &&
      (s->proposed_hierarchical[hproc] <= s->hierarchicalparameters[hproc].fmax);


  }

  return s->proposed;
}

static double single_hierarchical_misfit(void *arg, void *state)
{
  struct fmsingle *s = (struct fmsingle *)arg;
  int hierarchical;

  if (s->out) {
    hierarchical = (s->process >= s->nparameters);
    
    s->proposed_like = (s->lp_hierarchical)(s->user_arg, 
					    s->nparameters, 
					    s->proposed, 
					    hierarchical,
					    s->nhierarchicalparameters,
					    s->proposed_hierarchical,
					    &(s->proposed_logdetce));
    return s->proposed_like;

  } else {

    return FLT_MAX;

  }
}

static int single_hierarchical_accept(void *arg, 
				      double current,
				      double proposed)
{
  struct fmsingle *s = (struct fmsingle *)arg;
  double u;
  int accepted;

  if (s->out == 0) {
    return 0;
  }

  u = log(s->random());

  if (s->process < s->nparameters) {
    /* Normal peturbation */
    accepted = (u < (current - proposed));

  } else {
    /* Hierarchical peturbation */
    /* printf("    %f %f %f\n", */
    /* 	   s->current_logdetce, */
    /* 	   s->proposed_logdetce, */
    /* 	   s->current_logdetce - s->proposed_logdetce); */
    accepted = 
      (u < 
       (0.5 * (s->current_logdetce - s->proposed_logdetce) +  
	current - proposed));
  }

  if (accepted) {
    memcpy(s->current, s->proposed, sizeof(double) * s->nparameters);
    memcpy(s->current_hierarchical, s->proposed_hierarchical, sizeof(double) * s->nhierarchicalparameters);

    s->current_like = s->proposed_like;
    s->current_logdetce = s->proposed_logdetce;

    resultsetfm_accept(s->results, s->process);
  } 

  return accepted;
}

static int single_hierarchical_sample(void *arg,
					int i)
{
  struct fmsingle *s = (struct fmsingle *)arg;

  resultsetfm_sample_misfit(s->results,
			    i,
			    s->current_like);

  resultsetfm_sample(s->results,
		     i,
		     s->current);

  resultsetfm_sample_hierarchical(s->results,
				  i,
				  s->current_hierarchical);

  return 0;
}

#if defined(HAVE_MPI_H)

resultsetfm_t *
MPI_single_forwardmodel(int burnin,
			int total,
			rjmcmc_uniform_rand_t random,
			rjmcmc_normal_rand_t normal,
			int nparameters,
			const forwardmodelparameter_t *parameters,
			single_fm_likelihood_t lp,
			void *user_arg,
			int samples,
			double credible_interval,
			int results,
			int mpisize,
			int mpirank,
			int root,
			MPI_Comm comm)
{
  struct fmsingle s;
  rjmcmc_engine_cb_t cb;
  int i;

  s.results = resultsetfm_create(burnin,
				 total,
				 nparameters,
				 parameters,
				 0,
				 NULL,
				 samples,
				 credible_interval,
				 results);
  if (s.results == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create results\n");
    return NULL;
  }

  s.current = rjmcmc_create_array_1d(nparameters);
  if (s.current == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create model\n");
    return NULL;
  }

  s.current_hierarchical = NULL;
  s.current_logdetce = 0.0;
  s.current_like = 0.0;

  s.proposed = rjmcmc_create_array_1d(nparameters);
  if (s.proposed == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create model\n");
    return NULL;
  }

  s.proposed_hierarchical = NULL;
  s.proposed_logdetce = 0.0;
  s.proposed_like = 0.0;

  s.nparameters = nparameters;
  s.parameters = parameters;

  s.nhierarchicalparameters = 0;
  s.hierarchicalparameters = NULL;

  s.nsamples = samples;

  s.out = 0;
  s.process = 0;

  s.random = random;
  s.normal = normal;

  s.lp = lp;
  s.lp_hierarchical = NULL;

  s.user_arg = user_arg;
  
  cb.initialize_state = single_init;
  cb.select_process = single_select;
  cb.perturb_process = single_perturb;
  cb.compute_misfit = single_misfit;
  cb.accept = single_accept;
  cb.sample = single_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.current);
  rjmcmc_destroy_array_1d(s.proposed);

  MPI_resultsetfm_assemble_results(s.results,
				   mpisize,
				   mpirank,
				   root,
				   comm);

  return s.results;
}

resultsetfm_t *
MPI_single_forwardmodel_hierarchical(int burnin,
				     int total,
				     rjmcmc_uniform_rand_t random,
				     rjmcmc_normal_rand_t normal,
				     int nparameters,
				     const forwardmodelparameter_t *parameters,
				     int nhierarchicalparameters,
				     const forwardmodelparameter_t *hierarchicalparameters,
				     single_fm_likelihood_hierarchical_t lp,
				     void *user_arg,
				     int samples,
				     double credible_interval,
				     int results,
				     int mpisize,
				     int mpirank,
				     int root,
				     MPI_Comm comm)
{
  struct fmsingle s;
  rjmcmc_engine_cb_t cb;

  if (nhierarchicalparameters <= 0) {
    rjmcmc_error("single_forwardmodel_hierarchical: "
		 "invalid number of hierarchical parameters\n");
    return NULL;
  }

  s.results = resultsetfm_create(burnin,
				 total,
				 nparameters,
				 parameters,
				 nhierarchicalparameters,
				 hierarchicalparameters,
				 samples,
				 credible_interval,
				 results);

  if (s.results == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create results\n");
    return NULL;
  }

  s.current = rjmcmc_create_array_1d(nparameters);
  if (s.current == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create model\n");
    return NULL;
  }

  s.current_hierarchical = rjmcmc_create_array_1d(nhierarchicalparameters);
  if (s.current_hierarchical == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create model\n");
    return NULL;
  }

  s.current_logdetce = 0.0;
  s.current_like = 0.0;

  s.proposed = rjmcmc_create_array_1d(nparameters);
  if (s.proposed == NULL) {
    rjmcmc_error("single_forwardmodel: failed to create model\n");
    return NULL;
  }

  s.proposed_hierarchical = rjmcmc_create_array_1d(nhierarchicalparameters);
  if (s.proposed_hierarchical == NULL) {
    rjmcmc_error(
	    "single_forwardmodel_hierarchical: failed to create model\n");
    return NULL;
  }

  s.nparameters = nparameters;
  s.parameters = parameters;

  s.nhierarchicalparameters = nhierarchicalparameters;
  s.hierarchicalparameters = hierarchicalparameters;

  s.nsamples = samples;

  s.out = 0;
  s.process = 0;

  s.random = random;
  s.normal = normal;

  s.lp = NULL;
  s.lp_hierarchical = lp;
  s.user_arg = user_arg;
  
  cb.initialize_state = single_hierarchical_init;
  cb.select_process = single_hierarchical_select;
  cb.perturb_process = single_hierarchical_perturb;
  cb.compute_misfit = single_hierarchical_misfit;
  cb.accept = single_hierarchical_accept;
  cb.sample = single_hierarchical_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }  

  rjmcmc_destroy_array_1d(s.current);
  rjmcmc_destroy_array_1d(s.current_hierarchical);

  rjmcmc_destroy_array_1d(s.proposed);
  rjmcmc_destroy_array_1d(s.proposed_hierarchical);

  MPI_resultsetfm_assemble_results(s.results,
				   mpisize,
				   mpirank,
				   root,
				   comm);

  return s.results;
}

#endif /* HAVE_MPI_H */
