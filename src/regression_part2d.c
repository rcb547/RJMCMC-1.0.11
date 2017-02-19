
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/regression.h>

#include <rjmcmc/engine.h>
#include <rjmcmc/resultset2d.h>
#include <rjmcmc/part2d_regression_rj.h>
#include <rjmcmc/dataset2d.h>
#include <rjmcmc/rjmcmc_util.h>

static double part2d_init(void *arg);
static int part2d_select(void *arg);
static void *part2d_perturb(void *arg, int proc);
static double part2d_misfit(void *arg, void *state);
static int part2d_accept(void *arg, 
			   double current,
			   double proposed);
static int part2d_sample(void *arg,
			   int i);

struct part2d {

  resultset2d_t *results;

  part2d_regression_rj_t *current;
  double current_like;

  part2d_regression_rj_t *proposed;
  double proposed_like;

  int nprocesses;

  int out;
  int accepted;
  int process;

  double birth_prob;
  double death_prob;
  double move_prob;
  double value_prob;
  double order_prob;
  double lambda_prob;
  
  const dataset2d_t *dataset;
  double dpart;
  double dz;
  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;

  int xsamples;
  int ysamples;
  double **v;

  regression2d_cb_t user_callback;
  void *user_arg;
  double *user_x_centre;
  double *user_y_centre;

  int user_ensembles;
  const double *user_ensemble_x;
  const double *user_ensemble_y;

};

static void do_user_callback(struct part2d *s);

resultset2d_t *
part2d_regression(const dataset2d_t *dataset,
		  int burnin,
		  int total,
		  int min_part,
		  int max_part,
		  int xsamples,
		  int ysamples,
		  int zsamples,
		  double credible_interval,
		  double pv,
		  double pd,
		  rjmcmc_uniform_rand_t random,
		  rjmcmc_normal_rand_t normal,
		  int results,
		  regression2d_cb_t user_callback,
		  void *user_arg)
{
  struct part2d s;
  rjmcmc_engine_cb_t cb;

  /*
   * Allocate the results and state information
   */
  if (dataset->lambdastd > 0.0) {
    s.nprocesses = 5;
  } else {
    s.nprocesses = 4;
  }

  s.results = resultset2d_create(burnin,
				 total,
				 xsamples,
				 ysamples,
				 zsamples,
				 s.nprocesses,
				 max_part,
				 dataset->xmin,
				 dataset->xmax,
				 dataset->ymin,
				 dataset->ymax,
				 dataset->zmin,
				 dataset->zmax,
				 credible_interval,
				 results);
  if (s.results == NULL) {
    rjmcmc_error("part2d_regression: failed to create results\n");
    return NULL;
  }

  s.current = part2d_regression_rj_create(min_part,
					  max_part,
					  1,
					  dataset->xmin,
					  dataset->xmax,
					  dataset->ymin,
					  dataset->ymax,
					  pv,
					  pd,
					  pd, 
					  0.0);
  if (s.current == NULL) {
    rjmcmc_error("part2d_regression: failed to create current state\n");
    return NULL;
  }

  s.proposed = part2d_regression_rj_create(min_part,
					   max_part,
					   1,
					   dataset->xmin,
					   dataset->xmax,
					   dataset->ymin,
					   dataset->ymax,
					   pv,
					   pd,
					   pd, 
					   0.0);
  if (s.proposed == NULL) {
    rjmcmc_error("part2d_regression: failed to create proposed state\n");
    return NULL;
  }

  s.dataset = dataset;
  s.dz = dataset->zmax - dataset->zmin;
  s.dpart = max_part - min_part + 1.0;
  if (s.dz <= 0.0) {
    rjmcmc_error(
	    "part2d_regression: the z range of the dataset hasn't been "
	    "set correctly\n");
    return NULL;
  }

  s.random = random;
  s.normal = normal;
  
  s.xsamples = xsamples;
  s.ysamples = ysamples;
  s.v = rjmcmc_create_array_2d(xsamples, ysamples);
  if (s.v == NULL) {
    rjmcmc_error("part2d_regression: failed to create value array\n");
    return NULL;
  }

  s.user_callback = user_callback;
  s.user_arg = user_arg;
  
  if (s.user_callback != NULL) {
    s.user_x_centre = rjmcmc_create_array_1d(max_part);
    s.user_y_centre = rjmcmc_create_array_1d(max_part);
  } else {
    s.user_x_centre = NULL;
    s.user_y_centre = NULL;
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


  resultset2d_assemble_results(s.results);

  rjmcmc_destroy_array_2d(s.xsamples, s.v);

  rjmcmc_destroy_array_1d(s.user_x_centre);
  rjmcmc_destroy_array_1d(s.user_y_centre);

  part2d_regression_rj_destroy(s.current);
  part2d_regression_rj_destroy(s.proposed);

  
  return s.results;
}

static double part2d_init(void *arg)
{
  struct part2d *s = (struct part2d *)arg;

  if (part2d_regression_rj_initialize(s->current,
				      &(s->dataset),
				      1,
				      s->random,
				      s->normal) < 0) {
    rjmcmc_error("part2d_init: failed to initialize\n");
    return -1.0;
  }

  s->current_like = part2d_regression_rj_misfit(s->current,
						&(s->dataset),
						1);
  return s->current_like;
}

static int part2d_select(void *arg)
{
  struct part2d *s = (struct part2d *)arg;
  int i;

  i = rjmcmc_random_choose_int(0, s->nprocesses - 1, s->random);
  return i;
}

static void* part2d_perturb(void *arg, int proc)
{
  struct part2d *s = (struct part2d *)arg;

  s->process = proc;

  resultset2d_propose(s->results, s->process);

  switch (proc) {
  case 0: /* Birth */
    s->out = part2d_regression_rj_propose_birth(s->current,
						s->proposed,
						&(s->dataset),
						1,
						s->random,
						s->normal,
						&(s->birth_prob));
    break;

  case 1: /* Death */
    s->out = part2d_regression_rj_propose_death(s->current,
						s->proposed,
						&(s->dataset),
						1,
						s->random,
						s->normal,
						&(s->death_prob));
    break;

  case 2: /* Move */
    s->out = part2d_regression_rj_propose_move(s->current,
					       s->proposed,
					       &(s->dataset),
					       1,
					       s->random,
					       s->normal,
					       &(s->move_prob));
    break;

  case 3: /* Value */
    s->out = part2d_regression_rj_propose_value(s->current,
						s->proposed,
						&(s->dataset),
						1,
						s->random,
						s->normal,
						&(s->value_prob));
    break;

  case 4: /* Lambda */
    s->out = part2d_regression_rj_propose_lambda(s->current,
						s->proposed,
						&(s->dataset),
						1,
						s->random,
						s->normal,
						&(s->lambda_prob));
    break;

  default:
    return NULL;
  }

  return s->proposed;
}

static double part2d_misfit(void *arg, void *state)
{
  struct part2d *s = (struct part2d *)arg;

  if (s->out == 0) {
    s->proposed_like = FLT_MAX;
  } else {
    s->proposed_like = part2d_regression_rj_misfit(s->proposed,
						   &(s->dataset),
						   1);
  }

  return s->proposed_like;
}

static int part2d_accept(void *arg, 
			 double current_like,
			 double proposed_like)
{
  struct part2d *s = (struct part2d *)arg;
  double u;
  
  if (s->out == 0) {
    return 0;
  }

  u = log(s->random());

  switch (s->process) {
  case 0: /* Birth */
    s->accepted = (u < (current_like - proposed_like -
			log(s->dz) - log(s->birth_prob)));
    break;

  case 1: /* Death */
    s->accepted = (u < (current_like - proposed_like + 
			log(s->dz) + log(s->death_prob)));
    break;

  case 2: /* Move */
    s->accepted = (u < (current_like - proposed_like));
    break;

  case 3: /* Value */
    s->accepted = (u < (current_like - proposed_like));
    break;

  case 4: /* Lambda */
    /* printf("lambda: %g %g %g\n", s->lambda_prob, current_like, proposed_like); */
    s->accepted = (u < (s->lambda_prob + current_like - proposed_like));
    break;
    
  default:
    return 0;
  }

  if (s->accepted) {
    part2d_regression_rj_clone(s->proposed, s->current);
    s->current_like = s->proposed_like;

    resultset2d_accept(s->results, s->process);
  }

  return s->accepted;
}

static int part2d_sample(void *arg,
			 int i)
{
  int j;
  struct part2d *s = (struct part2d *)arg;

  int npart;
  double x;
  double y;

  if (part2d_regression_rj_evaluate(s->current,
				    0,
				    s->dataset->xmin,
				    s->dataset->xmax,
				    s->xsamples,
				    s->dataset->ymin,
				    s->dataset->ymax,
				    s->ysamples,
				    s->v) < 0) {
    rjmcmc_error("part2d_sample: failed to evaluate current state\n");
    return -1;
  }

  resultset2d_sample(s->results,
		     i,
		     (const double **)s->v);

  resultset2d_sample_misfit(s->results,
			    i,
			    s->current_like);

  resultset2d_sample_lambda(s->results,
			   i,
			   part2d_regression_rj_lambda(s->current));

  npart = part2d_regression_rj_partitions(s->current);
  resultset2d_sample_npartitions(s->results,
				 i,
				 npart);

  for (j = 0; j < npart; j ++) {
    if (part2d_regression_rj_partition_centre(s->current, j, &x, &y) >= 0) {
      resultset2d_sample_centre(s->results,
				x, y);
    } else {
      printf("Failed to get centre\n");
    }
  }

  do_user_callback(s);
  
  return 0;
}

static void do_user_callback(struct part2d *s)
{
  int i;
  int n;

  if (s->user_callback != NULL) {
    /*
     * Create the list of cell centres
     */
    n = part2d_regression_rj_partitions(s->current);
    for (i = 0; i < n; i ++) {
      part2d_regression_rj_partition_centre(s->current,
					    i,
					    s->user_x_centre + i,
					    s->user_y_centre + i);
    }

    s->user_callback(s->current,
		     s->user_x_centre,
		     s->user_y_centre,
		     n,
		     (regression2d_value_at_t)part2d_regression_rj_value_at,
		     part2d_regression_rj_lambda(s->current),
		     s->user_arg);
  }
}

