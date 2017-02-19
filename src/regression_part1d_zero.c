
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/regression.h>

#include <rjmcmc/engine.h>
#include <rjmcmc/resultset1d.h>
#include <rjmcmc/part1d_zero.h>
#include <rjmcmc/dataset1d.h>
#include <rjmcmc/rjmcmc_util.h>

static double part1d_init(void *arg);
static int part1d_select(void *arg);
static void *part1d_perturb(void *arg, int proc);
static double part1d_misfit(void *arg, void *state);
static int part1d_accept(void *arg, 
			 double current,
			 double proposed);
static int part1d_sample(void *arg,
			 int i);

struct part1d {

  resultset1d_t *results;

  part1d_zero_t *current;
  double current_like;

  part1d_zero_t *proposed;
  double proposed_like;

  int nprocesses;

  int out;
  int process;

  double birth_prob;
  double death_prob;
  double move_prob;
  double value_prob;
  double lambda_prob;

  const dataset1d_t *dataset;
  double dy;
  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;

  int xsamples;
  double *v;

  regression1d_cb_t user_callback;
  void *user_arg;
  double *partitions;

};

static int doublecmp(const void *_a, const void *_b)
{
  const double *a = (const double*)_a;
  const double *b = (const double*)_b;

  if (*a == *b) {
    return 0;
  } else if (*a < *b) {
    return -1;
  } else {
    return 1;
  }
}

static void
do_user_callback(struct part1d *p)
{
  int i;
  int n;
  
  if (p->user_callback != NULL) {
    n = part1d_zero_partitions(p->current);
    for (i = 0; i < n; i ++) {
      p->partitions[i] = part1d_zero_partition_position(p->current, i);
    }

    qsort(p->partitions, n, sizeof(double), doublecmp);
    
    (p->user_callback)(p->current,
		       p->partitions,
		       part1d_zero_partitions(p->current),
		       (regression1d_value_at_t)part1d_zero_value_at,
		       part1d_zero_lambda(p->current, 0),
		       p->user_arg);
  }
}

resultset1d_t *
part1d_zero_regression(const dataset1d_t *dataset,
		       int burnin,
		       int total,
		       int min_part,
		       int max_part,
		       int xsamples,
		       int ysamples,
		       double credible_interval,
		       double pd,
		       rjmcmc_uniform_rand_t random,
		       rjmcmc_normal_rand_t normal,
		       int results,
		       regression1d_cb_t callback,
		       void *user_arg)
{
  struct part1d s;
  rjmcmc_engine_cb_t cb;

  /*
   * Allocate the results and state information
   */
  if (dataset->lambdastd == 0.0) {
    s.nprocesses = 4;
  } else {
    s.nprocesses = 5;
  }

  s.results = resultset1d_create(burnin,
				 total,
				 xsamples,
				 ysamples,
				 s.nprocesses,
				 max_part,
				 0,
				 dataset->xmin,
				 dataset->xmax,
				 dataset->ymin,
				 dataset->ymax,
				 credible_interval,
				 results);
  if (s.results == NULL) {
    rjmcmc_error("part1d_regression: failed to create results\n");
    return NULL;
  }

  s.current = part1d_zero_create(min_part,
				 max_part,
				 1,
				 dataset->xmin,
				 dataset->xmax,
				 pd);
  if (s.current == NULL) {
    rjmcmc_error("part1d_regression: failed to create current state\n");
    return NULL;
  }
  
  s.proposed = part1d_zero_create(min_part,
				  max_part,
				  1,
				  dataset->xmin,
				  dataset->xmax,
				  pd);
  if (s.proposed == NULL) {
    rjmcmc_error("part1d_regression: failed to create proposed state\n");
    return NULL;
  }

  s.dataset = dataset;
  s.dy = dataset->ymax - dataset->ymin;
  s.random = random;
  s.normal = normal;
  
  s.xsamples = xsamples;
  s.v = rjmcmc_create_array_1d(xsamples);
  if (s.v == NULL) {
    rjmcmc_error("part1d_regression: failed to create value array\n");
    return NULL;
  }

  s.user_callback = callback;
  s.user_arg = user_arg;
  s.partitions = rjmcmc_create_array_1d(max_part);
  if (s.partitions == NULL) {
    rjmcmc_error(
	    "regression_part1d_natural: failed to create partitions array\n");
    return NULL;
  }

  /*
   * Set the engine callbacks 
   */
  cb.initialize_state = part1d_init;
  cb.select_process = part1d_select;
  cb.perturb_process = part1d_perturb;
  cb.compute_misfit = part1d_misfit;
  cb.accept = part1d_accept;
  cb.sample = part1d_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }

  rjmcmc_destroy_array_1d(s.v);
  rjmcmc_destroy_array_1d(s.partitions);

  part1d_zero_destroy(s.current);
  part1d_zero_destroy(s.proposed);

  resultset1d_assemble_results(s.results);

  return s.results;
}

static double part1d_init(void *arg)
{
  struct part1d *s = (struct part1d *)arg;

  if (part1d_zero_initialize(s->current,
				      &(s->dataset),
				      1,
				      s->random,
				      s->normal) < 0) {
    rjmcmc_error("part1d_init: failed to initialize\n");
    return -1.0;
  }

  s->current_like = part1d_zero_misfit(s->current,
				       &(s->dataset),
				       1);
  return s->current_like;
}

static int part1d_select(void *arg)
{
  struct part1d *s = (struct part1d *)arg;

  return (int)(s->random() * (double)s->nprocesses);
}

static void* part1d_perturb(void *arg, int proc)
{
  struct part1d *s = (struct part1d *)arg;

  s->process = proc;

  resultset1d_propose(s->results, s->process);

  switch (proc) {
  case 0: /* Birth */
    s->out = part1d_zero_propose_birth(s->current,
				       s->proposed,
				       &(s->dataset),
				       1,
				       s->random,
				       s->normal,
				       &(s->birth_prob));
    break;

  case 1: /* Death */
    s->out = part1d_zero_propose_death(s->current,
					     s->proposed,
					     &(s->dataset),
					     1,
					     s->random,
					     s->normal,
					     &(s->death_prob));
    break;

  case 2: /* Move */
    s->out = part1d_zero_propose_move(s->current,
				      s->proposed,
				      &(s->dataset),
				      1,
				      s->random,
				      s->normal,
				      &(s->move_prob));
    break;

  case 3: /* Value */
    s->out = part1d_zero_propose_value(s->current,
				       s->proposed,
				       &(s->dataset),
				       1,
				       s->random,
				       s->normal,
				       &(s->value_prob));
    break;

  case 4: /* Lambda */
    s->out = part1d_zero_propose_lambda(s->current,
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

static double part1d_misfit(void *arg, void *state)
{
  struct part1d *s = (struct part1d *)arg;

  if (s->out == 0) {
    s->proposed_like = FLT_MAX;
  } else {
    s->proposed_like = part1d_zero_misfit(s->proposed,
					  &(s->dataset),
					  1);
  }
  return s->proposed_like;
}

static int part1d_accept(void *arg, 
			 double current_like,
			 double proposed_like)
{
  struct part1d *s = (struct part1d *)arg;
  double u;
  int accepted;
  double dj;

  accepted = 0;
  if (s->out == 0) {
    return 0;
  }

  u = log(s->random());

  switch (s->process) {
  case 0: /* Birth */
    if (u < (current_like - proposed_like - log(s->dy) - log(s->birth_prob))) {
      accepted = 1;
    }
    break;

  case 1: /* Death */
    if (u < (current_like - proposed_like + log(s->dy) + log(s->death_prob))) {
      accepted = 1;
    }
    break;

  case 2: /* Move */
    if (u < (current_like - proposed_like)) {
      accepted = 1;
    }
    break;

  case 3: /* Value */
    if (u < (current_like - proposed_like)) {
      accepted = 1;
    }
    break;

  case 4: /* Lambda */
    if (u < (current_like - proposed_like + log(s->lambda_prob))) {
      accepted = 1;
    }
    break;
    
  default:
    return 0;
  }
  
  if (accepted) {
    part1d_zero_clone(s->proposed, s->current);
    s->current_like = s->proposed_like;
    
    resultset1d_accept(s->results, s->process);
  }

  return accepted;
}

static int part1d_sample(void *arg,
			   int i)
{
  struct part1d *s = (struct part1d *)arg;
  int j;
  int npart;

  if (part1d_zero_evaluate(s->current,
			   0,
			   s->dataset->xmin,
			   s->dataset->xmax,
			   s->xsamples,
			   s->v) < 0) {
    rjmcmc_error("part1d_sample: failed to evaluate current state\n");
    return -1;
  }

  resultset1d_sample(s->results,
		     i,
		     s->v);

  resultset1d_sample_misfit(s->results,
			    i,
			    s->current_like);

  resultset1d_sample_lambda(s->results,
			   i,
			   part1d_zero_lambda(s->current, 0));

  npart = part1d_zero_partitions(s->current);
  resultset1d_sample_npartitions(s->results,
				 i,
				 npart);
  
  for (j = 0; j < npart; j ++) {
    resultset1d_sample_partition_x(s->results,
				   i,
				   part1d_zero_partition_position(s->current,
								  j));
  }

  do_user_callback(s);

  return 0;
}

