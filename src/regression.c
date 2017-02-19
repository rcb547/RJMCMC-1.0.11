
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <rjmcmc/regression.h>

#include <rjmcmc/engine.h>
#include <rjmcmc/resultset1d.h>
#include <rjmcmc/single1d_regression.h>
#include <rjmcmc/dataset1d.h>
#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/curvefit.h>

static double single1d_init(void *arg);
static int single1d_select(void *arg);
static void *single1d_perturb(void *arg, int proc);
static double single1d_misfit(void *arg, void *state);
static int single1d_accept(void *arg, 
			   double current,
			   double proposed);
static int single1d_sample(void *arg,
			   int i);

struct single1d {

  resultset1d_t *results;

  single1d_regression_t *current;
  double current_like;

  single1d_regression_t *proposed;
  double proposed_like;

  int nprocesses;
  int max_order;

  int out;
  int process;

  double value_prob;
  double order_prob;
  double lambda_prob;
  
  const dataset1d_t *dataset;
  rjmcmc_uniform_rand_t random;
  rjmcmc_normal_rand_t normal;

  int xsamples;
  double *v;

  double *mean_misfit;
  double *detCm;
  double *autoprior;
  double **S;
  double *pk;
  double *kcdf;
  double *prodprior;
  double **mean;
  double **sigma;

  regression1d_cb_t user_callback;
  void *user_arg;

  double *boundaries;

};

typedef void (*resultassembler_t)(void *user, resultset1d_t *r);

static void
do_user_callback(struct single1d *s)
{
  if (s->user_callback != NULL) {
    (s->user_callback)(s->current,
		       NULL,
		       0,
		       (regression1d_value_at_t)single1d_regression_value_at,
		       single1d_regression_lambda(s->current),
		       s->user_arg);
  }
}

static resultset1d_t *
do_single1d_regression(const dataset1d_t *dataset,
		       const double *prior,
		       int useautoprior,
		       int burnin,
		       int total,
		       int max_order,
		       int xsamples,
		       int ysamples,
		       double credible_interval,
		       rjmcmc_uniform_rand_t random,
		       rjmcmc_normal_rand_t normal,
		       int results,
		       regression1d_cb_t user_callback,
		       void *user_arg,
		       resultassembler_t assemble,
		       void *assemble_arg)
{
  struct single1d s;
  rjmcmc_engine_cb_t cb;

  /*
   * Allocate the results and state information
   */
  if (dataset->lambdastd == 0.0) {
    s.nprocesses = 1;
  } else {
    s.nprocesses = 2;
  }
  s.max_order = max_order;

  s.results = resultset1d_create(burnin,
				 total,
				 xsamples,
				 ysamples,
				 s.nprocesses,
				 0,
				 s.max_order,
				 dataset->xmin,
				 dataset->xmax,
				 dataset->ymin,
				 dataset->ymax,
				 credible_interval,
				 results);
  if (s.results == NULL) {
    rjmcmc_error("single1d_regression: failed to create results\n");
    return NULL;
  }


  s.dataset = dataset;
  s.random = random;
  s.normal = normal;
  
  s.xsamples = xsamples;
  s.v = rjmcmc_create_array_1d(xsamples);
  if (s.v == NULL) {
    rjmcmc_error("single1d_regression: failed to create value array\n");
    return NULL;
  }

  s.mean = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (s.mean == NULL) {
    return NULL;
  }
  s.sigma = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (s.sigma == NULL) {
    return NULL;
  }

  s.mean_misfit = rjmcmc_create_array_1d(max_order + 1);
  if (s.mean_misfit == NULL) {
    return NULL;
  }

  s.detCm = rjmcmc_create_array_1d(max_order + 1);
  if (s.detCm == NULL) {
    return NULL;
  }
  
  s.autoprior = rjmcmc_create_array_1d(max_order + 1);
  if (s.autoprior == NULL) {
    return NULL;
  }

  s.S = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (s.S == NULL) {
    return NULL;
  }

  s.pk = rjmcmc_create_array_1d(max_order + 1);
  if (s.pk == NULL) {
    return NULL;
  }

  s.kcdf = rjmcmc_create_array_1d(max_order + 1);
  if (s.kcdf == NULL) {
    return NULL;
  }

  s.current = single1d_regression_create(max_order,
					 dataset->xmin,
					 dataset->xmax, 
					 s.kcdf);
  if (s.current == NULL) {
    rjmcmc_error("single1d_regression: failed to create current state\n");
    return NULL;
  }

  s.proposed = single1d_regression_create(max_order,
					  dataset->xmin,
					  dataset->xmax,
					  s.kcdf);
  if (s.proposed == NULL) {
    rjmcmc_error("single1d_regression: failed to create proposed state\n");
    return NULL;
  }

  if (single1d_evaluate_pk(s.current,
			   s.dataset,
			   prior,
			   3.0,
			   s.mean_misfit,
			   s.detCm,
			   s.autoprior,
			   s.S,
			   s.pk,
			   s.kcdf,
			   s.mean,
			   s.sigma) < 0) {
    rjmcmc_error(
	    "single1d_regression: failed to compute kcdf\n");
    return NULL;
  }

  s.user_callback = user_callback;
  s.user_arg = user_arg;
    
  /*
   * Set the engine callbacks 
   */
  cb.initialize_state = single1d_init;
  cb.select_process = single1d_select;
  cb.perturb_process = single1d_perturb;
  cb.compute_misfit = single1d_misfit;
  cb.accept = single1d_accept;
  cb.sample = single1d_sample;
  cb.arg = (void*)&s;

  if (rjmcmc_engine_run(&cb,
			burnin,
			total,
			1) < 0) {
    return NULL;
  }

  rjmcmc_destroy_array_1d(s.v);

  rjmcmc_destroy_array_2d(s.max_order + 1, s.mean);
  rjmcmc_destroy_array_2d(s.max_order + 1, s.sigma);
  rjmcmc_destroy_array_1d(s.mean_misfit);
  rjmcmc_destroy_array_1d(s.detCm);
  rjmcmc_destroy_array_1d(s.autoprior);
  rjmcmc_destroy_array_2d(s.max_order + 1, s.S);
  rjmcmc_destroy_array_1d(s.pk);
  rjmcmc_destroy_array_1d(s.kcdf);

  single1d_regression_destroy(s.current);
  single1d_regression_destroy(s.proposed);

  assemble(assemble_arg, s.results);

  return s.results;
}

static void normal_assemble(void *user, resultset1d_t *r)
{
  resultset1d_assemble_results(r);
}

resultset1d_t *
single1d_regression(const dataset1d_t *dataset,
		    int burnin,
		    int total,
		    int max_order,
		    int xsamples,
		    int ysamples,
		    double credible_interval,
		    rjmcmc_uniform_rand_t random,
		    rjmcmc_normal_rand_t normal,
		    int results,
		    regression1d_cb_t user_callback,
		    void *user_arg)
{
  return do_single1d_regression(dataset, 
				NULL,
				-1,
				burnin, 
				total, 
				max_order, 
				xsamples,
				ysamples,
				credible_interval,
				random,
				normal,
				results,
				user_callback,
				user_arg,
				normal_assemble,
				NULL);
}

resultset1d_t *
single1d_regression_with_prior(const dataset1d_t *dataset,
			       const double *prior,
			       int burnin,
			       int total,
			       int max_order,
			       int xsamples,
			       int ysamples,
			       double credible_interval,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal,
			       int results,
			       regression1d_cb_t user_callback,
			       void *user_arg)
{
  return do_single1d_regression(dataset, 
				prior,
				0,
				burnin, 
				total, 
				max_order, 
				xsamples,
				ysamples,
				credible_interval,
				random,
				normal,
				results,
				user_callback,
				user_arg,
				normal_assemble,
				NULL);
}

resultset1d_t *
single1d_direct_regression(const dataset1d_t *dataset,
			   const double *fixed_prior,
			   int max_order,
			   int xsamples,
			   rjmcmc_uniform_rand_t random,
			   rjmcmc_normal_rand_t normal)
{
  curvefit_result_t *cf;
  resultset1d_t *results;

  double *mean_misfit;
  double *detCm;
  double *prior;
  double **S;
  double *pk;
  double *kcdf;
  double **mean;
  double **sigma;
  double *weighted_curve;

  int j;
  int k;
  double x;
  double dx;

  cf = curvefit_create(max_order);
  if (cf == NULL) {
    return NULL;
  }

  results = resultset1d_create(0, 
			       0, 
			       xsamples, 
			       0,
			       0, 
			       0, 
			       max_order,
			       dataset->xmin,
			       dataset->xmax,
			       dataset->ymin,
			       dataset->ymax,
			       0.0,
			       RESULTSET1D_MEAN);
  if (results == NULL) {
    return NULL;
  }

  mean_misfit = rjmcmc_create_array_1d(max_order + 1);
  if (mean_misfit == NULL) {
    return NULL;
  }

  detCm = rjmcmc_create_array_1d(max_order + 1);
  if (detCm == NULL) {
    return NULL;
  }

  prior = rjmcmc_create_array_1d(max_order + 1);
  if (prior == NULL) {
    return NULL;
  }

  S = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (S == NULL) {
    return NULL;
  }

  pk = rjmcmc_create_array_1d(max_order + 1);
  if (pk == NULL) {
    return NULL;
  }

  kcdf = rjmcmc_create_array_1d(max_order + 1);
  if (kcdf == NULL) {
    return NULL;
  }

  mean = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (mean == NULL) {
    return NULL;
  }

  sigma = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  if (sigma == NULL) {
    return NULL;
  }

  weighted_curve = rjmcmc_create_array_1d(xsamples);
  if (weighted_curve == NULL) {
    return NULL;
  }
    
  if (curvefit_evaluate_pk(cf,
			   dataset,
			   0,
			   dataset->npoints - 1,
			   max_order,
			   fixed_prior,
			   3.0,
			   mean_misfit,
			   detCm,
			   prior,
			   S,
			   pk,
			   kcdf,
			   mean,
			   sigma) < 0) {
    rjmcmc_error("single1d_direct_regression: failed to evaluate pk\n");
    return NULL;
  }

  dx = (dataset->xmax - dataset->xmin)/(double)(xsamples - 1);
  for (k = 0; k <= max_order; k ++) {
    for (j = 0, x = dataset->xmin; j < xsamples; j ++, x += dx) {
      weighted_curve[j] += 
	rjmcmc_polynomial_value(mean[k], k + 1, x) * pk[k];
    }
  }

  resultset1d_sample(results, 0, weighted_curve);

  rjmcmc_destroy_array_1d(mean_misfit);
  rjmcmc_destroy_array_1d(detCm);
  rjmcmc_destroy_array_1d(prior);
  rjmcmc_destroy_array_2d(max_order + 1, S);
  rjmcmc_destroy_array_1d(pk);
  rjmcmc_destroy_array_1d(kcdf);
  rjmcmc_destroy_array_2d(max_order + 1, mean);
  rjmcmc_destroy_array_2d(max_order + 1, sigma);
  rjmcmc_destroy_array_1d(weighted_curve);

  return results;
}

#if defined(HAVE_MPI_H)

struct mpi_data {
  int size;
  int rank;
  int root;
  MPI_Comm comm;
};

static void MPI_assemble(void *user, resultset1d_t *r)
{
  struct mpi_data *d = (struct mpi_data *)user;
  MPI_resultset1d_assemble_results(r, 
				   d->size,
				   d->rank,
				   d->root,
				   d->comm);
}

resultset1d_t *
MPI_single1d_regression(const dataset1d_t *dataset,
			int burnin,
			int total,
			int max_order,
			int xsamples,
			int ysamples,
			double credible_interval,
			rjmcmc_uniform_rand_t random,
			rjmcmc_normal_rand_t normal,
			int results,
			regression1d_cb_t user_callback,
			void *user_arg,
			int mpisize,
			int mpirank,
			int root,
			MPI_Comm comm)
{
  struct mpi_data d;

  d.size = mpisize;
  d.rank = mpirank;
  d.root = root;
  d.comm = comm;

  return do_single1d_regression(dataset, 
				NULL,
				0,
				burnin, 
				total, 
				max_order, 
				xsamples,
				ysamples,
				credible_interval,
				random,
				normal,
				results,
				user_callback,
				user_arg,
				MPI_assemble,
				&d);
}

#endif /* HAVE_OPENMPI_MPI_H */

static double single1d_init(void *arg)
{
  struct single1d *s = (struct single1d *)arg;

  if (single1d_regression_initialize(s->current,
				     s->dataset,
				     s->random,
				     s->normal) < 0) {
    rjmcmc_error("single1d_init: failed to initialize\n");
    return -1.0;
  }

  s->current_like = single1d_regression_misfit(s->current,
					       s->dataset);

  s->process = -1;

  return s->current_like;
}

static int single1d_select(void *arg)
{
  struct single1d *s = (struct single1d *)arg;

  if (s->nprocesses == 1) {
    return 0;
  } else {
    return rjmcmc_random_choose_int(0, s->nprocesses - 1, s->random);
  }
}

static void* single1d_perturb(void *arg, int proc)
{
  struct single1d *s = (struct single1d *)arg;

  s->process = proc;

  resultset1d_propose(s->results, s->process);


  switch (proc) {
  case 0: /* Value */
    s->out = single1d_regression_propose_value(s->current,
					       s->proposed,
					       s->dataset,
					       s->random,
					       s->normal,
					       &(s->value_prob));
    break;

  case 1: /* Lambda */
    s->out = single1d_regression_propose_lambda(s->current,
						s->proposed,
						s->dataset,
						s->random,
						s->normal,
						&(s->lambda_prob));
    break;

  default:
    return NULL;
  }

  return s->proposed;
}

static double single1d_misfit(void *arg, void *state)
{
  struct single1d *s = (struct single1d *)arg;

  if (s->out == 0) {
    s->proposed_like = FLT_MAX;
  } else {
    s->proposed_like = single1d_regression_misfit(s->proposed,
						  s->dataset);
  }

  return s->proposed_like;
}

static int single1d_accept(void *arg, 
			   double current_like,
			   double proposed_like)
{
  struct single1d *s = (struct single1d *)arg;
  double u;
  int accepted;

  if (s->out == 0) {
    return 0;
  }

  u = log(s->random());

  switch (s->process) {
  case 0: /* Value */
    /* accepted = u < (current_like - proposed_like); */
    accepted = 1;
    break;

  case 1: /* Sigma */
    /* printf("sigma: %f -> %f : %f %f %f (%f < %f)\n", */
    /* 	   single1d_regression_sigma(s->current), */
    /* 	   single1d_regression_sigma(s->proposed), */
    /* 	   log(s->sigma_prob), */
    /* 	   current_like, */
    /* 	   proposed_like, */
    /* 	   u, */
    /* 	   (log(s->sigma_prob) + current_like - proposed_like)); */

    accepted = u < (log(s->lambda_prob) + current_like - proposed_like);
    break;
    
  default:
    return 0;
  }

  if (accepted) {
    single1d_regression_clone(s->proposed, s->current);
    s->current_like = s->proposed_like;


    resultset1d_accept(s->results, s->process);
  }

  return accepted;
}

static int single1d_sample(void *arg,
			   int i)
{
  struct single1d *s = (struct single1d *)arg;

  if (single1d_regression_evaluate(s->current,
				   s->dataset->xmin,
				   s->dataset->xmax,
				   s->xsamples,
				   s->v) < 0) {
    rjmcmc_error("single1d_sample: failed to evaluate current state\n");
    return -1;
  }

  resultset1d_sample(s->results,
		     i,
		     s->v);

  resultset1d_sample_misfit(s->results,
			    i,
			    s->current_like);

  resultset1d_sample_order(s->results,
			   i,
			   single1d_regression_order(s->current));
  
  resultset1d_sample_lambda(s->results,
			   i,
			   single1d_regression_lambda(s->current));
  
  do_user_callback(s);

  return 0;
}

