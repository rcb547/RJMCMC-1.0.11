#ifndef single1d_regression_h
#define single1d_regression_h

#include <rjmcmc/dataset1d.h>
#include <rjmcmc/rjmcmc_random.h>

typedef struct _single1d_regression single1d_regression_t;

single1d_regression_t *
single1d_regression_create(int max_order,
			   double xmin,
			   double xmax,
			   const double *kcdf);

void
single1d_regression_destroy(single1d_regression_t *s1d);

void 
single1d_regression_clone(const single1d_regression_t *src,
			  single1d_regression_t *dst);

double 
single1d_regression_value_at(single1d_regression_t *s,
			     double x);

double
single1d_regression_misfit(single1d_regression_t *s,
			   const dataset1d_t *dataset);

int
single1d_regression_evaluate(single1d_regression_t *s,
			     double xmin,
			     double xmax,
			     int nsamples,
			     double *samples);

int
single1d_regression_initialize(single1d_regression_t *s,
			       const dataset1d_t *dataset,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal);

int
single1d_regression_propose_value(const single1d_regression_t *current,
				  single1d_regression_t *proposed,
				  const dataset1d_t *dataset,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *value_prob);
			   
int
single1d_regression_propose_lambda(const single1d_regression_t *current,
				  single1d_regression_t *proposed,
				  const dataset1d_t *dataset,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *lambda_prob);

int
single1d_regression_order(const single1d_regression_t *current);

double
single1d_regression_lambda(const single1d_regression_t *current);

int 
single1d_compute_mean_misfit(single1d_regression_t *s,
			     int order,
			     const dataset1d_t *data,
			     double *mean_misfit,
			     double *detCm,
			     double *autoprior);

int
single1d_evaluate_pk(single1d_regression_t *s,
		     const dataset1d_t *data,
		     const double *fixedprior,
		     double auto_z,
		     double *mean_misfit,
		     double *detCm,
		     double *autoprior,
		     double **S,
		     double *pk,
		     double *kcdf,
		     double **mean,
		     double **sigma);
		      
		      

#endif /* single1d_regression_h */
