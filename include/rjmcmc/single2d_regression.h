#ifndef single2d_regression_h
#define single2d_regression_h

typedef struct _single2d_regression 2d_regression_t;

single2d_regression_t *
single2d_regression_create(int max_order,
			   double xmin,
			   double xmax,
			   double ymin,
			   double ymax);

void
single2d_regression_destroy(single2d_regression_t *s2d);

void 
single2d_regression_clone(const single2d_regression_t *src,
			  single2d_regression_t *dst);

double
single2d_regression_misfit(single2d_regression_t *s,
			   const dataset2d_t *dataset);

int
single2d_regression_initialize(single2d_regression_t *s,
			       const dataset2d_t *datasets,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal);

int
single2d_regression_propose_value(const single2d_regression_t *current,
				  single2d_regression_t *proposed,
				  const dataset2d_t *dataset,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *value_prob);
			   
int
single2d_regression_propose_order(const single2d_regression_t *current,
				  single2d_regression_t *proposed,
				  const dataset2d_t *dataset,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *order_prob);

int
single2d_regression_propose_sigma(const single2d_regression_t *current,
				  single2d_regression_t *proposed,
				  const dataset2d_t *dataset,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *sigma_prob);


#endif /* single2d_regression_h */
