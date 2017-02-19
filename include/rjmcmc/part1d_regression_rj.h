#ifndef part1d_regression_rj_h
#define part1d_regression_rj_h

#include "rjmcmc/dataset1d.h"
#include "rjmcmc/rjmcmc_random.h"

typedef struct _part1d_regression_rj part1d_regression_rj_t;

part1d_regression_rj_t *
part1d_regression_rj_create(int min_part,
			    int max_part,
			    int max_order,
			    int ndatasets,
			    double xmin, 
			    double xmax,
			    double pd);

void
part1d_regression_rj_destroy(part1d_regression_rj_t *p);

void
part1d_regression_rj_clone(const part1d_regression_rj_t *src,
			   part1d_regression_rj_t *dst);

int
part1d_regression_rj_evaluate(const part1d_regression_rj_t *current,
			      int di,
			      double xmin,
			      double xmax,
			      int xsamples,
			      double *v);

double
part1d_regression_rj_value_at(const part1d_regression_rj_t *current,
			      double x);

double
part1d_regression_rj_misfit(part1d_regression_rj_t *p,
			    const dataset1d_t **datasets,
			    int ndatasets);

int 
part1d_regression_rj_initialize(part1d_regression_rj_t *p,
				const dataset1d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal);

int
part1d_regression_rj_propose_birth(const part1d_regression_rj_t *current,
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *birth_prob);

int 
part1d_regression_rj_propose_death(const part1d_regression_rj_t *current,
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *death_prob);

int 
part1d_regression_rj_propose_value(const part1d_regression_rj_t *current, 
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *value_prob);

int 
part1d_regression_rj_propose_lambda(const part1d_regression_rj_t *current,
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *lambda_prob);

int 
part1d_regression_rj_propose_move(const part1d_regression_rj_t *current,
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *move_prob);
int 
part1d_regression_rj_propose_order(const part1d_regression_rj_t *current,
				   part1d_regression_rj_t *proposed,
				   const dataset1d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *order_prob);

int
part1d_regression_rj_partitions(const part1d_regression_rj_t *current);

double
part1d_regression_rj_lambda(const part1d_regression_rj_t *current,
			    int di);

int
part1d_regression_rj_order(const part1d_regression_rj_t *current,
			   int di);

int 
part1d_regression_rj_order_sum(const part1d_regression_rj_t *current,
			       int di);			  

double
part1d_regression_rj_partition_position(const part1d_regression_rj_t *current,
					int pi);

#endif /* part1d_regression_rj_h */
