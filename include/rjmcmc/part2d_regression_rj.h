#ifndef part2d_regression_rj_h
#define part2d_regression_rj_h

#include "rjmcmc/dataset2d.h"
#include "rjmcmc/rjmcmc_random.h"

typedef struct _part2d_regression_rj part2d_regression_rj_t;

part2d_regression_rj_t *
part2d_regression_rj_create(int min_part,
			    int max_part,
			    int ndatasets,
			    double xmin, 
			    double xmax,
			    double ymin,
			    double ymax,
			    double pv,
			    double pdx,
			    double pdy,
			    double pdxy);

void
part2d_regression_rj_destroy(part2d_regression_rj_t *p);

void
part2d_regression_rj_clone(const part2d_regression_rj_t *src,
			   part2d_regression_rj_t *dst);

double
part2d_regression_rj_misfit(part2d_regression_rj_t *p,
			    const dataset2d_t **datasets,
			    int ndatasets);

int
part2d_regression_rj_evaluate(part2d_regression_rj_t *current,
			      int di,
			      double xmin,
			      double xmax,
			      int xsamples,
			      double ymin,
			      double ymax,
			      int ysamples,
			      double **z);

int 
part2d_regression_rj_initialize(part2d_regression_rj_t *p,
				const dataset2d_t **datasets,
				int ndatasets,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal);

int
part2d_regression_rj_propose_birth(const part2d_regression_rj_t *current,
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *birth_prob);

int 
part2d_regression_rj_propose_death(const part2d_regression_rj_t *current,
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *death_prob);

int 
part2d_regression_rj_propose_value(const part2d_regression_rj_t *current, 
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *value_prob);

int 
part2d_regression_rj_propose_sigma(const part2d_regression_rj_t *current,
				   part2d_regression_rj_t *proposed,
				   const dataset2d_t **datasets,
				   int ndatasets,
				   rjmcmc_uniform_rand_t random,
				   rjmcmc_normal_rand_t normal,
				   double *sigma_prob);

int 
part2d_regression_rj_propose_move(const part2d_regression_rj_t *current,
				  part2d_regression_rj_t *proposed,
				  const dataset2d_t **datasets,
				  int ndatasets,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *move_prob);

int
part2d_regression_rj_partitions(const part2d_regression_rj_t *current);

double
part2d_regression_rj_lambda(const part2d_regression_rj_t *current);

double
part2d_regression_rj_value_at(const part2d_regression_rj_t *current,
			      double x,
			      double y);

int 
part2d_regression_rj_partition_centre(const part2d_regression_rj_t *current,
				      int pi,
				      double *x,
				      double *y);



#endif /* part2d_regression_rj_h */
