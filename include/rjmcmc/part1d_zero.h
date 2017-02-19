#ifndef part1d_zero_h
#define part1d_zero_h

#include "rjmcmc/dataset1d.h"
#include "rjmcmc/rjmcmc_random.h"

typedef struct _part1d_zero part1d_zero_t;

part1d_zero_t *
part1d_zero_create(int min_part,
		   int max_part,
		   int ndatasets,
		   double xmin, 
		   double xmax,
		   double pd);

void
part1d_zero_destroy(part1d_zero_t *p);

void
part1d_zero_clone(const part1d_zero_t *src,
			part1d_zero_t *dst);

double 
part1d_zero_value_at(const part1d_zero_t *current,
		     double x);

int
part1d_zero_evaluate(const part1d_zero_t *current,
		     int di,
		     double xmin,
		     double xmax,
		     int xsamples,
		     double *v);

double
part1d_zero_misfit(part1d_zero_t *p,
		   const dataset1d_t **datasets,
		   int ndatasets);

int 
part1d_zero_initialize(part1d_zero_t *p,
		       const dataset1d_t **datasets,
		       int ndatasets,
		       rjmcmc_uniform_rand_t random,
		       rjmcmc_normal_rand_t normal);

int
part1d_zero_propose_birth(const part1d_zero_t *current,
			  part1d_zero_t *proposed,
			  const dataset1d_t **datasets,
			  int ndatasets,
			  rjmcmc_uniform_rand_t random,
			  rjmcmc_normal_rand_t normal,
			  double *birth_prob);

int 
part1d_zero_propose_death(const part1d_zero_t *current,
			  part1d_zero_t *proposed,
			  const dataset1d_t **datasets,
			  int ndatasets,
			  rjmcmc_uniform_rand_t random,
			  rjmcmc_normal_rand_t normal,
			  double *death_prob);

int 
part1d_zero_propose_value(const part1d_zero_t *current, 
			  part1d_zero_t *proposed,
			  const dataset1d_t **datasets,
			  int ndatasets,
			  rjmcmc_uniform_rand_t random,
			  rjmcmc_normal_rand_t normal,
			  double *value_prob);

int 
part1d_zero_propose_lambda(const part1d_zero_t *current,
			   part1d_zero_t *proposed,
			   const dataset1d_t **datasets,
			   int ndatasets,
			   rjmcmc_uniform_rand_t random,
			   rjmcmc_normal_rand_t normal,
			   double *sigma_prob);

int 
part1d_zero_propose_move(const part1d_zero_t *current,
			       part1d_zero_t *proposed,
			       const dataset1d_t **datasets,
			       int ndatasets,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal,
			       double *move_prob);


int
part1d_zero_partitions(const part1d_zero_t *current);

double
part1d_zero_partition_position(const part1d_zero_t *current,
			       int pi);

double
part1d_zero_lambda(const part1d_zero_t *current,
		   int di);
			  

#endif /* part1d_zero_h */
