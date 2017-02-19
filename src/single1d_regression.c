
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <rjmcmc/single1d_regression.h>

#include <rjmcmc/curvefit.h>
#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_debug.h>

struct _single1d_regression {

  int max_order;
  double xmin;
  double xmax;

  const double *kcdf;
  
  int order;
  double *a;

  double lambda;

  curvefit_result_t *cf;
  
};

static int update(single1d_regression_t *s,
		  int order,
		  const dataset1d_t *data,
		  rjmcmc_uniform_rand_t random,
		  rjmcmc_normal_rand_t normal);

single1d_regression_t *
single1d_regression_create(int max_order,
			   double xmin,
			   double xmax,
			   const double *kcdf)
{
  single1d_regression_t *r;

  r = (single1d_regression_t*)malloc(sizeof(single1d_regression_t));
  if (r == NULL) {
    return NULL;
  }

  r->max_order = max_order;
  r->xmin = xmin;
  r->xmax = xmax;

  r->a = rjmcmc_create_array_1d(max_order + 1);
  if (r->a == NULL) {
    return NULL;
  }

  r->cf = curvefit_create(max_order);
  if (r->cf == NULL) {
    return NULL;
  }

  r->kcdf = kcdf;

  return r;
}

void
single1d_regression_destroy(single1d_regression_t *s)
{
  if (s != NULL) {

    curvefit_destroy(s->cf);
    rjmcmc_destroy_array_1d(s->a);

    free(s);
  }
}

void
single1d_regression_clone(const single1d_regression_t *src,
			  single1d_regression_t *dst)
{
  int i;
  
  if (src == NULL) {
    rjmcmc_error("single1d_regression_clone: NULL src\n");
    return;
  }

  if (dst == NULL) {
    rjmcmc_error("single1d_regression_clone: NULL dst\n");
    return;
  }

  if (src->max_order != dst->max_order) {
    rjmcmc_error("single1d_regression_clone: max. order mismatch\n");
    return;
  }

  for (i = 0; i <= src->max_order; i ++) {
    dst->a[i] = src->a[i];
  }
  dst->order = src->order;
  dst->lambda = src->lambda;
}

double 
single1d_regression_value_at(single1d_regression_t *s,
			     double x)
{
  int i;
  double yi;
  double xi;

  xi = 1.0;
  yi = 0.0;
  for (i = 0; i <= s->order; i ++) {
    yi += s->a[i] * xi;
    xi *= x;
  }

  return yi;
}

double
single1d_regression_misfit(single1d_regression_t *s,
			   const dataset1d_t *dataset)
{
  int i;

  double y;
  double n;
  double nf;
  double dy;
  double sum;

  sum = 0.0;

  if (dataset->lambdastd > 0.0) {
    nf = s->lambda*s->lambda;
  } else {
    nf = 1.0;
  }

  for (i = 0; i < dataset->npoints; i ++) {
    y = single1d_regression_value_at(s, dataset->points[i].x);
    
    dy = y - dataset->points[i].y;
    n = dataset->points[i].n;
    sum += (dy*dy)/(2.0*n*n*nf);
  }

  return sum;
}

int
single1d_regression_evaluate(single1d_regression_t *s,
			     double xmin,
			     double xmax,
			     int nsamples,
			     double *samples)
{
  int i;

  double xi;

  for (i = 0; i < nsamples; i ++) {
    xi = xmin + (double)i * (xmax - xmin)/(double)(nsamples - 1);
    samples[i] = single1d_regression_value_at(s, xi);
  }
   
  return 0;
}

int
single1d_regression_initialize(single1d_regression_t *s,
			       const dataset1d_t *data,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal)
{
  /* Choose random order to start with */
  s->order = rjmcmc_random_choose_int(0, s->max_order, random);

  update(s, s->order, data, random, normal);

  if (data->lambdastd > 0.0) {
    s->lambda = 
      rjmcmc_random_choose_double(data->lambdamin,
				  data->lambdamax,
				  random);
  }

  return 0;
}

int
single1d_regression_propose_value(const single1d_regression_t *current,
				  single1d_regression_t *proposed,
				  const dataset1d_t *dataset,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *value_prob)
{
  single1d_regression_clone(current, proposed);

  proposed->order = rjmcmc_random_choose_interval(proposed->kcdf,
						  proposed->max_order + 1,
						  random);
  
  update(proposed, proposed->order, dataset, random, normal);

  *value_prob = 1.0;
  return -1;
}
			   
int
single1d_regression_propose_lambda(const single1d_regression_t *current,
				  single1d_regression_t *proposed,
				  const dataset1d_t *dataset,
				  rjmcmc_uniform_rand_t random,
				  rjmcmc_normal_rand_t normal,
				  double *lambda_prob)
{
  double new_s;

  single1d_regression_clone(current, proposed);
  
  new_s = proposed->lambda + normal() * dataset->lambdastd;
  if (new_s < dataset->lambdamin ||
      new_s > dataset->lambdamax) {
    return 0;
  }

  proposed->lambda = new_s;

  *lambda_prob = pow(current->lambda/proposed->lambda,
  		    (double)dataset->npoints);

  return -1;
}

static int update(single1d_regression_t *s,
		  int order,
		  const dataset1d_t *data,
		  rjmcmc_uniform_rand_t random,
		  rjmcmc_normal_rand_t normal)
{
  double curve_prob;
  
  if (curvefit_compute(data,
		       0,
		       data->npoints - 1,
		       order,
		       s->cf) < 0 ||
      curvefit_sample(s->cf,
		      normal,
		      s->a,
		      order + 1,
		      &curve_prob) < 0) {

    return -1;
  }

  return 0;
    
}

int
single1d_regression_order(const single1d_regression_t *current)
{
  return current->order;
}

double
single1d_regression_lambda(const single1d_regression_t *current)
{
  return current->lambda;
}

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
		     double **sigma)
{
  return curvefit_evaluate_pk(s->cf,
			      data,
			      0,
			      data->npoints - 1,
			      s->max_order,
			      fixedprior,
			      auto_z,
			      mean_misfit,
			      detCm,
			      autoprior,
			      S,
			      pk,
			      kcdf,
			      mean,
			      sigma);
}
