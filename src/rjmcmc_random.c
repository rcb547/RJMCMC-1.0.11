
#include <stdlib.h>
#include <math.h>

#include "rjmcmc/rjmcmc_random.h"

#include "rjmcmc/wellrng.h"

static wellrng_t *unif = NULL;
static wellrng_t *norm = NULL;

void rjmcmc_seed(int s)
{
  if (unif != NULL) {
    wellrng_destroy(unif);
  }
  unif = wellrng_init_simple(s);

  if (norm != NULL) {
    wellrng_destroy(norm);
  }
  norm = wellrng_init_simple(s);
}

double rjmcmc_uniform(void)
{
  if (unif == NULL) {
    unif = wellrng_init_simple(0);
  }

  return wellrng_sample(unif);
}

double rjmcmc_normal(void)
{
  double u;
  double v;
  double s;

  if (norm == NULL) {
    norm = wellrng_init_simple(0);
  }

  do {
    u = 2.0*wellrng_sample(norm) - 1.0;
    v = 2.0*wellrng_sample(norm) - 1.0;

    s = u*u + v*v;
  } while (s == 0.0 || s >= 1.0);

  return u * sqrt(-2.0 * log(s)/s);
}

double
rjmcmc_random_choose_double(double low,
			    double high,
			    rjmcmc_uniform_rand_t custom_rand)
{
  return low + (high - low) * custom_rand();
}

int rjmcmc_random_choose_int(int low,
			     int high,
			     rjmcmc_uniform_rand_t custom_rand)
{
  int r = low + (int)(custom_rand() * (double)(high - low + 1));

  return r;
}

int rjmcmc_random_choose_interval(const double *cdf,
				  int n,
				  rjmcmc_uniform_rand_t custom_rand)
{
  double u = custom_rand();
  int i;

  for (i = 0; i < n; i ++) {
    if (u < cdf[i]) {
      return i;
    }
  }

  return -1;
}
