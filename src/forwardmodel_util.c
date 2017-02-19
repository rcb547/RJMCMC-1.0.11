
#include <math.h>
#include <float.h>

#include <rjmcmc/forwardmodel_util.h>
#include <rjmcmc/rjmcmc_debug.h>

double forwardmodel_misfit_sigma_power(const double *phi,
				       int n,
				       double sigma,
				       double r)
{
  double s;
  double rp1;

  int i;

  if (n <= 2) {
    rjmcmc_error("forwardmodel_misfit_sigma_power: n too small\n");
    return DBL_MAX;
  }

  /*
   * Initialise the sum with the 2 end points
   */
  s = 
    phi[0] * (phi[0] - r*phi[1]) + 
    phi[n - 1] * (phi[n - 1] - r*phi[n - 2]);

  rp1 = 1.0 + r*r;
  for (i = 1; i < (n - 1); i ++) {

    s += phi[i] * (-r*phi[i - 1] + rp1*phi[i] - r*phi[i + 1]);

  }

  s /= (sigma*sigma*(1.0 - r*r));
  return s;
}

double forwardmodel_log_determinant_sigma_power(int n,
						double sigma,
						double r)
{
  return (2.0 * (double)n) * log(sigma) + 
    ((double)n - 1.0) * log(1.0 - r*r);
}

