#ifndef curvefit_h
#define curvefit_h

/** \file curvefit.h

\brief 1D Curve Fitting routines

The ::curvefit_result_t type is used for computing best fit and
randomly perturbed polynomials for a given ::dataset1d_t .
*/

#include "rjmcmc/dataset1d.h"
#include "rjmcmc/rjmcmc_random.h"

typedef struct _curvefit_result curvefit_result_t;

/** \brief Create a new ::curvefit_result_t

\param maxorder Allow 0 .. maxorder to fit a dataset.
*/
curvefit_result_t *
curvefit_create(int maxorder);

/** \brief Free memory allocated for ::curvefit_result_t obecjt
 */
void 
curvefit_destroy(curvefit_result_t *cf);

/** \brief Computing polynomial fitting parameters

Computes the various fitting parameters for a given dataset or 
subset thereof. 

\param d The dataset
\param di First index into dataset
\param dj Last index into dataset
\param order The order of the polynomial to fit
\param cf The ::curvefit_result_t to store the data (must have been allocated
previously using ::curvefit_create
*/
int
curvefit_compute(const dataset1d_t *d,
		 int di, /* First index */
		 int dj, /* Last index */
		 int order,
		 curvefit_result_t *cf);

/** \brief Compute polynomial fitting parameters with an error scale

This function is equivalent to ::curvefit_compute with the addition
of an error scaling term lambda.

\param d The dataset
\param lambda The error scaling term, must be greater than 0.
\param di First index into dataset
\param dj Last index into dataset
\param order The order of the polynomial to fit
\param cf The ::curvefit_result_t to store the data (must have been allocated
previously using ::curvefit_create
*/
int
curvefit_compute_lambda(const dataset1d_t *d,
			double lambda,
			int di,
			int dj,
			int order,
			curvefit_result_t *cf);

/** \brief Generate a random polynomial

From a previous computed ::curvefit_result_t object, generate a random
polynomial.

\param cf The :curvefit_result_t object.
\param normal A function pointer to the Gaussian random number generator.
\param coeff Output array for polynomial coefficients.
\param coeff_len Size of the coeff array
\param prob A single double output value representing the probability of this
curve.
*/
int 
curvefit_sample(curvefit_result_t *cf,
		rjmcmc_normal_rand_t normal,
		double *coeff,
		int coeff_len,
		double *prob);

/** \brief Retrieve the mean (best fit) polynomial

From a previous computed ::curvefit_result_t object, retrieve the mean
or best fit polynomial to the data.

\param cf The :curvefit_result_t object.
\param coeff Output array for polynomial coefficients.
\param coeff_len Size of the coeff array
*/
int
curvefit_sample_mean(curvefit_result_t *cf,
		     double *coeff,
		     int coeff_len);

/** \brief Retrieve the sigma diagonals

Retrieves the square root of the diagonals of the covariance
matrix, ie standard deviation.

\param cf The :curvefit_result_t object.
\param sigma Output array for diagonal elements.
\param sigma_len sigma array length.
*/
int
curvefit_sample_sigma(curvefit_result_t *cf,
		      double *sigma,
		      int sigma_len);

/** \brief Compute the determinant of the Covariance matrix

Computes the determinant of the covariance matrix using the Cholesky 
decomposition.

\param cf The :curvefit_result_t object.
\param detCm Output double for determinant
\param order Order of the previous computation.
*/
int
curvefit_sample_detCm(curvefit_result_t *cf,
		      double *detCm,
		      int order);

/** \brief Compute the complete fit over all orders

Computes the fit of polynomials of varying order to a given dataset. Used
to determine the relative likelihood of support for polynomials of different
order to a given dataset.

\param cf The ::curvefit_result_t object
\param data The dataset
\param di First index into dataset
\param dj Last index into dataset
\param max_order The maximum order polynomial to fit (cannot be greater than 
that of th the curvefit object).
\param fixed_prior Use a fixed prior for the different orders (set to NULL to
use the automatic prior).
\param auto_z Z score for automatic prior, ie 3 = 3 standard deviations prior.
\param mean_misfit Output array of the misfit of the mean curve to the data 
for each order.
\param detCm Output array of determinant of the Covariance matrix for each
order.
\param S Covariance Matrix output
\param pk Output array of probability of each order.
\param kcdf Cummulative probability of each order.
\param mean Output Mean curves for each order.
\param sigma Output Standard deviation for each order.

*/		   
int
curvefit_evaluate_pk(curvefit_result_t *cf,
		     const dataset1d_t *data,
		     int di, 
		     int dj,
		     int max_order,
		     const double *fixed_prior,
		     double auto_z,
		     double *mean_misfit,
		     double *detCm,
		     double *prior,
		     double **S,
		     double *pk,
		     double *kcdf,
		     double **mean,
		     double **sigma);

/* Internal functions (exposed for testing) */

double **cf_L(curvefit_result_t *cf);
double **cf_Z(curvefit_result_t *cf);
double **cf_S(curvefit_result_t *cf);
double *cf_mu(curvefit_result_t *cf);

#endif /* curvefit_h */
