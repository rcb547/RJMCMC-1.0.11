#ifndef regression_h
#define regression_h

/** \file regression.h

\brief Single, 1D Partitioned and 2D Partitioned Regression

The regression functions run a MCMC simulation to fit data using
regression metrics. The likelihood is calculated as a gaussian
probability distribution, ie:

\f[

\Pr(d|m) = \frac{1}{\prod_{i=1}^n \sqrt{2 \pi \cdot \sigma_i^2}} 

\exp{\left( \sum_{i = 1}^n \frac{-\left(g\left(m\right)_i - d_i\right)^2}{2\sigma_i^2}\right)}

\f]

Where \f$d\f$ is the data, \f$m\f$ is the proposed model and \f$\sigma\f$
is the per data point error (uncorrelated). Hierarchical regression
is also supported for all of the regression methods by configuring the
lambda parameters in the 1D and 2D dataset structures (by setting
the lambda perturbation deviation to a non-zero value). This changes
the likelihood slightly to include the \f$\lambda\f$ term as follows:

\f[

\Pr(d|m) = \frac{1}{\prod_{i=1}^n \sqrt{2 \pi \cdot (\lambda \cdot \sigma_i)^2}} 

\exp{\left( \sum_{i = 1}^n \frac{-\left(g\left(m\right)_i - d_i\right)^2}{2(\lambda \cdot \sigma_i)^2}\right)}

\f]

In each of the methods, the results returned are configurable by
setting the results parameter using flags. See the ::resultset1d_t 
and ::resultset2d_t types. The default result is just the mean 
of the ensembles.

*/

#include <rjmcmc/rjmcmc_config.h>

#include <rjmcmc/resultset1d.h>
#include <rjmcmc/dataset1d.h>

#include <rjmcmc/resultset2d.h>
#include <rjmcmc/dataset2d.h>

#include <rjmcmc/rjmcmc_random.h>

typedef double (*regression1d_value_at_t)(void *state,
					  double x);

typedef void (*regression1d_cb_t)(void *state,
				  double *boundaries,
				  int nboundaries,
				  regression1d_value_at_t value_at,
				  double lambda,
				  void *user_arg);

/**

\brief Single Partition Regression

The ::single1d_regression function performs a regression simulation on
a single partition of the given dataset. It will attempt to fit trial
polynomials of order 0 to order_max and sampling acceptable
polynomials at each order.  The method of selecting the order of the
polynomial (or the probability distribution of the order of the data)
is outlined in \cite sambridge2006A.

The results are returned in a ::resultset1d_t structure that can be
interrogated with various functions in ::resultset1d.h.

For example of the use of this function, see \ref 1d/single/regression/cubic/cubic.c 

\param dataset The dataset to run the regression on. See the ::dataset1d_load_known, ::dataset1d_load_estimated, and ::dataset1d_load_fixed functions.
\param burnin The number of burn in iterations.
\param total The total number of iterations.
\param max_order The maximum order of the polynomial to allow when fitting a trial curve 
to the data.
\param xsamples The number of samples to use for discretization of the horizontal axis.
\param ysamples The number of samples to use for discretization of the vertical axis
when computing mode, median, and credible intervals.
\param credible_interval The credible interval to use for results expressed as a 
ratio, i.e. 0.95 for the 95\% credible interval.
\param random A function pointer to uniform random number generator to use. See the
::rjmcmc_uniform function.
\param normal A function pointer to the normal random number generator for generating
normally distributed random numbers with 0 mean and a standard deviation of 1. See the
::rjmcmc_normal function.
\param results A bit mask of results to store. See ::resultset1d_result_t.
\param user_callback An optional callback function to do your own sampling of the generated
curves. Set to NULL if you don't need this.
\param user_arg A user pointer to pass to the callback function.

*/
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
		    void *user_arg);
		   
/**

\brief Single Partition Regression with a custom prior

Deprecated. Used for testing.

*/
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
			       void *user_arg);

/**

\brief Single Partition Regression with direct integration

The ::single1d_regression samples the fit of polynomials drawn from a
distribution but it is possible to integrate these out directly and avoid 
the stochastic sampling process. The benefit of this is that the regression 
is performed quicker. The downside is that results such as credible intervals
are unavailable.

*/
resultset1d_t *
single1d_direct_regression(const dataset1d_t *dataset,
			   const double *fixed_prior,
			   int max_order,
			   int xsamples,
			   rjmcmc_uniform_rand_t random,
			   rjmcmc_normal_rand_t normal);


/**

\brief Multiple partition arbitrary order regression

This method fits up to max_order polynomials into a partitioned
dataset. The number of partitions and order of polynomials within each
partition is determined by the data. If higher order polynomials are
permitted (ie > 5) then there can be "ringing" near any
discontinuities in the data as a higher order polynomial may fit a
discontinuity as well as a 2 discontinuous polynomials. For this
reason care should be taken with setting the max_order parameter.

\param dataset The dataset to run the regression on. See the ::dataset1d_load_known, ::dataset1d_load_estimated, and ::dataset1d_load_fixed functions.
\param burnin The number of burn in iterations.
\param total The total number of iterations.
\param min_part The minimum number of partitions to allow.
\param max_part The maximum number of partitions to allow.
\param max_order The maximum order of the fitting polynomial. 
to the data.
\param xsamples The number of samples to use for discretization of the horizontal axis.
\param ysamples The number of samples to use for discretization of the vertical axis
when computing mode, median, and credible intervals.
\param credible_interval The credible interval to use for results expressed as a 
ratio, i.e. 0.95 for the 95\% credible interval.
\param pd The standard deviation to use when peturbing the x-position of a
partition. The new position is moved by N(0, pd).
\param random A function pointer to uniform random number generator to use. See the
::rjmcmc_uniform function.
\param normal A function pointer to the normal random number generator for generating
normally distributed random numbers with 0 mean and a standard deviation of 1. See the
::rjmcmc_normal function.
\param results A bit mask of results to store. See ::resultset1d_result_t.
\param callback A user function to call for every sample for collecting
extra statistics. Set to NULL if this is not needed.
\param user_arg The user argument to pass to the callback function.

*/

resultset1d_t *
part1d_regression(const dataset1d_t *dataset,
		  int burnin,
		  int total,
		  int min_part,
		  int max_part,
		  int max_order,
		  int xsamples,
		  int ysamples,
		  double credible_interval,
		  double pd,
		  rjmcmc_uniform_rand_t random,
		  rjmcmc_normal_rand_t normal,
		  int results,
		  regression1d_cb_t callback,
		  void *user_arg);

/**

\brief Multiple partition zeroth order regression

This method computes a regression using multiple partitions with zeroth order
curves in each partition. The level of the curve is determined from the 
data within each partition (i.e. sampled from the mean and standard deviation).
The exception is when there are too few data points within a partition, in 
which case the level of the curve is determined by uniformly sampling from
the range of the entire dataset. For the detailed theory of this method,
see the Theory section \ref part1dzeroreg.



\param dataset The dataset to run the regression on. See the ::dataset1d_load_known, ::dataset1d_load_estimated, and ::dataset1d_load_fixed functions.
\param burnin The number of burn in iterations.
\param total The total number of iterations.
\param min_part The minimum number of partitions to allow.
\param max_part The maximum number of partitions to allow.
to the data.
\param xsamples The number of samples to use for discretization of the horizontal axis.
\param ysamples The number of samples to use for discretization of the vertical axis
when computing mode, median, and credible intervals.
\param credible_interval The credible interval to use for results expressed as a 
ratio, i.e. 0.95 for the 95\% credible interval.
\param pd The standard deviation to use when peturbing the x-position of a
partition. The new position is moved by N(0, pd).
\param random A function pointer to uniform random number generator to use. See the
::rjmcmc_uniform function.
\param normal A function pointer to the normal random number generator for generating
normally distributed random numbers with 0 mean and a standard deviation of 1. See the
::rjmcmc_normal function.
\param results A bit mask of results to store. See ::resultset1d_result_t.
\param callback A user function to call for every sample for collecting
extra statistics. Set to NULL if this is not needed.
\param user_arg The user argument to pass to the callback function.

*/

resultset1d_t *
part1d_zero_regression(const dataset1d_t *dataset,
		       int burnin,
		       int total,
		       int min_part,
		       int max_part,
		       int xsamples,
		       int ysamples,
		       double credible_interval,
		       double pd,
		       rjmcmc_uniform_rand_t random,
		       rjmcmc_normal_rand_t normal,
		       int results,
		       regression1d_cb_t callback,
		       void *user_arg);

/**

\brief Multiple partition joined line segments regression

This method uses joined line segments between partitions so that a
continuous function is formed to fit the data. This is useful for 
determining changes in gradient in the data.

\param dataset The dataset to run the regression on. See the ::dataset1d_load_known, ::dataset1d_load_estimated, and ::dataset1d_load_fixed functions.
\param burnin The number of burn in iterations.
\param total The total number of iterations.
\param min_part The minimum number of partitions to allow.
\param max_part The maximum number of partitions to allow.
to the data.
\param xsamples The number of samples to use for discretization of the horizontal axis.
\param ysamples The number of samples to use for discretization of the vertical axis
when computing mode, median, and credible intervals.
\param credible_interval The credible interval to use for results expressed as a 
ratio, i.e. 0.95 for the 95\% credible interval.
\param pv The standard deviation to use when peturbing the y value of a partition. This is needed since we can no longer use the data to determine this. A rule
of thumb for selecting this value is to make it around 5 to 10 percent of the
range of the data.
\param pd The standard deviation to use when peturbing the x-position of a
partition. The new position is moved by N(0, pd).
\param random A function pointer to uniform random number generator to use. See the
::rjmcmc_uniform function.
\param normal A function pointer to the normal random number generator for generating
normally distributed random numbers with 0 mean and a standard deviation of 1. See the
::rjmcmc_normal function.
\param results A bit mask of results to store. See ::resultset1d_result_t.
\param callback A user function to call for every sample for collecting
extra statistics. Set to NULL if this is not needed.
\param user_arg The user argument to pass to the callback function.
*/

resultset1d_t *
part1d_natural_regression(const dataset1d_t *dataset,
			  int burnin,
			  int total,
			  int min_part,
			  int max_part,
			  int xsamples,
			  int ysamples,
			  double credible_interval,
			  double pv,
			  double pd,
			  rjmcmc_uniform_rand_t random,
			  rjmcmc_normal_rand_t normal,
			  int results,
			  regression1d_cb_t callback,
			  void *user_arg);

typedef double (*regression2d_value_at_t)(void *state,
					  double x,
					  double y);

typedef void (*regression2d_cb_t)(void *state,
				  double *x,
				  double *y,
				  int ncells,
				  regression2d_value_at_t value_at,
				  double lambda,
				  void *user_arg);

/*

\brief Multiple partition 2D Regression

This method performs regression on 2D datasets using zeroth order
voronoi cells. The partitions are choosen by their cell centres
rather than edges as in the 1D regression routines. 

\param dataset The dataset to run the regression on. See the ::dataset1d_load_known, ::dataset1d_load_estimated, and ::dataset1d_load_fixed functions.
\param burnin The number of burn in iterations.
\param total The total number of iterations.
\param min_part The minimum number of partitions to allow.
\param max_part The maximum number of partitions to allow.
to the data.
\param xsamples The number of samples to use for discretization of the x axis.
\param ysamples The number of samples to use for discretization of the y axis.
\param zsamples The number of samples to use for discretization of the z axis
when computing mode, median, and credible intervals.
\param credible_interval The credible interval to use for results expressed as a 
ratio, i.e. 0.95 for the 95\% credible interval.
\param pv The standard deviation to use when peturbing the y value of a partition. A rule of thumb for selecting this value is to make it around 5 to 10 percent of the range of the data.
\param pd The standard deviation to use when peturbing the (x, y) position of a
partition. The new position is moved by (N(0, pd), N(0, pd)).
\param random A function pointer to uniform random number generator to use. See the
::rjmcmc_uniform function.
\param normal A function pointer to the normal random number generator for generating
normally distributed random numbers with 0 mean and a standard deviation of 1. See the
::rjmcmc_normal function.
\param results A bit mask of results to store. See ::resultset2d_result_t.
\param callback A user function to call for every sample for collecting
extra statistics. Set to NULL if this is not needed.
\param user_arg The user argument to pass to the callback function.

*/

resultset2d_t *
part2d_regression(const dataset2d_t *dataset,
		  int burnin,
		  int total,
		  int min_part,
		  int max_part,
		  int xsamples,
		  int ysamples,
		  int zsamples,
		  double credible_interval,
		  double pv,
		  double pd,
		  rjmcmc_uniform_rand_t random,
		  rjmcmc_normal_rand_t normal,
		  int results,
		  regression2d_cb_t user_callback,
		  void *user_arg);

/** \example 1d/single/regression/cubic/cubic.c

 */

/** \example 1d/partitioned/regression/multiquad/multiquad.c

 */

/** \example 1d/partitioned/regression/multistep/multistep.c

 */

/** \example 1d/partitioned/regression/sawtooth/sawtooth.c

 */

/** \example 1d/partitioned/regression/zeromultistep/zeromultistep.c

 */

/** \example 2d/partitioned/regression/callbackex/callbackex.c

 */

/** \example 2d/partitioned/regression/square/square.c

 */

/** \example 2d/partitioned/regression/disc/disc.c

 */

#endif /* regression_h */
