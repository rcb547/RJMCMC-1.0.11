#ifndef forwardmodel_h
#define forwardmodel_h

#include <rjmcmc/rjmcmc_config.h>

#include <rjmcmc/resultsetfm.h>
#include <rjmcmc/resultset1dfm.h>
#include <rjmcmc/resultset2dfm.h>
#include <rjmcmc/bbox2d.h>
#include <rjmcmc/rjmcmc_random.h>

/** \file forwardmodel.h

\brief Simple Forward Model Routines

The Forward Model routines allow the user to supply their own likelihood
functions where a forward model is required to transform the data an error
model other than uncorrelated gaussian noise is required.

These methods use callback functions to perform the forward model
calculations.  See the ::single_fm_likelihood_t,
::single_fm_likelihood_hierarchical_t, ::part1d_fm_likelihood_t,
::part1d_fm_hierarchical_likelihood_t, ::part2d_fm_likelihood_t and
::part2d_fm_hierarchical_likelihood_t types for the signatures of
these functions.

The examples under the 1d/single/fm, 1d/partitioned/fm and
2d/partitioned/fm directories are instructive in the use of these
methods.

*/

typedef double (*single_fm_likelihood_t)(void *user_arg,
					 int n,
					 const double *values);

typedef double (*single_fm_likelihood_hierarchical_t)(void *user_arg,
						      int nvalues,
						      const double *values,
						      int hierarchical_proposal,
						      int nhierarchicalvalues,
						      const double *hierarchicalvalues,
						      double *logdetCe);

/**

\brief Single Partition Forward Model

A single partition 1D forward model is essentially attempting to solve for
a number of parameters. A range of these parameters is supplied via the 
parameters parameter.

The user supplied callback function, lp, will be passed at each iteration
a set of parameters in an array. The callback should calculate negative 
log likelihood of the model (note that multiplicative and additive constants can be left out as these will cancel)
. For an uncorrelated gaussian, this simplifies to (and removing the additive constant):

\f[

\log({\mathrm{likelihood}}) = \sum_{i = 1}^n \frac{\left(g\left(m\right)_i - d_i\right)^2}{2 \sigma_i^2}

\f]

From the actual likelihood:

\f[

\Pr(d|m) = \frac{1}{\prod_{i=1}^n \sqrt{2 \pi \cdot \sigma_i^2}} 

\exp{\left( \sum_{i = 1}^n \frac{-\left(g\left(m\right)_i - d_i\right)^2}{2\sigma_i^2}\right)}

\f]


Where \f$g(m)_i\f$ is the data/model generated from the parameters
passed to the function, \f$d_i\f$ is the data/model to compare, and .
\f$\sigma_i\f$ is the variance of the noise or uncertainty in the
data. If the variance is not known and set to unity, the simulation
will tend to overfit the data (unless the variance is actually near
1). In this case it is recommended that the
::single_forwardmodel_hierarchical function be used instead as it will
also perturb an estimate of \f$\sigma\f$. Note that we don't need to
calculate the scaling constant for the gaussian probability in this 
case since it is cancelled out in the MCMC step. This is not the case
for the hierarchical version ::single_forwardmodel_hierarchical.

\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nparameters The number of parameters in the forward model.
\param parameters The array of parameter information.
\param lp The hierarchical likelihood function.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param samples Results such as median, mode, credible intervals are computed using a histogram of values, this parameter sets the number of bins for this histogram. A value of 100 suits most cases.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param results The results desired as a bitwise or'ing of ::resultsetfm_result_t values.

For an example of using this function, see the example \ref 1d/single/fm/simpleimage/simpleimage.c

*/

resultsetfm_t *
single_forwardmodel(int burnin,
		    int total,
		    rjmcmc_uniform_rand_t random,
		    rjmcmc_normal_rand_t normal,
		    int nparameters,
		    const forwardmodelparameter_t *parameters,
		    single_fm_likelihood_t lp,
		    void *user_arg,
		    int samples,
		    double credible_interval,
		    int results);

/**

\brief Single Partition Forward Model with Hierarchical Noise

This forward model adds a custom number of hierarchical parameters
that affect the scaling constant of the likelihood and therefore the
log(likelihood) that needs to be returned from the likelihood function.

The general gaussian probability is:

\f[

\Pr() = \frac{1}{(2 \pi)^{n/2} |C_e|} \exp{-\frac{1}{2}(g(m) - d)^T C_e^{-1} (g(m) - d)}

\f]

Where \f$ C_e \f$ is the covariance matrix that is generated from the 
hierarchical parameters. In the simplest case you would have 1 hierarchical
parameter with the diagonal of the covariance matrix set to this value.

The callback function now needs to calculate the entire likelihood
value but only when the hierarchical parameters are perturbed. To
facilitate this and to allow some efficiency gains in not having to
continuously evaluate the determinant of the covariance matrix, the
callback function only needs to return the positive value of the
exponential term above, ie \f$\frac{1}{2}(g(m) - d)^T C_e^{-1} (g(m) -
d)\f$, and the determinant of the covariance when needed. This is
indicated with the hierarchical parameter to the callback function. In
psuedo code, the callback function should follow the following
pattern:

    loglikelihood = calculate_log_likelihood();
    if (hierarchical == 1) {
        *logdetce = calculate_log_detCe();
    } 

    return loglikelihood;

\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nparameters The number of parameters in the forward model.
\param parameters The array of parameter information.
\param nhierarchicalparameters The number of hierarchical parameters to sample
\param hierarchicalparameters 
\param lp The hierarchical likelihood function. Note that this has a different signature than
the ::single_forwardmodel function parameter of the same name.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param nsamples The number of data points used in calculating the likelihood.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param results The results desired as a bitwise or'ing of ::resultsetfm_result_t values.

For an example of using this function, see the example \ref 1d/single/fm/simpleimage/simpleimage.c

*/
resultsetfm_t *
single_forwardmodel_hierarchical(int burnin,
				 int total,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 int nparameters,
				 const forwardmodelparameter_t *parameters,
				 int nhierarchicalparameters,
				 const forwardmodelparameter_t *hierarchicalparameters,
				 single_fm_likelihood_hierarchical_t lp,
				 void *user_arg,
				 int nsamples,
				 double credible_interval,
				 int results);

typedef struct _part1d_fm_likelihood_state part1d_fm_likelihood_state_t;

typedef const double *
(*part1d_fm_value_at_t)(part1d_fm_likelihood_state_t *state,
			double x);

typedef double (*part1d_fm_likelihood_t)(void *userarg,
					 int npartitions,
					 const double *partition_boundaries,
					 int nglobalparameters,
					 const double *global_parameters,
					 part1d_fm_likelihood_state_t *state,
					 part1d_fm_value_at_t value_at,
					 part1d_fm_value_at_t gradient_at);

typedef struct _part1d_fm_hierarchical_likelihood_state 
  part1d_fm_hierarchical_likelihood_state_t;

typedef const double *
(*part1d_fm_hierarchical_value_at_t)(part1d_fm_hierarchical_likelihood_state_t *state,
				     double x);

typedef double (*part1d_fm_hierarchical_likelihood_t)(
  void *userarg,
  int npartitions,
  const double *partition_boundaries,
  int nglobalparameters,
  const double *global_parameters,
  int hierarchical,
  int nhierarchicalparameters,
  const double *heirarchical_parameters,
  part1d_fm_hierarchical_likelihood_state_t *state,
  part1d_fm_hierarchical_value_at_t value_at,
  part1d_fm_hierarchical_value_at_t gradient_at,
  double *logdetce);


/**

\brief A partitioned 1d forward model.

This method performs a transdimensional mcmc forward model with zero order fitting within each
partition. For a detailed description of this method, see \cite bodinThesis.

The return value calculated should be the negative log likelihood as
in the ::single_forwardmodel method which is essentially the weighted
sum of square errors.

The callback function supplied to the lp parameter is where the forward model is implemented.
This function must match the signature of the ::part1d_fm_likelihood_t type. The parameters
passed to the this callback give you opaque access to the proposed model for which the log
likelihood needs to be computed.

\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param minpart The minimum number of partitions, this should be 2.
\param maxpart The maximum number of partitions.
\param minx The lower bound value of the x coordinate
\param maxx The upper bound value of the x coordinate.
\param xsamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param ysamples Results such as median, mode, credible intervals are computed using a histogram of values, this parameter sets the number of bins for this histogram. A value of 100 suits most cases.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param pd The standard deviation of the perturb partition boundary (move) step.
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nglobalparameters The number of global parameters.
\param global_parameters An array of global parameters (their min, max and perturbation 
standard deviations).
\param nlocalparameters The number of local parameters.
\param local_parameters An array of local parameters (their min, max and perturbation 
standard deviations).
\param lp The hierarchical likelihood function.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param results The results desired as a bitwise or'ing of ::resultset1dfm_result_t values.

For an example of using this function, see the example \ref 1d/single/fm/partitioned/functionfit/functionfit.c

*/
resultset1dfm_t *
part1d_forwardmodel(int burnin,
		    int total,
		    int minpart,
		    int maxpart,
		    double minx,
		    double maxx,
		    int xsamples,
		    int ysamples,
		    double credible_interval,
		    double pd,
		    rjmcmc_uniform_rand_t random,
		    rjmcmc_normal_rand_t normal,
		    int nglobalparameters,
		    const forwardmodelparameter_t *global_parameters,
		    int nlocalparameters,
		    const forwardmodelparameter_t *local_parameters,
		    part1d_fm_likelihood_t lp,
		    void *user_arg,
		    int results);

/** 
 
\brief Partitioned 1d forwardmodel with hierarchical parameters.

This method is the same as the ::part1d_forwardmodel except it also solves for 
hierarchical parameters that are used to construct the covariance matrix. See
\cite bodin2012A for an example application of hierarchical solving.
 
\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param minpart The minimum number of partitions, this should be 2.
\param maxpart The maximum number of partitions.
\param minx The lower bound value of the x coordinate
\param maxx The upper bound value of the x coordinate.
\param xsamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param ysamples Results such as median, mode, credible intervals are computed using a histogram of values, this parameter sets the number of bins for this histogram. A value of 100 suits most cases.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param pd The standard deviation of the perturb partition boundary (move) step.
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nglobalparameters The number of global parameters.
\param global_parameters An array of global parameters (their min, max and perturbation 
standard deviations).
\param nlocalparameters The number of local parameters.
\param local_parameters An array of local parameters (their min, max and perturbation 
standard deviations).
\param nhierarchicalparameters The number of hierarchical parameters.
\param hierarchical_parameters An array of hierarchical parameters (their min, max and perturbation 
standard deviations).
\param lp The hierarchical likelihood function.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param results The results desired as a bitwise or'ing of ::resultset1dfm_result_t values.

For an example of using this function, see the example \ref 1d/single/fm/partitioned/regression/regression.c

 */
resultset1dfm_t *
part1d_forwardmodel_hierarchical(int burnin,
				 int total,
				 int minpart,
				 int maxpart,
				 double minx,
				 double maxx,
				 int xsamples,
				 int ysamples,
				 double credible_interval,
				 double pd,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 int nglobalparameters,
				 const forwardmodelparameter_t *
				 global_parameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *
				 local_parameters,
				 int nhierarchicalparameters,
				 const forwardmodelparameter_t *
				 hierarchical_parameters,
				 part1d_fm_hierarchical_likelihood_t lp,
				 void *user_arg,
				 int results);

/**
  
\brief Natural partitioned forward model 

This method is similar to the ::part1d_forwardmodel method except that instead of fitting
zero order values in each partition, a continuous function is constructed using joined 
line segments between partition boundaries. For an example application of this method
see \cite hopcroft2007A.

\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param minpart The minimum number of partitions, this should be 2.
\param maxpart The maximum number of partitions.
\param minx The lower bound value of the x coordinate
\param maxx The upper bound value of the x coordinate.
\param xsamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param ysamples Results such as median, mode, credible intervals are computed using a histogram of values, this parameter sets the number of bins for this histogram. A value of 100 suits most cases.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param pd The standard deviation of the perturb partition boundary (move) step.
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nglobalparameters The number of global parameters.
\param global_parameters An array of global parameters (their min, max and perturbation 
standard deviations).
\param nlocalparameters The number of local parameters.
\param local_parameters An array of local parameters (their min, max and perturbation 
standard deviations).
\param lp The hierarchical likelihood function.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param results The results desired as a bitwise or'ing of ::resultset1dfm_result_t values.


 */
resultset1dfm_t *
part1d_forwardmodel_natural(int burnin,
			    int total,
			    int minpart,
			    int maxpart,
			    double xmin,
			    double xmax,
			    int xsamples,
			    int ysamples,
			    double credible_interval,
			    double pd,
			    rjmcmc_uniform_rand_t random,
			    rjmcmc_normal_rand_t normal,
			    int nglobalparameters,
			    const forwardmodelparameter_t *global_parameters,
			    int nlocalparameters,
			    const forwardmodelparameter_t *local_parameters,
			    part1d_fm_likelihood_t lp,
			    void *user_arg,
			    int results);

/**

\brief Natural partitioned forward model with hierarchical parameters.

This method is the same as the ::part1d_forwardmodel_hierarchical except it uses joined line
segments between partition boundaries rather than zero order fits within each partition
to fit the data.
 
\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param minpart The minimum number of partitions, this should be 2.
\param maxpart The maximum number of partitions.
\param minx The lower bound value of the x coordinate
\param maxx The upper bound value of the x coordinate.
\param xsamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param ysamples Results such as median, mode, credible intervals are computed using a histogram of values, this parameter sets the number of bins for this histogram. A value of 100 suits most cases.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param pd The standard deviation of the perturb partition boundary (move) step.
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nglobalparameters The number of global parameters.
\param global_parameters An array of global parameters (their min, max and perturbation 
standard deviations).
\param nlocalparameters The number of local parameters.
\param local_parameters An array of local parameters (their min, max and perturbation 
standard deviations).
\param nhierarchicalparameters The number of hierarchical parameters.
\param hierarchical_parameters An array of hierarchical parameters (their min, max and perturbation 
standard deviations).
\param lp The hierarchical likelihood function.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param results The results desired as a bitwise or'ing of ::resultset1dfm_result_t values.

*/

resultset1dfm_t *
part1d_forwardmodel_natural_hierarchical(int burnin,
					 int total,
					 int minpart,
					 int maxpart,
					 double xmin,
					 double xmax,
					 int xsamples,
					 int ysamples,
					 double credible_interval,
					 double pd,
					 rjmcmc_uniform_rand_t random,
					 rjmcmc_normal_rand_t normal,
					 int nglobalparameters,
					 const forwardmodelparameter_t *
					 global_parameters,
					 int nlocalparameters,
					 const forwardmodelparameter_t *
					 local_parameters,
					 int nhierarchicalparameters,
					 const forwardmodelparameter_t *
					 hierarchical_parameters,
					 part1d_fm_hierarchical_likelihood_t lp,
					 void *user_arg,
					 int results);

void *
part1d_forwardmodel_natural_hierarchical_init(int burnin,
					      int total,
					      int minpart,
					      int maxpart,
					      double xmin,
					      double xmax,
					      int xsamples,
					      int ysamples,
					      double credible_interval,
					      double pd,
					      rjmcmc_uniform_rand_t random,
					      rjmcmc_normal_rand_t normal,
					      int nglobalparameters,
					      const forwardmodelparameter_t *
					      global_parameters,
					      int nlocalparameters,
					      const forwardmodelparameter_t *
					      local_parameters,
					      int nhierarchicalparameters,
					      const forwardmodelparameter_t *
					      hierarchical_parameters,
					      part1d_fm_hierarchical_likelihood_t lp,
					      void *user_arg,
					      int results);

int
part1d_forwardmodel_natural_hierarchical_step(void *state);

resultset1dfm_t *
part1d_forwardmodel_natural_hierarchical_finish(void *state);

/**

\brief Cubic partitioned 1d forward model.

**Experimental** code similar to ::part1d_forwardmodel_natural except that instead of 
joining partitions with straight line segments, piecewise cubic segments are used 
with the gradients set to zero at the partition boundary. 

*/
resultset1dfm_t *
part1d_forwardmodel_zero_cubic(int burnin,
			       int total,
			       int minpart,
			       int maxpart,
			       double xmin,
			       double xmax,
			       int xsamples,
			       int ysamples,
			       double credible_interval,
			       double pd,
			       rjmcmc_uniform_rand_t random,
			       rjmcmc_normal_rand_t normal,
			       int nglobalparameters,
			       const forwardmodelparameter_t *global_parameters,
			       int nlocalparameters,
			       const forwardmodelparameter_t *local_parameters,
			       part1d_fm_likelihood_t lp,
			       void *user_arg,
			       int results);


typedef struct _part2d_fm_likelihood_state part2d_fm_likelihood_state_t;

typedef const double *
(*part2d_fm_value_at_t)(part2d_fm_likelihood_state_t *state,
			double x,
			double y);

typedef double (*part2d_fm_likelihood_t)(void *userarg,
					 int nglobalparameters,
					 const double *global_parameters,
					 part2d_fm_likelihood_state_t *state,
					 part2d_fm_value_at_t value_at,
					 part2d_fm_value_at_t gradient_at,
					 const bbox2d_t *bound);

typedef double (*part2d_fm_hierarchical_likelihood_t)(
  void *userarg,
  int nglobalparameters,
  const double *global_parameters,
  int hierarchical,
  int nhierarchicalparameters,
  const double *hierarchical_parameters,
  part2d_fm_likelihood_state_t *state,
  part2d_fm_value_at_t value_at,
  part2d_fm_value_at_t gradient_at,
  const bbox2d_t *bound,
  double *logdetce);

/**

\brief A partitioned 2d forward model.

This method performs a transdimensional mcmc forward model with zero order fitting within each
partition. In the 2d case, the partitions are defined as voronoi cells. For a detailed description of this method, see \cite bodinThesis.

The return value calculated should be the negative log likelihood as
in the ::single_forwardmodel method which is essentially the weighted
sum of square errors.

The callback function supplied to the lp parameter is where the forward model is implemented.
This function must match the signature of the ::part2d_fm_likelihood_t type. The parameters
passed to the this callback give you opaque access to the proposed model for which the log
likelihood needs to be computed.

\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param minpart The minimum number of partitions, this should be 2.
\param maxpart The maximum number of partitions.
\param initpart The initial number of partitions.
\param minx The lower bound value of the x coordinate
\param maxx The upper bound value of the x coordinate.
\param miny The lower bound value of the y coordinate
\param maxy The upper bound value of the y coordinate.
\param xsamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param ysamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param zsamples Results such as median, mode, credible intervals are computed using a histogram of values, this parameter sets the number of bins for this histogram. A value of 100 suits most cases.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param pdx The standard deviation of the perturb partition boundary (move) step in the x direction.
\param pdy The standard deviation of the perturb partition boundary (move) step in the y direction.
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nglobalparameters The number of global parameters.
\param global_parameters An array of global parameters (their min, max and perturbation 
standard deviations).
\param nlocalparameters The number of local parameters.
\param local_parameters An array of local parameters (their min, max and perturbation 
standard deviations).
\param lp The likelihood function.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param results The results desired as a bitwise or'ing of ::resultset1dfm_result_t values.

For an example of using this function, see the example \ref 2d/partitioned/fm/regression/regression.h

*/

resultset2dfm_t *
part2d_forwardmodel(int burnin,
		    int total,
		    int thin,
		    int minpart,
		    int maxpart,
		    int initpart,
		    double minx,
		    double maxx,
		    double miny,
		    double maxy,
		    int xsamples,
		    int ysamples,
		    int zsamples,
		    double credible_interval,
		    double pdx,
		    double pdy,
		    rjmcmc_uniform_rand_t random,
		    rjmcmc_normal_rand_t normal,
		    int nglobalparameters,
		    const forwardmodelparameter_t *global_parameters,
		    int nlocalparameters,
		    const forwardmodelparameter_t *local_parameters,
		    part2d_fm_likelihood_t lp,
		    void *user_arg,
		    int results);


resultset2dfm_t *
part2d_forwardmodel_restartable(const char *model_file_in,
				const char *model_file_out,
				int burnin,
				int total,
				int thin,
				int minpart,
				int maxpart,
				int initpart,
				double xmin,
				double xmax,
				double ymin,
				double ymax,
				int xsamples,
				int ysamples,
				int zsamples,
				double credible_interval,
				double pdx,
				double pdy,
				rjmcmc_uniform_rand_t random,
				rjmcmc_normal_rand_t normal,
				int nglobalparameters,
				const forwardmodelparameter_t *
				global_parameters,
				int nlocalparameters,
				const forwardmodelparameter_t *
				local_parameters,
				part2d_fm_likelihood_t lp,
				void *user_arg,
				int results);

/** 

\brief A partitioned 2d forward model with hierarchical parameters.

This method is the same as the ::part2d_forwardmodel method except it also solves for 
hierarchical parameters that are used to construct the covariance matrix.

\param burnin The number of burn in iterations.
\param total The total number of iterations (must be greater than burnin).
\param thin The number of samples to skip when averaging (0 means none).
\param minpart The minimum number of partitions, this should be 2.
\param maxpart The maximum number of partitions.
\param initpart The initial number of partitions.
\param minx The lower bound value of the x coordinate
\param maxx The upper bound value of the x coordinate.
\param miny The lower bound value of the y coordinate
\param maxy The upper bound value of the y coordinate.
\param xsamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param ysamples When producing mean curves etc, the number of discretization steps used for 
constructing the curve. A value of 100 suits most cases.
\param zsamples Results such as median, mode, credible intervals are computed using a histogram of values, this parameter sets the number of bins for this histogram. A value of 100 suits most cases.
\param credible_interval The credible interval expressed as a percentage between 0 and 1.
\param pdx The standard deviation of the perturb partition boundary (move) step in the x direction.
\param pdy The standard deviation of the perturb partition boundary (move) step in the y direction.
\param random User supplied uniform random number generator.
\param normal User supplied normal random number generator.
\param nglobalparameters The number of global parameters.
\param global_parameters An array of global parameters (their min, max and perturbation 
standard deviations).
\param nlocalparameters The number of local parameters.
\param local_parameters An array of local parameters (their min, max and perturbation 
standard deviations).
\param nhierarchicalparameters The number of hierarchical parameters.
\param hierarchical_parameters An array of hierarchical parameters (their min, max and perturbation 
standard deviations).
\param lp The likelihood function.
\param user_arg A pointer to user data that will be supplied to the likelihood function.
\param results The results desired as a bitwise or'ing of ::resultset1dfm_result_t values.

For an example of using this function, see the example \ref 2d/partitioned/fm/regression/regression.h

*/

resultset2dfm_t *
part2d_forwardmodel_hierarchical(int burnin,
				 int total,
				 int thin,
				 int minpart,
				 int maxpart,
				 int initpart,
				 double minx,
				 double maxx,
				 double miny,
				 double maxy,
				 int xsamples,
				 int ysamples,
				 int zsamples,
				 double credible_interval,
				 double pdx,
				 double pdy,
				 rjmcmc_uniform_rand_t random,
				 rjmcmc_normal_rand_t normal,
				 int nglobalparameters,
				 const forwardmodelparameter_t *
				 global_parameters,
				 int nlocalparameters,
				 const forwardmodelparameter_t *
				 local_parameters,
				 int nhierarchicalparameters,
				 const forwardmodelparameter_t *
				 hierarchical_parameters,
				 part2d_fm_hierarchical_likelihood_t lp,
				 void *user_arg,
				 int results);

/** \example 1d/single/fm/simplef/simplef.f90

 */

/** \example 1d/single/fm/simpleimage/simpleimage.c

 */

/** \example 1d/single/fm/spherefit/spherefit.c

 */

/** \example 1d/partitioned/fm/functionfit/functionfit.c

 */

/** \example 1d/partitioned/fm/functionfitf/functionfitf.f90

 */

/** \example 1d/partitioned/fm/regression/regression.c

 */


/** \example 2d/partitioned/fm/regression/regression.c

 */


#endif /* forwardmodel_h */

