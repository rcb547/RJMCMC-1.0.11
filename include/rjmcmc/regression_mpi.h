#ifndef regression_mpi_h
#define regression_mpi_h

#include <rjmcmc/regression.h>

#if defined(HAVE_MPI_H)
/** 
 
\brief MPI version 

*/

resultset1d_t *
MPI_single1d_regression(const dataset1d_t *dataset,
			int burnin,
			int total,
			int thin,
			int max_order,
			int xsamples,
			int ysamples,
			double credible_interval,
			rjmcmc_uniform_rand_t random,
			rjmcmc_normal_rand_t normal,
			int results,
			regression1d_cb_t user_callback,
			void *user_arg,
			int mpisize,
			int mpirank,
			int root,
			MPI_Comm comm);
#endif /* HAVE_OPENMPI_MPI_H */

#endif /* regression_mpi_h */
