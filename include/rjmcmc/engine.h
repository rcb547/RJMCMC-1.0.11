#ifndef rjmcmc_engine_h
#define rjmcmc_engine_h

#include <rjmcmc/rjmcmc.h>

struct _rjmcmc_engine_cb {
  int burnin;
  int total;
  int sample_rate;
  int i;

  double current_like;

  /* Create an initial proposal to start from and return it's misfit */
  double (*initialize_state)(void *arg);

  /* Select the index of the process/parameter to execute/perturb */
  int (*select_process)(void *arg);

  /* Execute/Perturb the given process and return the new state */
  void *(*perturb_process)(void *arg, int proc);

  /* Compute the misfit on the state */
  double (*compute_misfit)(void *arg, void *state);

  /* Based on the likelihoods of the current and proposed misfits return
     whether the proposed state is accepted */
  int (*accept)(void *arg, 
		double current,
		double proposed);

  /* Sample the progress of the MCMC simulation */
  int (*sample)(void *arg, 
		int i);

  void *arg;
};
typedef struct _rjmcmc_engine_cb rjmcmc_engine_cb_t;
  
int rjmcmc_engine_run(rjmcmc_engine_cb_t *cb,
		      int burn_iterations,
		      int total_iterations,
		      int sample_rate);

int rjmcmc_engine_init(rjmcmc_engine_cb_t *cb,
		       int burn_iterations,
		       int total_iterations,
		       int sample_rate);

int rjmcmc_engine_restart(rjmcmc_engine_cb_t *cb,
			  int burn_iterations,
			  int total_iterations,
			  int sample_rate,
			  int iteration,
			  double likelihood);

int rjmcmc_engine_step(rjmcmc_engine_cb_t *cb);

int rjmcmc_engine_finish(rjmcmc_engine_cb_t *cb);
			     


#endif /* rjmcmc_engine_h */
