
#include <stdlib.h>

#include "rjmcmc/engine.h"

#include <rjmcmc/rjmcmc_debug.h>

int rjmcmc_engine_run(rjmcmc_engine_cb_t *cb,
		      int burn_iterations,
		      int total_iterations,
		      int sample_rate)
{
  int i;
  int p;

  if (rjmcmc_engine_init(cb,
			 burn_iterations,
			 total_iterations,
			 sample_rate) < 0) {
    return -1;
  }

  while ((i = rjmcmc_engine_step(cb)) == 1) {
  }

  if (i < 0) {
    return -1;
  }

  /* Loop completed normally */
  return 0;
}

int rjmcmc_engine_init(rjmcmc_engine_cb_t *cb,
		      int burn_iterations,
		      int total_iterations,
		      int sample_rate)
{
  if (cb == NULL) {
    rjmcmc_error("rjmcmc_engine_init: null callback\n");
    return -1;
  }

  if (burn_iterations >= total_iterations) {
    rjmcmc_error("rjmcmc_engine_init: "
		 "number of iterations must be greater than burnin\n");
    return -1;
  }

  cb->current_like = cb->initialize_state(cb->arg);
  if (cb->current_like <= 0.0) {
    rjmcmc_error("rjmcmc_engine_init: "
		 "invalid initial misfit value\n");
    return -1;
  }

  cb->burnin = burn_iterations;
  cb->total = total_iterations;
  cb->sample_rate = sample_rate;
  cb->i = 0;

  return 0;
}

int rjmcmc_engine_restart(rjmcmc_engine_cb_t *cb,
			  int burn_iterations,
			  int total_iterations,
			  int sample_rate,
			  int iteration,
			  double likelihood)
{
  int i;
  int p;

  if (cb == NULL) {
    rjmcmc_error("rjmcmc_engine_init: null callback\n");
    return -1;
  }

  if (burn_iterations >= total_iterations) {
    rjmcmc_error("rjmcmc_engine_restart: "
		 "number of iterations must be greater than burnin\n");
    return -1;
  }

  if (iteration >= total_iterations) {
    rjmcmc_error("rjmcmc_engine_restart: "
		 "already finished.");
    return -1;
  }

  cb->current_like = likelihood;

  cb->burnin = burn_iterations;
  cb->total = total_iterations;
  cb->sample_rate = sample_rate;
  cb->i = iteration;

  while ((i = rjmcmc_engine_step(cb)) == 1) {
  }

  if (i < 0) {
    return -1;
  }

  /* Loop completed normally */
  return 0;
}


int rjmcmc_engine_step(rjmcmc_engine_cb_t *cb)
{
  void *state;
  int p;
  double prop_like;
    

  p = cb->select_process(cb->arg);
  if (p < 0) {
    rjmcmc_error("rjmcmc_engine_run: "
		 "invalid process\n");
    return -1;
  }

  state = cb->perturb_process(cb->arg, p);
  if (state != NULL) {
    
    prop_like = cb->compute_misfit(cb->arg, state);
    if (prop_like <= 0.0) {
      rjmcmc_error("rjmcmc_engine_run: "
		   "invalid misfit value\n");
      return -1;
    }
    
    if (cb->accept(cb->arg, cb->current_like, prop_like)) {
      cb->current_like = prop_like;
    }
  }
  
  if (cb->sample(cb->arg, cb->i) < 0) {
    rjmcmc_error("rjmcmc_engine_run: sampling error\n");
    return -1;
  }

  cb->i ++;
  if (cb->i >= cb->total) {
    return 0;
  }

  return 1;
}

int 
rjmcmc_engine_finish(rjmcmc_engine_cb_t *cb)
{
  return 0;
}

