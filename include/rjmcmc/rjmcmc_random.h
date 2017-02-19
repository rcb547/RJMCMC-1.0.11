#ifndef rjmcmc_random_h
#define rjmcmc_random_h

typedef double (*rjmcmc_uniform_rand_t)(void);
typedef double (*rjmcmc_normal_rand_t)(void);

void rjmcmc_seed(int s);
double rjmcmc_uniform(void);
double rjmcmc_normal(void);

double
rjmcmc_random_choose_double(double low,
			    double high,
			    rjmcmc_uniform_rand_t rand);

int 
rjmcmc_random_choose_int(int low,
			 int high,
			 rjmcmc_uniform_rand_t rand);

int rjmcmc_random_choose_interval(const double *cdf,
				  int n,
				  rjmcmc_uniform_rand_t rand);

#endif /* rjmcmc_random_h */
