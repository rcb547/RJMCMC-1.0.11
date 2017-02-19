#ifndef wellrng_h
#define wellrng_h

 
/*
This is an implementation of the psuedo random number generator described
in:

  F. Panneton, P. L'Ecuyer and M. Matsumoto, 
  "Improved Long-Period Generators Based on Linear Recurrences Modulo 2", 
  submitted to ACM TOMS.

Specifically, it implements the tempered version WELL44497b and is based
on code written by the authors available from here:

  http://www.iro.umontreal.ca/~panneton/WELLRNG.html
*/


struct wellrng;
typedef struct wellrng wellrng_t;

wellrng_t *
wellrng_init_simple(unsigned int seed);

unsigned int
wellrng_seed_size(void);

wellrng_t *
wellrng_init_direct(unsigned int *seed);

double 
wellrng_sample(wellrng_t *w);

void
wellrng_destroy(wellrng_t *w);

#endif /* wellrng_h */
