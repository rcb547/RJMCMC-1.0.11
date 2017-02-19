
#include <stdio.h>
#include <sys/time.h>

#include <rjmcmc/position_map2d.h>
#include <rjmcmc/rjmcmc_random.h>

int main(int argc, char *argv[])
{

  #define NPOINTS 100
  #define SAMPLES 1000000

  position_map2d_t *p;
  bbox2d_t bound;
  double px[NPOINTS];
  double py[NPOINTS];

  double x;
  double y;
  double nx;
  double ny;

  int li;
  int pi;

  int i;

  struct timeval starttime, stoptime;

  int sec, usec;

  p = position_map2d_create(NPOINTS + 4, -1.0, 1.0, -1.0, 1.0);
  if (p == NULL) {
    fprintf(stderr, "error: failed to create position map\n");
    return -1;
  }

  /*
   * Create points
   */
  /* printf("\n"); */
  for (i = 0; i < NPOINTS; i ++) {

    px[i] = rjmcmc_uniform() * 2.0 - 1.0;
    py[i] = rjmcmc_uniform() * 2.0 - 1.0;

    /* printf("%d %f %f\n", i, px[i], py[i]); */
    if (position_map2d_insert(p, px[i], py[i], &bound) < 0) {
      fprintf(stderr, "error: failed to add point\n");
      return -1;
    }
  }
  
  if (position_map2d_validate(p) < 0) {
    fprintf(stderr, "error: invalid position map\n");
    return -1;
  }

  /*
   * 
   */

  gettimeofday(&starttime, NULL);
  for (i = 0; i < SAMPLES; i ++) {
    x = rjmcmc_uniform() * 2.0 - 1.0;
    y = rjmcmc_uniform() * 2.0 - 1.0;

    pi = position_map2d_nearest(p,
				x,
				y,
				&nx,
				&ny,
				0);

  }
  gettimeofday(&stoptime, NULL);

  sec = (stoptime.tv_sec - starttime.tv_sec);
  usec = (stoptime.tv_usec - starttime.tv_usec);
  if (usec < 0) {
    usec = usec + 1000000;
    sec ++;
  }
  printf("%d.%06d\n", sec, usec);
    
  return 0;

}
