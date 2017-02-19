
#include <stdio.h>
#include <stdlib.h>

#include <rjmcmc/position_map1d.h>

int main(int argc, char *argv[])
{
  int i;
  int np;
  double u;
  double x;
  int pi;

  position_map1d_t *p;
  position_map1d_t *porig;

  int maxpoints = 200;

  p = position_map1d_create(maxpoints, -1.0, 1.0);
  if (p == NULL) {
    fprintf(stderr, "error: failed to create position map\n");
    return -1;
  }

  porig = position_map1d_create(maxpoints, -1.0, 1.0);
  if (p == NULL) {
    fprintf(stderr, "error: failed to create position map\n");
    return -1;
  }

  
  srand(0);

  for (i = 0; i < 10000; i ++) {

    u = (double)rand()/(double)RAND_MAX;
    np = position_map1d_npartitions(p);
    position_map1d_clone(p, porig);

    if ((np == 2 || u < 0.75) && (np < maxpoints)) {

      x = (double)random()/(double)RAND_MAX * 2.0 - 1.0;
      
      if (position_map1d_insert(p, x, np) < 0) {
	fprintf(stderr, "error: failed to insert point %f\n", x);
	return -1;
      }

    } else {
      
      pi = 2 + (int)((double)(np - 2) * (double)rand()/(double)RAND_MAX);
      
      x = position_map1d_position_of_index(p, pi);
      if (position_map1d_delete(p, x, pi) < 0) {
	fprintf(stderr, "error: failed to remove point\n");
	return -1;
      }
    }
  }
  
  printf("%d partitions\n",  position_map1d_npartitions(p));

  return 0;
}


