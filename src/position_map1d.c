
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "rjmcmc/position_map1d.h"

#include "rjmcmc/rjmcmc_config.h"
#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_defines.h"

struct _position_map1d 
{
  int max_partitions;

  int npartitions;

  double *pos;
  int *ind;
};

position_map1d_t *
position_map1d_create(int max_partitions, double minx, double maxx)
{
  position_map1d_t *p;

  if (max_partitions < 2) {
    rjmcmc_error("position_map1d_create: invalid no. partitions\n");
    return NULL;
  }

  if (minx >= maxx) {
    rjmcmc_error("position_map1d_create: maxx must be greater than minx\n");
    return NULL;
  }

  p = malloc(sizeof(position_map1d_t));
  if (p == NULL) {
    rjmcmc_error("position_map1d_create: failed to allocate memory\n");
    return NULL;
  }

  p->max_partitions = max_partitions;
  p->npartitions = 2;

  p->pos = rjmcmc_create_array_1d(max_partitions);
  p->ind = rjmcmc_create_int_array_1d(max_partitions);

  p->pos[0] = minx;
  p->ind[0] = 0;

  p->pos[1] = maxx;
  p->ind[1] = 1;

  return p;
}

void
position_map1d_destroy(position_map1d_t *p)
{
  if (p != NULL) {
    rjmcmc_destroy_array_1d(p->pos);
    rjmcmc_destroy_int_array_1d(p->ind);
    free(p);
  }
}


void
position_map1d_clone(const position_map1d_t *src,
		     position_map1d_t *dst)
{
  int i;

  RJMCMC_NULLCHECKVOID(src, "position_map1d_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "position_map1d_clone: null dst\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "position_map1d_clone: size mismatch\n");
  
  dst->npartitions = src->npartitions;
  for (i = 0; i < src->npartitions; i ++) {
    dst->pos[i] = src->pos[i];
    dst->ind[i] = src->ind[i];
  }
}

int
position_map1d_npartitions(const position_map1d_t *p)
{
  RJMCMC_NULLCHECKINT(p, "position_map1d_npartitions: null partition map\n");

  return p->npartitions;
}


int
position_map1d_insert(position_map1d_t *p,
		      double x,
		      int iy)
{
  int i;
  int j;

  RJMCMC_NULLCHECKINT(p, "position_map1d_insert: null map\n");
  RJMCMC_CONDITIONCHECKINT(p->npartitions >= p->max_partitions,
			   "position_map1d_insert: max partitions\n");

  RJMCMC_CONDITIONCHECKINT(x <= p->pos[0],
			   "position_map1d_insert: invalid position (left)\n");
  RJMCMC_CONDITIONCHECKINT(x >= p->pos[p->npartitions - 1],
			   "position_map1d_insert: invalid position (right)\n");

  for (i = 1; i < p->npartitions; i ++) {
    if (x < p->pos[i]) {
      for (j = p->npartitions - 1; j >= i; j --) {
	p->pos[j + 1] = p->pos[j];
	p->ind[j + 1] = p->ind[j];
      }

      p->pos[i] = x;
      p->ind[i] = iy;
      p->npartitions ++;

      return 0;
    }
  }

  rjmcmc_error("position_map1d_insert: failed to find inverval\n");
  return -1;
}

int 
position_map1d_delete(position_map1d_t *p,
		      double x,
		      int iy)
{
  int i;
  int di;

  RJMCMC_NULLCHECKINT(p, "position_map1d_delete: null map\n");
  RJMCMC_CONDITIONCHECKINT(p->npartitions <= 2,
			   "position_map1d_delete: min partitions\n");

  di = -1;
  for (i = 0; i < p->npartitions; i ++) {
    if (p->pos[i] == x) {
      di = i;
      break;
    }
  }

  if (di < 0) {
    rjmcmc_error("position_map1d_delete: failed to find point\n");
    return -1;
  }

  for (i = di; i < (p->npartitions - 1); i ++) {
    p->pos[i] = p->pos[i + 1];
    p->ind[i] = p->ind[i + 1];
  }

  p->npartitions --;

  for (i = 0; i < p->npartitions; i ++) {
    if (p->ind[i] > iy) {
      p->ind[i] --;
    } else if (p->ind[i] == iy) {
      rjmcmc_error("position_map1d_delete: invalid entry\n");
      return -1;
    }
  }

  return position_map1d_valid(p);
}

int
position_map1d_move(position_map1d_t *p,
		    double x,
		    double new_x)
{
  int i;
  int mi;

  double tx;
  int ti;

  RJMCMC_NULLCHECKINT(p, "position_map1d_t: null map\n");
  
  mi = -1;

  

  for (i = 1; i < p->npartitions; i ++) {
    if (p->pos[i] == x) {
      mi = i;
      break;
    }
  }

  if (mi < 0) {
    rjmcmc_error("position_map1d_move: failed to find old point\n");
    return -1;
  }

  p->pos[mi] = new_x;

  if (new_x < x) {
    while (p->pos[mi] < p->pos[mi - 1]) {

      tx = p->pos[mi - 1];
      ti = p->ind[mi - 1];

      p->pos[mi - 1] = p->pos[mi];
      p->ind[mi - 1] = p->ind[mi];

      p->pos[mi] = tx;
      p->ind[mi] = ti;

      mi --;
    }
  } else {

    while (p->pos[mi] > p->pos[mi + 1]) {

      tx = p->pos[mi + 1];
      ti = p->ind[mi + 1];

      p->pos[mi + 1] = p->pos[mi];
      p->ind[mi + 1] = p->ind[mi];

      p->pos[mi] = tx;
      p->ind[mi] = ti;

      mi ++;
    }
  }
      
  return 0;
}

int 
position_map1d_small_move(position_map1d_t *p,
			  double x,
			  double new_x)
{
  return position_map1d_move(p,
			     x,
			     new_x);
}

int 
position_map1d_nearest(position_map1d_t *p,
		       double x,
		       double *nx)
{
  int i;
  
  int mi;
  double md;
  double d;

  RJMCMC_NULLCHECKINT(p, "position_map1d_nearest: null map\n");

  mi = -1;
  md = 1e37;

  if (p->npartitions < 1) {
    fprintf(stderr, "position_map1d_nearest: no partitions\n");
    return -1;
  }

  for (i = 0; i < p->npartitions; i ++) {
    if (x >= p->pos[i]) {
      d = x - p->pos[i];

      if (d < md) {
	md = d;
	mi = i;
      }
    }
  }

  if (mi < 0) {
    return -1;
  }
  
  if (nx != NULL) {
    *nx = p->pos[mi];
  }

  return p->ind[mi];
}

double 
position_map1d_position_of_index(position_map1d_t *p,
				 int iy)
{
  int i;

  if (p == NULL) {
    rjmcmc_error("position_map1d_position_of_index: null map\n");
    return 0.0;
  }

  for (i = 0; i < p->npartitions; i ++) {
    if (p->ind[i] == iy) {
      return p->pos[i];
    }
  }

  rjmcmc_error("position_map1d_position_of_index: failed to find interval\n");
  return 0.0;
}

int 
position_map1d_next_index(position_map1d_t *p,
			  double x)
{
  int i;
  
  RJMCMC_NULLCHECKINT(p, "position_map1d_next_index: null map\n");
  
  for (i = 1; i < p->npartitions; i ++) {
    if (p->pos[i] > x) {
      return p->ind[i];
    }
  }

  rjmcmc_error("position_map1d_next_index: failed to find interval\n");
  return -1;
}

double
position_map1d_next_position(position_map1d_t *p,
			     double x)
{
  int i;

  if (p == NULL) {
    rjmcmc_error("position_map1d_next_position: null map\n");
    return 0.0;
  }

  for (i = 1; i < p->npartitions; i ++) {
    if (p->pos[i] > x) {
      return p->pos[i];
    }
  }

  rjmcmc_error("position_map1d_next_position: failed to find interval\n");
  return 0.0;
}

int
position_map1d_predecessor_of_point(position_map1d_t *p,
				    double x)
{
  int i;
  int j;

  RJMCMC_NULLCHECKINT(p, "position_map1d_predecessor_of_point: null map\n");

  if (x >= p->pos[p->npartitions - 1]) {
    return 1;
  }

  for (i = 0, j = 1; j < p->npartitions; i ++, j ++) {
    if (p->pos[i] <= x && p->pos[j] > x) {
      return p->ind[i];
    }
  }

  return -1;
}

int
position_map1d_predecessor_of_index(position_map1d_t *p,
				    int iy)
{
  int i;

  RJMCMC_NULLCHECKINT(p, "position_map1d_predecessor_of_index: null map\n");

  if (iy == 0) {
    fprintf(stderr, 
	    "position_map1d_predecessor_of_index: invalid index\n");
    return -1;
  }

  for (i = 1; i < p->npartitions; i ++) {
    if (p->ind[i] == iy) {
      return p->ind[i - 1];
    }
  }

  return -1;
}

int 
position_map1d_successor_of_point(position_map1d_t *p,
				  double x)
{
  int i;
  int j;

  RJMCMC_NULLCHECKINT(p, "position_map1d_successor_of_point: null map\n");

  for (i = 0, j = 1; j < p->npartitions; i ++, j ++) {
    if (p->pos[i] <= x && p->pos[j] > x) {
      return p->ind[j];
    }
  }

  return -1;
}

int 
position_map1d_successor_of_index(position_map1d_t *p,
				  int iy)
{
  int i;

  RJMCMC_NULLCHECKINT(p, "position_map1d_success_of_index: null map\n");

  if (iy == 1) {
    fprintf(stderr, 
	    "position_map1d_predecessor_of_index: invalid index\n");
    return -1;
  }

  for (i = 0; i < (p->npartitions - 1); i ++) {
    if (p->ind[i] == iy) {
      return p->ind[i + 1];
    }
  }

  return -1;
}

int
position_map1d_traverse_intervals(position_map1d_t *p,
				  int (*interval_cb)(void *user_arg,
						     double xmin,
						     double xmax,
						     int iy,
						     int riy),
				  void *user_arg)
{
  int i;

  RJMCMC_NULLCHECKINT(p, "position_map1d_traverse_intervals: null map\n");
  RJMCMC_NULLCHECKINT(interval_cb, "position_map1d_traverse_intervals: null cb\n");

  for (i = 1; i < p->npartitions; i ++) {
    if (interval_cb(user_arg, 
		    p->pos[i - 1],
		    p->pos[i],
		    p->ind[i - 1],
		    p->ind[i]) < 0) {
      return -1;
    }
  }

  return 0;
}

int 
position_map1d_fill_list(position_map1d_t *p,
			 double *positions,
			 int *npartitions)
{
  int i;
  int np;

  RJMCMC_NULLCHECKINT(p, "position_map1d_fill_list: null map\n");
  RJMCMC_NULLCHECKINT(positions, "position_map1d_fill_list: null list\n");
  RJMCMC_NULLCHECKINT(npartitions, "position_map1d_fill_list: null partition count\n");

  np = p->npartitions;
  if (*npartitions < np) {
    np = *npartitions;
  }

  for (i = 0; i < np; i ++) {
    positions[i] = p->pos[i];
  }
  *npartitions = np;

  return 0;
}

void
position_map1d_dump(position_map1d_t *p)
{
  int i;

  for (i = 0; i < p->npartitions; i ++) {
    printf("%f %d\n", p->pos[i], p->ind[i]);
  }
}

int
position_map1d_valid(position_map1d_t *p)
{
  int i;
  double lastx;

  RJMCMC_NULLCHECKINT(p, "position_map1d_valid: null map\n");
  RJMCMC_CONDITIONCHECKINT(p->ind[0] != 0,
			   "position_map1d_valid: invalid first index\n");
  RJMCMC_CONDITIONCHECKINT(p->ind[p->npartitions - 1] != 1,
			   "position_map1d_valid: invalid last index\n");

  lastx = p->pos[0];
  for (i = 1; i < p->npartitions; i ++) {
    if (p->pos[i] < lastx) {
      fprintf(stderr, 
	      "position_map1d_valid: out of order %d %f %f\n",
	      i, lastx, p->pos[i]);
      return -1;
    }

    lastx = p->pos[i];
  }

  return 0;
}
