
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <rjmcmc/position_map2d.h>
#include <rjmcmc/position_map2d_linear.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_defines.h>

struct _position_map2d_linear {
  int max_partitions;

  int npartitions;
  double *x;
  double *y;

  int v[3];
  double vdist2[3];
};

static double dist2(const position_map2d_linear_t *p,
		    int i,
		    double x,
		    double y);

static void add_triangle_vertex(position_map2d_linear_t *p,
				int v,
				double d2);

position_map2d_linear_t *
position_map2d_linear_create(int max_partitions,
			     double xmin,
			     double xmax,
			     double ymin,
			     double ymax)
{
  position_map2d_linear_t *p;

  p = (position_map2d_linear_t*)malloc(sizeof(position_map2d_linear_t));
  if (p == NULL) {
    return NULL;
  }

  p->max_partitions = max_partitions + 4;
  p->npartitions = 4;
  p->x = rjmcmc_create_array_1d(p->max_partitions);
  if (p->x == NULL) {
    return NULL;
  }

  p->y = rjmcmc_create_array_1d(p->max_partitions);
  if (p->y == NULL) {
    return NULL;
  }

  p->x[0] = xmin;
  p->y[0] = ymin;

  p->x[1] = xmin;
  p->y[1] = ymax;

  p->x[2] = xmax;
  p->y[2] = ymax;

  p->x[3] = xmax;
  p->y[3] = ymin;

  return p;
}

void
position_map2d_linear_destroy(position_map2d_linear_t *p)
{
  if (p != NULL) {

    rjmcmc_destroy_array_1d(p->x);
    rjmcmc_destroy_array_1d(p->y);
    
    free(p);
  }
}

void
position_map2d_linear_clone(const position_map2d_linear_t *src,
			    position_map2d_linear_t *dst)
{
  int i;

  RJMCMC_NULLCHECKVOID(src, "position_map2d_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "position_map2d_clone: null dst\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "position_map2d_clone: size mismatch\n");
  
  dst->npartitions = src->npartitions;
  for (i = 0; i < src->npartitions; i ++) {
    dst->x[i] = src->x[i];
    dst->y[i] = src->y[i];
  }
}

int
position_map2d_linear_insert(position_map2d_linear_t *p,
			     double x, double y,
			     bbox2d_t *bound)
{
  int i;

  if (p->npartitions == (p->max_partitions - 1)) {
    rjmcmc_error("position_map2d_insert: map full\n");
    return -1;
  }

  i = p->npartitions;
  p->x[i] = x;
  p->y[i] = y;

  p->npartitions ++;

  bbox2d_set(bound, p->x[0], p->x[2], p->y[0], p->y[1]);

  return i;
}

int 
position_map2d_linear_delete(position_map2d_linear_t *p,
			     int iy,
			     bbox2d_t *bound)
{
  int i;

  if (iy < 4) {
    rjmcmc_error("position_map2d_delete: can't delete corner points\n");
    return -1;
  }

  if (iy >= p->npartitions) {
    rjmcmc_error("position_map2d_delete: out of range %d >= %d\n",
		 iy,
		 p->npartitions);
    return -1;
  }

  for (i = iy + 1; i < p->npartitions; i ++) {
    p->x[i - 1] = p->x[i];
    p->y[i - 1] = p->y[i];
  }

  p->npartitions --;
  bbox2d_set(bound, p->x[0], p->x[2], p->y[0], p->y[1]);

  return 0;
}

int
position_map2d_linear_move(position_map2d_linear_t *p,
			   int iy,
			   double new_x, double new_y,
			   bbox2d_t *bound)
{
  if (iy < 4) {
    rjmcmc_error("position_map2d_move: can't move corner points\n");
    return -1;
  }

  if (iy >= p->npartitions) {
    rjmcmc_error("position_map2d_move: out of range %d >= %d\n", 
		 iy, p->npartitions);
    return -1;
  }

  p->x[iy] = new_x;
  p->y[iy] = new_y;

  bbox2d_set(bound, p->x[0], p->x[2], p->y[0], p->y[1]);
  return 0;
}

int 
position_map2d_linear_nearest(position_map2d_linear_t *p,
			      double x,
			      double y,
			      double *nx,
			      double *ny,
			      int include_boundary_points)
{
  double mind;
  int mini;
  double dx;
  double dy;
  double d;
  int i;
  int si;

  mind = DBL_MAX;
  mini = -1;

  if (include_boundary_points) {
    si = 0;
  } else {
    si = 4;
  }

  for (i = si; i < p->npartitions; i ++) {
    dx = (p->x[i] - x);
    dy = (p->y[i] - y);

    d = dx*dx + dy*dy;
    if (d < mind) {
      mind = d;
      mini = i;
    }
  }

  if (mini >= 0) {
    *nx = p->x[mini];
    *ny = p->y[mini];
  }

  return mini;
}

int 
position_map2d_linear_position_of_index(position_map2d_linear_t *p,
					int iy,
					double *x,
					double *y)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_position_of_index: null map\n");

  if (iy < 0 || iy >= p->npartitions) {
    return -1;
  }

  *x = p->x[iy];
  *y = p->y[iy];

  return 0;
}

int
position_map2d_linear_enclosing_triangle(position_map2d_linear_t *p,
					 double x,
					 double y,
					 int *ta,
					 int *tb,
					 int *tc,
					 double *ba,
					 double *bb,
					 double *bc)
{
  int i;

  for (i = 0; i < 3; i ++) {
    p->v[i] = -1;
    p->vdist2[i] = 1.0e37;
  }

  for (i = 0; i < p->npartitions; i ++) {
    add_triangle_vertex(p, i, dist2(p, i, x, y));
  }

  return 0;
}

int
position_map2d_linear_polygon_bound(position_map2d_linear_t *p,
			     int iy,
			     bbox2d_t *bound)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_polygon_bound: null map\n");

  bound->xmin = p->x[0];
  bound->xmax = p->x[2];
  bound->ymin = p->y[0];
  bound->ymax = p->y[1];

  return 0;
}

int
position_map2d_linear_validate(position_map2d_linear_t *p)
{
  return 0;
}


static double dist2(const position_map2d_linear_t *p,
		    int i,
		    double x,
		    double y)
{
  double dx;
  double dy;

  dx = p->x[i] - x;
  dy = p->y[i] - y;

  return dx*dx + dy*dy;
}

static void add_triangle_vertex(position_map2d_linear_t *p,
				int v,
				double d2)
{
  int i;
  int j;

  for (i = 0; i < 3; i ++) {
    if (d2 < p->v[i]) {
      break;
    }
  }

  if (i < 3) {
    for (j = i + 1; j < 3; j ++) {
      p->v[j] = p->v[j - 1];
      p->vdist2[j] = p->vdist2[j - 1];
    }

    p->v[i] = v;
    p->vdist2[i] = d2;
  }
}
