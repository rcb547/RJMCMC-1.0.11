
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <rjmcmc/position_map2d_quadtree.h>
#include <rjmcmc/quadtree.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_defines.h>

static int quadtree_levels = 4;

struct _position_map2d_quadtree {
  int max_partitions;

  quadtree_t *d;
  int si;
};

position_map2d_quadtree_t *
position_map2d_quadtree_create(int max_partitions,
			       double xmin,
			       double xmax,
			       double ymin,
			       double ymax)
{
  position_map2d_quadtree_t *p;
  bbox2d_t bound;

  p = (position_map2d_quadtree_t*)malloc(sizeof(position_map2d_quadtree_t));
  if (p == NULL) {
    return NULL;
  }

  p->max_partitions = max_partitions;
  p->d = quadtree_create(max_partitions,
			 quadtree_levels,
			 xmin,
			 xmax,
			 ymin,
			 ymax);
  if (p->d == NULL) {
    return NULL;
  }

  p->si = 0;

  return p;
}

void
position_map2d_quadtree_destroy(position_map2d_quadtree_t *p)
{
  if (p != NULL) {

    quadtree_destroy(p->d);
    free(p);

  }
}

void
position_map2d_quadtree_clone(const position_map2d_quadtree_t *src,
			      position_map2d_quadtree_t *dst)
{
  RJMCMC_NULLCHECKVOID(src, "position_map2d_quadtree_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "position_map2d_quadtree_clone: null dst\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "position_map2d_quadtree_clone: size mismatch\n");

  quadtree_clone(src->d, dst->d);
  dst->si = src->si;
 
}

int
position_map2d_quadtree_insert(position_map2d_quadtree_t *p,
			       double x, double y,
			       bbox2d_t *bound)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_quadtree_insert: null map\n");
  
  return quadtree_add(p->d, x, y, bound);
}

int 
position_map2d_quadtree_delete(position_map2d_quadtree_t *p,
		      int iy,
		      bbox2d_t *bound)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_quadtree_delete: null map\n");
  
  if (quadtree_delete(p->d, iy, bound) < 0) {
    return -1;
  }

  p->si = 0;
  return 0;
}

int
position_map2d_quadtree_move(position_map2d_quadtree_t *p,
			     int iy,
			     double new_x, double new_y,
			     bbox2d_t *bound)
{
  bbox2d_t add_bound;

  RJMCMC_NULLCHECKINT(p, "position_map2d_quadtree_move: null map\n");

  /*
   * Delete old point
   */
  if (quadtree_delete(p->d, iy, bound) < 0) {
    return -1;
  }

  /* 
   * Add new point
   */
  if (quadtree_add(p->d, new_x, new_y, &add_bound) < 0) {
    return -1;
  }

  bbox2d_bound_expand(bound, &add_bound);

  /* 
   * Remap indices.
   */
  return quadtree_shift_replace(p->d, iy);
}

int 
position_map2d_quadtree_nearest(position_map2d_quadtree_t *p,
				double x,
				double y,
				double *nx,
				double *ny,
				int include_boundary_points)
{
  int pi;

  RJMCMC_NULLCHECKINT(p, "position_map2d_quadtree_nearest: null_map\n");


  pi = quadtree_nearest(p->d, include_boundary_points, x, y);
  if (pi < 0) {
    fprintf(stderr, "position_map2d_quadtree_nearest: failed to find nearest\n");
    return -1;
  }


  if (quadtree_point_of_index(p->d, pi, nx, ny) < 0) {
    fprintf(stderr, "position_map2d_quadtree_nearest: failed to get point\n");
    return -1;
  }

  return pi;
}

int
position_map2d_quadtree_enclosing_triangle(position_map2d_quadtree_t *p,
					   double x,
					   double y,
					   int *ta,
					   int *tb,
					   int *tc,
					   double *ba,
					   double *bb,
					   double *bc)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_quadtree_enclosing_triangle: null map\n");

  return -1;
}

int 
position_map2d_quadtree_position_of_index(position_map2d_quadtree_t *p,
					  int iy,
					  double *x,
					  double *y)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_quadtree_position_of_index: null_map\n");

  return quadtree_point_of_index(p->d, iy, x, y);
}

int
position_map2d_quadtree_polygon_bound(position_map2d_quadtree_t *p,
				      int iy,
				      bbox2d_t *bound) 
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_quadtree_polygon_bound: null map\n");

  return quadtree_polygon_bound(p->d, iy, bound);
}

int
position_map2d_quadtree_validate(position_map2d_quadtree_t *p)
{
  return 0;
}

