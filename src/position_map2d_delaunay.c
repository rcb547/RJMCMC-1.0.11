
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <rjmcmc/position_map2d_delaunay.h>
#include <rjmcmc/delaunay2d.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/rjmcmc_defines.h>

struct _position_map2d_delaunay {
  int max_partitions;

  delaunay2d_t *d;
  int si;
};

position_map2d_delaunay_t *
position_map2d_delaunay_create(int max_partitions,
			       double xmin,
			       double xmax,
			       double ymin,
			       double ymax)
{
  position_map2d_delaunay_t *p;

  p = (position_map2d_delaunay_t*)malloc(sizeof(position_map2d_delaunay_t));
  if (p == NULL) {
    return NULL;
  }

  p->max_partitions = max_partitions;
  p->d = delaunay2d_create(max_partitions + 4,
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
position_map2d_delaunay_destroy(position_map2d_delaunay_t *p)
{
  if (p != NULL) {

    delaunay2d_destroy(p->d);
    free(p);

  }
}

void
position_map2d_delaunay_clone(const position_map2d_delaunay_t *src,
		     position_map2d_delaunay_t *dst)
{
  RJMCMC_NULLCHECKVOID(src, "position_map2d_delaunay_clone: null src\n");
  RJMCMC_NULLCHECKVOID(dst, "position_map2d_delaunay_clone: null dst\n");
  RJMCMC_CONDITIONCHECKVOID(src->max_partitions != dst->max_partitions,
			    "position_map2d_delaunay_clone: size mismatch\n");

  delaunay2d_clone(src->d, dst->d);
  dst->si = src->si;
 
}

int
position_map2d_delaunay_insert(position_map2d_delaunay_t *p,
		      double x, double y,
		      bbox2d_t *bound)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_delaunay_insert: null map\n");
  
  return delaunay2d_add(p->d, x, y, bound);
}

int 
position_map2d_delaunay_delete(position_map2d_delaunay_t *p,
		      int iy,
		      bbox2d_t *bound)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_delaunay_delete: null map\n");
  
  if (delaunay2d_delete(p->d, iy, bound) < 0) {
    return -1;
  }

  p->si = 0;
  return 0;
}

int
position_map2d_delaunay_move(position_map2d_delaunay_t *p,
			     int iy,
			     double new_x, double new_y,
			     bbox2d_t *bound)
{
  bbox2d_t add_bound;

  RJMCMC_NULLCHECKINT(p, "position_map2d_delaunay_move: null map\n");

  /*
   * Delete old point
   */
  if (delaunay2d_delete(p->d, iy, bound) < 0) {
    fprintf(stderr, "position_map2d_delaunay_move: failed to delete old point\n");
    return -1;
  }

  if (position_map2d_delaunay_validate(p) < 0) {
    fprintf(stderr, "position_map2d_delaunay_move: invalid after delete (%d)\n", iy);
    delaunay2d_print_points(p->d);
    delaunay2d_print_triangles(p->d);
    exit(-1);
    return -1;
  }

  /* 
   * Add new point
   */
  if (delaunay2d_add(p->d, new_x, new_y, &add_bound) < 0) {
    fprintf(stderr, "position_map2d_delaunay_move: failed to insert new point\n");
    fprintf(stderr, "trying to add point: %g %g\n", new_x, new_y);
    return -1;
  }

  bbox2d_bound_expand(bound, &add_bound);

  /* 
   * Remap indices.
   */
  return delaunay2d_shift_replace(p->d, iy);
}

int 
position_map2d_delaunay_nearest(position_map2d_delaunay_t *p,
				double x,
				double y,
				double *nx,
				double *ny,
				int include_boundary_points)
{
  int pi;

  RJMCMC_NULLCHECKINT(p, "position_map2d_delaunay_nearest: null_map\n");


  pi = delaunay2d_nearest_from(p->d, p->si, include_boundary_points, x, y);
  if (pi < 0) {
    fprintf(stderr, "position_map2d_delaunay_nearest: failed to find nearest\n");
    return -1;
  }
  p->si = pi;


  if (delaunay2d_point_of_index(p->d, pi, nx, ny) < 0) {
    fprintf(stderr, "position_map2d_delaunay_nearest: failed to get point\n");
    return -1;
  }

  return pi;
}

int
position_map2d_delaunay_enclosing_triangle(position_map2d_delaunay_t *p,
					   double x,
					   double y,
					   int *ta,
					   int *tb,
					   int *tc,
					   double *ba,
					   double *bb,
					   double *bc)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_delaunay_enclosing_triangle: null map\n");

  return delaunay2d_find_enclosing_triangle(p->d, 0, x, y, ta, tb, tc, ba, bb, bc);
}

int 
position_map2d_delaunay_position_of_index(position_map2d_delaunay_t *p,
					  int iy,
					  double *x,
					  double *y)
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_delaunay_position_of_index: null_map\n");

  return delaunay2d_point_of_index(p->d, iy, x, y);
}

int
position_map2d_delaunay_polygon_bound(position_map2d_delaunay_t *p,
				      int iy,
				      bbox2d_t *bound) 
{
  RJMCMC_NULLCHECKINT(p, "position_map2d_delaunay_polygon_bound: null map\n");
  
  return delaunay2d_polygon_bound(p->d, iy, bound);
}

int
position_map2d_delaunay_validate(position_map2d_delaunay_t *p)
{
  /* delaunay2d_save(p->d, "validate.d"); */
  /* delaunay2d_save_geo(p->d, "validate.geo"); */
  /* delaunay2d_save_cc_geo(p->d, "validate_cc.geo"); */

  if (delaunay2d_validate_neighbours(p->d) < 0) {
    return -1;
  }

  if (delaunay2d_validate_circumcircles(p->d) < 0) {
    return -1;
  }

  if (delaunay2d_validate_delaunay(p->d) < 0) {
    return -1;
  }

  if (delaunay2d_validate_edges(p->d) < 0) {
    return -1;
  }

  return 0;
}

int
position_map2d_delaunay_save(position_map2d_delaunay_t *p,
			     const char *filename)
{
  return delaunay2d_save(p->d, filename);
}
