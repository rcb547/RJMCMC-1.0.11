#ifndef position_map2d_delaunay_h
#define position_map2d_delaunay_h

#include <rjmcmc/bbox2d.h>

typedef struct _position_map2d_delaunay position_map2d_delaunay_t;

position_map2d_delaunay_t *
position_map2d_delaunay_create(int max_partitions,
			       double xmin,
			       double xmax,
			       double ymin,
			       double ymax);

void
position_map2d_delaunay_destroy(position_map2d_delaunay_t *p);

void
position_map2d_delaunay_clone(const position_map2d_delaunay_t *src,
			      position_map2d_delaunay_t *dst);

int
position_map2d_delaunay_insert(position_map2d_delaunay_t *p,
			       double x, double y,
			       bbox2d_t *bound);

int 
position_map2d_delaunay_delete(position_map2d_delaunay_t *p,
			       int iy,
			       bbox2d_t *bound);

int
position_map2d_delaunay_move(position_map2d_delaunay_t *p,
			     int iy,
			     double new_x, double new_y,
			     bbox2d_t *bound);

int 
position_map2d_delaunay_nearest(position_map2d_delaunay_t *p,
				double x,
				double y,
				double *nx,
				double *ny,
				int include_boundary_points);

int
position_map2d_delaunay_enclosing_triangle(position_map2d_delaunay_t *p,
					   double x,
					   double y,
					   int *ta,
					   int *tb,
					   int *tc,
					   double *ba,
					   double *bb,
					   double *bc);

int 
position_map2d_delaunay_position_of_index(position_map2d_delaunay_t *p,
					  int iy,
					  double *x,
					  double *y);

int
position_map2d_delaunay_polygon_bound(position_map2d_delaunay_t *p,
				      int iy,
				      bbox2d_t *bound);

int
position_map2d_delaunay_validate(position_map2d_delaunay_t *p);

int
position_map2d_delaunay_save(position_map2d_delaunay_t *p,
			     const char *filename);

#endif /* position_map2d_delaunay_h */
