#ifndef delaunay_h
#define delaunay_h

/** \file delaunay2d.h

\brief Delaunay 2D Triangulation point tracking

The Delaunay 2D routines are used for fast neighbour tracking in 2D. It is
accessed through the ::position_map2d_t interface when dynamically selected.

*/
#include <rjmcmc/bbox2d.h>

typedef struct _delaunay2d delaunay2d_t;

/** \brief Create an initial Delaunay tesselation.

Creates an initial ::delaunay2d_t object with four points representing
the bounding box and 2 triangles tesselating the region.

\param maxpoints The maximum number of points allowed to managed.
\param xmin Min x-coordinate
\param xmax Max x-coordinate
\param ymin Min y-coordinate
\param ymax Max y-coordinate
*/
delaunay2d_t *
delaunay2d_create(int maxpoints,
		  double xmin,
		  double xmax,
		  double ymin,
		  double ymax);

void 
delaunay2d_destroy(delaunay2d_t *d);

delaunay2d_t *
delaunay2d_load(const char *filename);

int
delaunay2d_save(const delaunay2d_t *d,
		const char *filename);

int
delaunay2d_clone(const delaunay2d_t *src,
		 delaunay2d_t *dst);

int 
delaunay2d_add(delaunay2d_t *d, 
	       double x,
	       double y,
	       bbox2d_t *bound);

int 
delaunay2d_delete(delaunay2d_t *d,
		  int pi,
		  bbox2d_t *bound);

int 
delaunay2d_shift_replace(delaunay2d_t *d,
			 int pi);

int 
delaunay2d_find_enclosing_triangle(const delaunay2d_t *d,
				   int t0,
				   double px,
				   double py,
				   int *pa,
				   int *pb, 
				   int *pc,
				   double *ba,
				   double *bb,
				   double *bc);

int
delaunay2d_ct_update(delaunay2d_t *d,
		     const double *point_z,
		     int n);

int 
delaunay2d_ct_value_at(delaunay2d_t *d,
		       int t0,
		       double px,
		       double py,
		       double *z);

int 
delaunay2d_nearest(delaunay2d_t *d,
		   int include_corners,
		   double px,
		   double py);

int 
delaunay2d_nearest_from(delaunay2d_t *d,
			int i,
			int include_corners,
			double px,
			double py);

int 
delaunay2d_point_of_index(const delaunay2d_t *d,
			  int i,
			  double *px,
			  double *py);

int 
delaunay2d_index_of_point(const delaunay2d_t *d,
			  double x,
			  double y);

int
delaunay2d_npoints(const delaunay2d_t *d);

int
delaunay2d_polygon_bound(const delaunay2d_t *d,
			 int i,
			 bbox2d_t *bound);

int 
delaunay2d_validate_circumcircles(const delaunay2d_t *d);

int
delaunay2d_validate_delaunay(const delaunay2d_t *d);

int
delaunay2d_validate_neighbours(const delaunay2d_t *d);

int
delaunay2d_validate_nonintersecting(const delaunay2d_t *d);

int
delaunay2d_validate_edges(const delaunay2d_t *d);

int 
delaunay2d_validate_delaunay(const delaunay2d_t *d);

void
delaunay2d_print_triangles(const delaunay2d_t *d);

void
delaunay2d_print_points(const delaunay2d_t *d);

void
delaunay2d_print_edges(const delaunay2d_t *d);

int
delaunay2d_save_geo(const delaunay2d_t *d, const char *filename);

int
delaunay2d_save_cc_geo(const delaunay2d_t *d, const char *filename);

int
triangle_circumcircle(double x1, double y1,
		      double x2, double y2,
		      double x3, double y3,
		      double *cx, double *cy,
		      double *cr2);

int 
point_in_triangle(double px,
		  double py,
		  double x1, double y1,
		  double x2, double y2,
		  double x3, double y3);

#endif /* delaunay2d_h */
