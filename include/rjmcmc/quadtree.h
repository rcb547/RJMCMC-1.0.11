#ifndef quadtree_h
#define quadtree_h

#include <rjmcmc/bbox2d.h>

typedef struct _quadtree quadtree_t;

quadtree_t *
quadtree_create(int maxpoints,
		int maxlevels,
		double xmin,
		double xmax,
		double ymin,
		double ymax);

void 
quadtree_destroy(quadtree_t *d);

quadtree_t *
quadtree_load(const char *filename);

int
quadtree_save(const quadtree_t *d,
	      const char *filename);

int
quadtree_clone(const quadtree_t *src,
	       quadtree_t *dst);

int 
quadtree_add(quadtree_t *d, 
	     double x,
	     double y,
	     bbox2d_t *bound);

int 
quadtree_delete(quadtree_t *d,
		int pi,
		bbox2d_t *bound);

int 
quadtree_shift_replace(quadtree_t *d,
		       int pi);

int 
quadtree_nearest(quadtree_t *d,
		 int include_corners,
		 double px,
		 double py);

int 
quadtree_point_of_index(const quadtree_t *d,
			int i,
			double *px,
			double *py);

int 
quadtree_index_of_point(const quadtree_t *d,
			double x,
			double y);

int
quadtree_npoints(const quadtree_t *d);

int
quadtree_polygon_bound(const quadtree_t *d,
		       int i,
		       bbox2d_t *bound);

#endif /* quadtree_h */
