
#ifndef voronoi2d_h
#define voronoi2d_h

typedef struct _voronoi2d voronoi2d_t;

voronoi2d_t *
voronoi2d_create(int maxpoints);

void
voronoi2d_destroy(voronoi2d_t *v);

int
voronoi2d_initialize(voronoi2d_t *v,
		     double *x,
		     double *y,
		     int npoints);

int
voronoi2d_index_of_point(voronoi2d_t *v,
			 double x,
			 double y);


#endif /* voronoi2d_h */
