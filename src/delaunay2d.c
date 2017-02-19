
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <rjmcmc/delaunay2d.h>
#include <rjmcmc/rjmcmc_defines.h>

#define TRIANGLE_INCR 1024
#define DLIST_INCR 16
#define MAXEDGES 32

//#define DEBUG_FLIPS 1

/* static const double BARYCENTRE_EPSILON = 1.0e-12; */
/* static const double BARYCENTRE_ONE_EPSILON = 1.0000000001; */

static const double BARYCENTRE_EPSILON = 0.0;
static const double BARYCENTRE_ONE_EPSILON = 1.0;

struct _point {
  double x;
  double y;

  /* Cached distance to point */
  int ci;
  double pd2;

  /* Partial derivatives */
  double z;
  double dx;
  double dy;
};
typedef struct _point point_t;

struct _triangle {
  int v[3]; // int a, b, c;
  int n[3]; // int nab, nbc, nca;

  double detT;
  double cx;
  double cy;
  double cr2;

  double ct[20];
};
typedef struct _triangle triangle_t;

struct _edgelist {
  int n;
  int edge[MAXEDGES];
};
typedef struct _edgelist edgelist_t;

struct _delaunay2d {

  int maxpoints;

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double d2max;

  point_t *p;
  edgelist_t *e;
  int np;

  triangle_t *t;
  int nt;
  int st;

  int *nlist;
  int *tlist;
  int *vlist;
  int listsize;

  int *ftlist;
  int *felist;
  int nf;
  int sf;
  
  int ci;
};
  
static int add_point(delaunay2d_t *d,
		     double x,
		     double y);

static int add_triangle(delaunay2d_t *d,
			int a, 
			int b,
			int c,
			int nab, 
			int nbc, 
			int nca);

static int add_edge(delaunay2d_t *d,
		    int p,
		    int np);

static void triangle_copy(const triangle_t *src,
			  triangle_t *dst);

static int triangle_replace_neighbour(delaunay2d_t *d,
				      int ti,
				      int oldn,
				      int newn);

static int triangle_verify_neighbour(delaunay2d_t *d,
				     int ti,
				     int ei);

static int triangle_update(delaunay2d_t *d,
			   int ti);

static int validate_triangle_pair(delaunay2d_t *d,
				  int ti,
				  int edge,
				  int bidirectional);

static double compute_detT(const delaunay2d_t *d,
			   int a,
			   int b,
			   int c);

static int barycentre_of_point(const delaunay2d_t *d,
			       int ti,
			       double x,
			       double y,
			       double *ba,
			       double *bb,
			       double *bc);

static int circumcircle_of_triangle(const delaunay2d_t *d,
				    int ti,
				    double *cx,
				    double *cy,
				    double *r2);

static int neighbour_index(const delaunay2d_t *d,
			   int neighbour,
			   int tri);

static int *resize_int_array(int *a, int size, int new_size);

static int delete_triangle(delaunay2d_t *d,
			   int ti,
			   int *new_ti);

static int delete_point(delaunay2d_t *d,
			int pi);

static int delete_edge(delaunay2d_t *d,
		       int p,
		       int np);

static int shift_down_edge(delaunay2d_t *d,
			   int p);

static int flip_edge(delaunay2d_t *d,
		     int e1a, 
		     int e1b,
		     int e2a,
		     int e2b);

static void rotate_int_array(int *a, int size, int s);

static int delaunay2d_fill_hole(delaunay2d_t *d, 
				int nd);

static int clear_flip(delaunay2d_t *d);
static int queue_flip(delaunay2d_t *d, 
		      int ti,
		      int ei);
static int update_flip(delaunay2d_t *d,
		       int old_ti,
		       int new_ti);

static int fill_fan_about_point(delaunay2d_t *d,
				int pi,
				int ti,
				bbox2d_t *bound);

static double dist2(delaunay2d_t *d,
		    int pi,
		    double px,
		    double py);


int
delaunay2d_all_visible(const delaunay2d_t *d,
		       int *vlist,
		       int nvertices);

delaunay2d_t *
delaunay2d_create(int maxpoints,
		  double xmin,
		  double xmax,
		  double ymin,
		  double ymax)
{
  delaunay2d_t *d;

  RJMCMC_CONDITIONCHECKPTR(maxpoints < 4, "delaunay2d_create: maxpoints must be at least 4\n");
  RJMCMC_CONDITIONCHECKPTR(xmin >= xmax, "delaunay2d_create: x range invalid\n");
  RJMCMC_CONDITIONCHECKPTR(ymin >= ymax, "delaunay2d_create: y range invalid\n");
  
  d = malloc(sizeof(delaunay2d_t));
  RJMCMC_NULLCHECKPTR(d, "delaunay2d_create: failed to allocate memory\n");

  d->maxpoints = maxpoints;

  d->xmin = xmin;
  d->xmax = xmax;
  d->ymin = ymin;
  d->ymax = ymax;

  d->d2max = (xmax - xmin)*(xmax - xmin) + (ymax - ymin)*(ymax - ymin);

  d->p = malloc(sizeof(point_t) * maxpoints);
  RJMCMC_NULLCHECKPTR(d->p, "delaunay2d_create: failed to allocate memory for points\n");
  memset(d->p, 0, sizeof(point_t) * maxpoints);

  d->e = malloc(sizeof(edgelist_t) * maxpoints);
  RJMCMC_NULLCHECKPTR(d->e, "delaunay2d_create: failed to allocate memory for edge list\n");
  d->np = 0;

  d->t = malloc(sizeof(triangle_t) * TRIANGLE_INCR);
  RJMCMC_NULLCHECKPTR(d->t, "delaunay2d_create: failed to allocate memory for triangles\n");
  d->st = TRIANGLE_INCR;
  d->nt = 0;

  RJMCMC_INTCHECKPTR(add_point(d, xmin, ymin), "delaunay2d_create: failed to add point\n");
  RJMCMC_INTCHECKPTR(add_point(d, xmin, ymax), "delaunay2d_create: failed to add point\n");
  RJMCMC_INTCHECKPTR(add_point(d, xmax, ymax), "delaunay2d_create: failed to add point\n");
  RJMCMC_INTCHECKPTR(add_point(d, xmax, ymin), "delaunay2d_create: failed to add point\n");

  RJMCMC_INTCHECKPTR(add_triangle(d, 0, 1, 2, -1, -1, 1), "delaunay2d_create: failed to add triangle\n");
  RJMCMC_INTCHECKPTR(add_triangle(d, 0, 2, 3, 0, -1, -1), "delaunay2d_create: failed to add triangle\n");

  d->listsize = DLIST_INCR;
  d->nlist = malloc(sizeof(int) * d->listsize);
  d->tlist = malloc(sizeof(int) * d->listsize);
  d->vlist = malloc(sizeof(int) * d->listsize);

  RJMCMC_NULLCHECKPTR(d->nlist, "delaunay2d_create: failed to allocate memory for nlist\n");
  RJMCMC_NULLCHECKPTR(d->tlist, "delaunay2d_create: failed to allocate memory for tlist\n");
  RJMCMC_NULLCHECKPTR(d->vlist, "delaunay2d_create: failed to allocate memory for vlist\n");

  d->sf = DLIST_INCR;
  d->ftlist = malloc(sizeof(int) * d->sf);
  d->felist = malloc(sizeof(int) * d->sf);
  d->nf = 0;

  d->ci = 1;
  
  return d;
}
			
void delaunay2d_destroy(delaunay2d_t *d)
{
  if (d != NULL) {
    free(d->p);
    free(d->e);
    free(d->t);

    free(d->nlist);
    free(d->tlist);
    free(d->vlist);

    free(d->ftlist);
    free(d->felist);

    free(d);
  }
}

delaunay2d_t *
delaunay2d_load(const char *filename)
{
  FILE *fp;
  delaunay2d_t *d;
  int maxpoints;
  int npoints;

  double xmin;
  double xmax;
  double ymin;
  double ymax;

  double x;
  double y;

  int a, b, c;
  int size_triangles;
  int ntriangles;
  int nab, nbc, nca;

  int i;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "delaunay2d_load: failed to open file\n");
    return NULL;
  }

  if (fscanf(fp, "%d %d\n", &maxpoints, &npoints) != 2) {
    fprintf(stderr, "delaunay2d_load: failed to parse points header\n");
    return NULL;
  }

  if (fscanf(fp, "%lf %lf %lf %lf\n", &xmin, &xmax, &ymin, &ymax) != 4) {
    fprintf(stderr, "delaunay2d_load: failed to parse bbox header\n");
    return NULL;
  }
  
  d = delaunay2d_create(maxpoints, xmin, xmax, ymin, ymax);
  if (d == NULL) {
    fprintf(stderr, "delaunay2d_load: failed to create delaunay structure\n");
    return NULL;
  }

  for (i = 0; i < npoints; i ++) {
    if (fscanf(fp, "%lf %lf\n", &x, &y) != 2) {
      fprintf(stderr, "delaunay2d_load: failed to parse point\n");
      return NULL;
    }

    d->p[i].x = x;
    d->p[i].y = y;
  }
  d->np = npoints;
  

  if (fscanf(fp, "%d %d\n", &size_triangles, &ntriangles) != 2) {
    fprintf(stderr, "delaunay2d_load: failed to parse triangles header\n");
    return NULL;
  }

  d->nt = ntriangles;
  for (i = 0; i < ntriangles; i ++) {
    if (fscanf(fp, "%d %d %d %d %d %d\n", &a, &b, &c, &nab, &nbc, &nca) != 6) {
      fprintf(stderr, "delaunay2d_load: failed to parse triangle\n");
      return NULL;
    }

    d->t[i].v[0] = a;
    d->t[i].v[1] = b;
    d->t[i].v[2] = c;

    d->t[i].n[0] = nab;
    d->t[i].n[1] = nbc;
    d->t[i].n[2] = nca;

    if (triangle_update(d, i) < 0) {
      fprintf(stderr, "delaunay2d_load: failed to update triangle (%d %d %d)\n", a, b, c);
      return NULL;
    }
  }
  
  fclose(fp);
  return d;
}

int
delaunay2d_save(const delaunay2d_t *d,
		const char *filename)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "delaunay2d_save: failed to create file\n");
    return -1;
  }

  fprintf(fp, "%d %d\n", d->maxpoints, d->np);
  fprintf(fp, "%f %f %f %f\n", d->xmin, d->xmax, d->ymin, d->ymax);

  for (i = 0; i < d->np; i ++) {
    fprintf(fp, "%f %f\n", d->p[i].x, d->p[i].y);
  }

  fprintf(fp, "%d %d\n", d->st, d->nt);
  for (i = 0; i < d->nt; i ++) {
    fprintf(fp, "%d %d %d %d %d %d\n",
	    d->t[i].v[0], d->t[i].v[1], d->t[i].v[2], 
	    d->t[i].n[0], d->t[i].n[1], d->t[i].n[2]);
  }

  fclose(fp);
  return 0;
}

int
delaunay2d_clone(const delaunay2d_t *src,
		 delaunay2d_t *dst)
{
  int i;
  int j;

  RJMCMC_CONDITIONCHECKINT(src == NULL, "delaunay2d_clone: null src\n");
  RJMCMC_CONDITIONCHECKINT(dst == NULL, "delaunay2d_clone: null dst\n");
  RJMCMC_CONDITIONCHECKINT(src->maxpoints != dst->maxpoints, "delaunay2d_clone: maxpoints mismatch\n");

  dst->np = src->np;
  for (i = 0; i < src->np; i ++) {
    dst->p[i].x = src->p[i].x;
    dst->p[i].y = src->p[i].y;

    dst->e[i].n = src->e[i].n;
    for (j = 0; j < src->e[i].n; j ++) {
      dst->e[i].edge[j] = src->e[i].edge[j];
    }
  }

  dst->nt = src->nt;
  if (dst->st < src->st) {
    free(dst->t);
    dst->t = malloc(sizeof(triangle_t) * src->st);
    RJMCMC_CONDITIONCHECKINT(dst->t == NULL, "delaunay2d_clone: failed to resize triangles\n");
    dst->st = src->st;
  }
  
  for (i = 0; i < src->nt; i ++) {
    dst->t[i].v[0] = src->t[i].v[0];
    dst->t[i].v[1] = src->t[i].v[1];
    dst->t[i].v[2] = src->t[i].v[2];
    
    dst->t[i].n[0] = src->t[i].n[0];
    dst->t[i].n[1] = src->t[i].n[1];
    dst->t[i].n[2] = src->t[i].n[2];

    dst->t[i].detT = src->t[i].detT;

    dst->t[i].cx = src->t[i].cx;
    dst->t[i].cy = src->t[i].cy;
    dst->t[i].cr2 = src->t[i].cr2;
  }

  return 0;
}

int delaunay2d_add(delaunay2d_t *d, 
		   double x,
		   double y,
		   bbox2d_t *bound)
{
  int ti;

  int pa;
  int pb;
  int pc;

  double ba;
  double bb;
  double bc;

  int a;
  int b;
  int c;
   
  /* int nab; */
  int nbc;
  int nca;

  int nab1;
  int nbc1;
  int nca1;

  int nab2;
  int nbc2;
  int nca2;

  int pi;

  int tbc;
  int tca;

  int oedge;

  int a1, b1, c1;
  int a2, b2, c2;
  int n1, n2, n3, n4;
  int n;
  int t1, t2, t3, t4;

  int ei;
  
  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_add: null delaunay\n");
  RJMCMC_CONDITIONCHECKINT(bound == NULL, "delaunay2d_add: null bound\n");

  ti = delaunay2d_find_enclosing_triangle(d,
					  0,
					  x,
					  y,
					  &pa,
					  &pb, 
					  &pc,
					  &ba,
					  &bb,
					  &bc);
  if (ti < 0) {
    fprintf(stderr, "delaunay2d_add: failed to find triangle to add point\n");
    return -1;
  }

  /* fprintf(stderr, "delaunay2d_add: bc: %d %g %g %g (%f %f) %d\n", ti, ba, bb, bc, x, y, (int)(fabs(bb) < BARYCENTRE_EPSILON)); */
  /*   /\* fprintf(stderr, "\nPre add\n"); *\/ */
  /*   /\* delaunay2d_print_points(d); *\/ */
  /*   /\* delaunay2d_print_triangles(d); *\/ */

  bbox2d_initialize(bound, x, y);

  /* Determine edge index of the new point, -1 mean within the triangle and not on an edge */
  ei = -1; 
  if (fabs(ba) <= BARYCENTRE_EPSILON) {
    ei = 1;
  } else if (fabs(bb) <= BARYCENTRE_EPSILON) {
    ei = 2;
  } else if (fabs(bc) <= BARYCENTRE_EPSILON) {
    ei = 0;
  }

  if (ei >= 0) {

    n = d->t[ti].n[ei];

    oedge = neighbour_index(d,
			    ti,
			    n);
    if (oedge < 0) {
      RJMCMC_CONDITIONCHECKINT(d->t[ti].n[1] >= 0, 
		  "delaunay2d_add: failed to find neighbour edge\n");
	
      pi = d->np;
      t2 = d->nt;

      /* No neighbour, so split triangle in half */
      RJMCMC_CONDITIONCHECKINT(add_point(d, x, y) < 0, 
		  "delaunay2d_add: failed to add point\n");

      bbox2d_expand(bound, d->p[d->t[ti].v[0]].x, d->p[d->t[ti].v[0]].y);
      bbox2d_expand(bound, d->p[d->t[ti].v[1]].x, d->p[d->t[ti].v[1]].y);
      bbox2d_expand(bound, d->p[d->t[ti].v[2]].x, d->p[d->t[ti].v[2]].y);

      a1 = d->t[ti].v[ei];
      b1 = pi;
      c1 = d->t[ti].v[(ei + 2) % 3];

      nab1 = -1;
      nbc1 = t2;
      nca1 = d->t[ti].n[(ei + 2) % 3];

      a2 = pi;
      b2 = d->t[ti].v[(ei + 1) % 3];
      c2 = d->t[ti].v[(ei + 2) % 3];
      
      nab2 = -1;
      nbc2 = d->t[ti].n[(ei + 1) % 3];
      nca2 = ti;
 

      RJMCMC_CONDITIONCHECKINT(delete_edge(d, d->t[ti].v[ei], d->t[ti].v[(ei + 1) % 3]), 
			       "delaunay2d_add: failed to remove edge\n");

      /* Replace containing triangle (a, b, np) */
      d->t[ti].v[0] = a1;
      d->t[ti].v[1] = b1;
      d->t[ti].v[2] = c1;

      d->t[ti].n[0] = nab1;
      d->t[ti].n[1] = nbc1;
      d->t[ti].n[2] = nca1;
      
      RJMCMC_CONDITIONCHECKINT(triangle_update(d, ti) < 0,
		  "delaunay2d_add: failed to update new triangle 1\n");
      
      RJMCMC_CONDITIONCHECKINT(add_edge(d, d->t[ti].v[ei], pi),
			       "delaunay2d_add: failed to add edge ba = 0, -1\n");

      /* Add new triangle */
      RJMCMC_CONDITIONCHECKINT(add_triangle(d, a2, b2, c2, nab2, nbc2, nca2) < 0,
		  "delaunay2d_add: failed to add triangle a, p, c\n");

      if (nbc2 >= 0) {
	RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, nbc2, ti, t2) < 0,
		    "delaunay2d_add: failed to replace neighbour\n");
      }

      /* Flip edges as required */
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, ti, 2, 0) < 0,
		  "delaunay2d_add: failed to flip triangle\n");
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, t2, 1, 0) < 0,
		  "delaunay2d_add: failed to flip new triangle\n");

    } else {

      /* Edge is shared */

      a1 = d->t[ti].v[ei];
      b1 = d->t[ti].v[(ei + 1) % 3];
      c1 = d->t[ti].v[(ei + 2) % 3];

      a2 = d->t[n].v[oedge];
      b2 = d->t[n].v[(oedge + 1) % 3];
      c2 = d->t[n].v[(oedge + 2) % 3];

      bbox2d_expand(bound, d->p[a1].x, d->p[a1].y);
      bbox2d_expand(bound, d->p[b1].x, d->p[b1].y);
      bbox2d_expand(bound, d->p[c1].x, d->p[c1].y);
      bbox2d_expand(bound, d->p[c2].x, d->p[c2].y);

      n1 = d->t[ti].n[(ei + 2) % 3];
      n2 = d->t[n].n[(oedge + 1) % 3];
      n3 = d->t[n].n[(oedge + 2) % 3];
      n4 = d->t[ti].n[(ei + 1) % 3];

      t1 = ti;
      t2 = n;
      t3 = d->nt;
      t4 = t3 + 1;

      pi = d->np;
      RJMCMC_CONDITIONCHECKINT(add_point(d, x, y) < 0, 
		  "delaunay2d_add: failed to add point\n");

      /* First triangle (replace) */
      d->t[ti].v[0] = a1;
      d->t[ti].v[1] = pi;
      d->t[ti].v[2] = c1;
      d->t[ti].n[0] = t2;
      d->t[ti].n[1] = t4;
      d->t[ti].n[2] = n1;

      RJMCMC_CONDITIONCHECKINT(triangle_update(d, ti) < 0,
		  "delaunay2d_add: failed to update triangle 1\n");

      /* 2nd triangle (replace) */
      d->t[t2].v[0] = c2;
      d->t[t2].v[1] = pi;
      d->t[t2].v[2] = b2;
      d->t[t2].n[0] = t3;
      d->t[t2].n[1] = t1;
      d->t[t2].n[2] = n2;

      RJMCMC_CONDITIONCHECKINT(triangle_update(d, t2) < 0,
		  "delaunay2d_add: failed to update triangle 2\n");

      RJMCMC_CONDITIONCHECKINT(delete_edge(d, a1, b1), 
			       "delaunay2d_add: failed to remove edge bb = 0, 0\n");
      RJMCMC_CONDITIONCHECKINT(add_edge(d, a1, pi),
			       "delaunay2d_add: failed to add edge bb = 0, 0\n");

      /* 3rd triangle (new) */
      RJMCMC_CONDITIONCHECKINT(add_triangle(d, a2, pi, c2, t4, t2, n3) < 0, 
		  "delaunay2d_add: failed to add new triangle\n");
      if (n3 >= 0) {
	RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, n3, t2, t3) < 0,
		    "delaunay2d_add: failed to replace neighbour\n");
      }

      /* 4th triangle (new) */
      RJMCMC_CONDITIONCHECKINT(add_triangle(d, c1, pi, b1, t1, t3, n4) < 0, 
		  "delaunay2d_add: failed to add new triangle\n");
      if (n3 >= 0) {
	RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, n4, t1, t4) < 0,
		    "delaunay2d_add: failed to replace neighbour\n");
      }
      
      /* Flip edges as required */
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, ti, 2, 0) < 0,
		  "delaunay2d_add: failed to flip triangle\n");
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, t2, 2, 0) < 0,
		  "delaunay2d_add: failed to flip new triangle\n");
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, t3, 2, 0) < 0,
		  "delaunay2d_add: failed to flip new triangle\n");
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, t4, 2, 0) < 0,
		  "delaunay2d_add: failed to flip new triangle\n");
	
    }
  } else {

    /*
     * Point completely within triangle, ie not on an edge.
     */

    a = d->t[ti].v[0];
    b = d->t[ti].v[1];
    c = d->t[ti].v[2];

    bbox2d_expand(bound, d->p[a].x, d->p[a].y);
    bbox2d_expand(bound, d->p[b].x, d->p[b].y);
    bbox2d_expand(bound, d->p[c].x, d->p[c].y);
    
    /* nab = d->t[ti].n[0]; */
    nbc = d->t[ti].n[1];
    nca = d->t[ti].n[2];
    
    /* fprintf(stderr, "\nPre add\n"); */
    /* delaunay2d_print_points(d); */
    /* delaunay2d_print_triangles(d); */

    pi = d->np;
    if (add_point(d, x, y) < 0) {
      return -1;
    }
    
    tbc = d->nt;
    tca = tbc + 1;
    
    /* 
     *  Replace containing triangle (a, b, np) 
     */

    /* d->t[ti].v[0] = a; already correct */
    /* d->t[ti].v[1] = b; already correct */
    d->t[ti].v[2] = pi;
    
    /*d->t[ti].n[0] = nab; already correct */
    d->t[ti].n[1] = tbc;
    d->t[ti].n[2] = tca;

    RJMCMC_CONDITIONCHECKINT(triangle_update(d, ti) < 0,
		"delaunay2d_add: failed to update new triangle 1\n");
    
    /* Add 1st Triangle (b, c, np) */
    RJMCMC_CONDITIONCHECKINT(add_triangle(d, b, c, pi, nbc, tca, ti) < 0, 
		"delaunay2d_add: failed to add triangle b,c,np\n");

    if (nbc >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, nbc, ti, tbc) < 0, 
		  "delaunay2d_add: failed to update neighbour for nbc\n");
    }

    /* Add 2nd Triangle (c, a, np) */
    RJMCMC_CONDITIONCHECKINT(add_triangle(d, c, a, pi, nca, ti, tbc) < 0, 
		"delaunay2d_add: failed to add triangle c,a,np\n");

    if (nca >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, nca, ti, tca) < 0, 
		  "delaunay2d_add: failed to update neighbour for nca");
    }

    /* RJMCMC_CONDITIONCHECKINT(delaunay2d_validate_neighbours(d) < 0, */
    /* 		"delaunay2d_add: neighbours mismatched (pre flip)\n"); */

    /* fprintf(stderr, "\nPre flipping %d\n", ti); */
    /* delaunay2d_print_points(d); */
    /* delaunay2d_print_triangles(d); */
    
    /*
     * Flip edges as required 
     */
    RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, ti, 0, 0) < 0,
		"delaunay2d_add: failed to flip new triangle 1\n");
    RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, tbc, 0, 0) < 0,
		"delaunay2d_add: failed to flip new triangle 2\n");
    RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, tca, 0, 0) < 0,
		"delaunay2d_add: failed to flip new triangle 3\n");
  }

  return 0;
}

int delaunay2d_delete(delaunay2d_t *d,
		      int pi,
		      bbox2d_t *bound)
{
  int i;
  int j;
  int ti;
  int nd;

  int new_ti;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_delete: null delaunay\n");
  RJMCMC_CONDITIONCHECKINT(bound == NULL, "delaunay2d_delete: null bound\n");
  RJMCMC_CONDITIONCHECKINT(pi < 4, "delaunay2d_delete: cannot delete corner points\n");
  RJMCMC_CONDITIONCHECKINT(pi >= d->np, "delaunay2d_delete: invalid point index\n");

  /*
   * Find the first triangle that has the point as one of its vertices
   */
  ti = -1;
  for (i = 0; i < d->nt; i ++) {

    if ((d->t[i].v[0] == pi) ||
	(d->t[i].v[1] == pi) ||
	(d->t[i].v[2] == pi)) {
      ti = i;
      break;
    }
  }

  RJMCMC_CONDITIONCHECKINT(ti < 0, "delaunay2d_delete: unable to find first triangle\n");
  nd = fill_fan_about_point(d, pi, ti, bound);

  /*
   * Remove fan edges
   */
  for (j = 0; j < nd; j ++) {
    if (delete_edge(d, pi, d->vlist[j]) < 0) {
      fprintf(stderr, "delaunay2d_delete: failed to remove edge %d %d\n",
	      pi, d->vlist[j]);
    }
  }
			     

  RJMCMC_CONDITIONCHECKINT(clear_flip(d) < 0, 
	      "delaunay2d_delete: failed to clear flip list\n");
  
  RJMCMC_CONDITIONCHECKINT(delaunay2d_fill_hole(d, nd) < 0,
	      "delaunay2d_delete: failed to fill hole\n");

  nd = 2;
  for (j = 0; j < nd; j ++) {
    /*
     * Mark triangles as deleted temporarily
     */
    ti = d->tlist[j];
    d->t[ti].v[0] = -1;
    d->t[ti].v[1] = -1;
    d->t[ti].v[2] = -1;
  }

  for (j = 0; j < nd; j ++) {
    if (d->tlist[j] < d->nt) {
      RJMCMC_CONDITIONCHECKINT(delete_triangle(d, d->tlist[j], &new_ti) < 0,
		  "delaunay2d_delete: failed to remove triangle\n");

      
      if (new_ti >= 0) {
	RJMCMC_CONDITIONCHECKINT(update_flip(d, new_ti, d->tlist[j]) < 0,
		    "delaunay2d_delete: failed to update flip list\n");
      }
    }

    d->nt --;
  }

  /* remove point */
  RJMCMC_CONDITIONCHECKINT(delete_point(d, pi) < 0,
	      "delaunay2d_delete: failed to remove point\n");
  d->np --;

  /* Flip triangles as required */
  /* for (i = 0; i < d->nf; i ++) { */
  for (i = d->nf - 1; i >= 0; i --) {
    RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, d->ftlist[i], d->felist[i], -1) < 0,
		"delaunay2d_delete: failed to validate triangle pair\n");
  }
  d->nf = 0;

  return 0;
}

int 
delaunay2d_shift_replace(delaunay2d_t *d,
			 int pi)
{
  int i;
  int j;
  double px;
  double py;

  int n;
  int edge[MAXEDGES];

  /*
   * This functions swaps the last point with the index pi. It needs to update the triangles
   * and edges to maintain a consistent triangulation.
   */

  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_shift_replace: null delaunay\n");
  RJMCMC_CONDITIONCHECKINT(pi < 4, "delaunay2d_shift_replace: cannot delete corner points\n");
  RJMCMC_CONDITIONCHECKINT(pi >= d->np, "delaunay2d_shift_replace: invalid point index\n");

  px = d->p[d->np - 1].x;
  py = d->p[d->np - 1].y;
  n = d->e[d->np - 1].n;
  for (j = 0; j < n; j ++) {
    edge[j] = d->e[d->np - 1].edge[j];
  }

  for (i = (d->np - 1); i > pi; i --) {
    d->p[i].x = d->p[i - 1].x;
    d->p[i].y = d->p[i - 1].y; 

    d->e[i].n = d->e[i - 1].n;
    for (j = 0; j < d->e[i - 1].n; j ++) {
      d->e[i].edge[j] = d->e[i - 1].edge[j];
    }
  }

  d->p[pi].x = px;
  d->p[pi].y = py;

  d->e[pi].n = n;
  for (j = 0; j < n; j ++) {
    d->e[pi].edge[j] = edge[j];
  }

  /*
   * Update triangle indexing
   */
  for (i = 0; i < d->nt; i ++) {
    for (j = 0; j < 3; j ++) {
      if (d->t[i].v[j] == (d->np - 1)) {
	d->t[i].v[j] = pi;
      } else if (d->t[i].v[j] >= pi) {
	d->t[i].v[j] ++;
      }
    }
  }

  /*
   * Update edge indexing
   */
  for (i = 0; i < d->np; i ++) {
    for (j = 0; j < d->e[i].n; j ++) {
      if (d->e[i].edge[j] == (d->np - 1)) {
	d->e[i].edge[j] = pi;
      } else if (d->e[i].edge[j] >= pi) {
	d->e[i].edge[j] ++;
      }
    }
  }

  return 0;
}


int delaunay2d_find_enclosing_triangle(const delaunay2d_t *d,
				       int t0,
				       double px,
				       double py,
				       int *pa,
				       int *pb, 
				       int *pc,
				       double *ba,
				       double *bb,
				       double *bc)
{
  double _ba, _bb, _bc;
  int ti;
  int oldti;
  int i;
  int nsteps;
  
  ti = t0;
  nsteps = 0;
  while (nsteps < d->nt) {
    oldti = ti;
    i = barycentre_of_point(d,
			    ti,
			    px,
			    py,
			    &_ba,
			    &_bb,
			    &_bc);
    if (i < 0) {
      return -1;
    }

    if (_ba < -BARYCENTRE_EPSILON) {
      /* Point in direction of edge bc */
      ti = d->t[ti].n[1];
    } else if (_bb < -BARYCENTRE_EPSILON) {
      /* Point in direction of edge ca */
      ti = d->t[ti].n[2];
    } else if (_bc < -BARYCENTRE_EPSILON) {
      /* Point in direction of edge ab */
      ti = d->t[ti].n[0];
    } else {
      if ((_ba <= BARYCENTRE_ONE_EPSILON) &&
	  (_bb <= BARYCENTRE_ONE_EPSILON) && 
	  (_bc <= BARYCENTRE_ONE_EPSILON)) {
	/* Inside this triangle */
	
	*ba = _ba;
	*bb = _bb;
	*bc = _bc;
	
	*pa = d->t[ti].v[0];
	*pb = d->t[ti].v[1];
	*pc = d->t[ti].v[2];
	
	return ti;
      } else {
	/* Unknown error */
	fprintf(stderr, "delaunay2d_find_enclosing_triangle: invalid barycentre coordinates: %f %f %f\n", 
		_ba, _bb, _bc);
	return -1;
      }
    }

    if (ti < 0) {
      /* Point is outside our region (perhaps check for this on initial entry with min/max parameters */
      fprintf(stderr, "delaunay2d_find_enclosing_triangle: edge found %d (%g %g %g)\n", ti, _ba, _bb, _bc);
      fprintf(stderr, "                                  : %f %f (%f %f %f %f)\n", 
	      px, py, d->xmin, d->xmax, d->ymin, d->ymax);
      fprintf(stderr, "                                  : %d (%d %d %d) (%d %d %d)\n",
	      oldti, d->t[oldti].v[0], d->t[oldti].v[1], d->t[oldti].v[2],
	      d->t[oldti].n[0], d->t[oldti].n[1], d->t[oldti].n[2]);
      if (delaunay2d_validate_neighbours(d) < 0) {
	fprintf(stderr, "                                  : invalid neighbours\n");
      }
      if (delaunay2d_validate_circumcircles(d) < 0) {
	fprintf(stderr, "                                  : invalid circumcircles\n");
      }
      if (delaunay2d_validate_nonintersecting(d) < 0) {
	fprintf(stderr, "                                  : invalid non-intersecting\n");
      }
      return -1;
    }

    nsteps ++;
  }

    
  fprintf(stderr, "delaunay2d_find_enclosing_triangle: error finding triangle\n");
  delaunay2d_print_points(d);
  delaunay2d_print_triangles(d);
  delaunay2d_validate_neighbours(d);
  delaunay2d_validate_circumcircles(d);
  delaunay2d_validate_nonintersecting(d);
  return -1;
}

int 
delaunay2d_nearest(delaunay2d_t *d,
		   int include_corners,
		   double px,
		   double py)
{
  return delaunay2d_nearest_from(d, 0, include_corners, px, py);
}

int delaunay2d_nearest_from(delaunay2d_t *d,
			    int p0,
			    int include_corners,
			    double px,
			    double py)
{
  int i;
  int pi;
 
  double best_dist2;
  int best_pi;

  double d2;

  /*
   * New point so increment the cache index
   */
  d->ci ++;

  pi = p0;
  best_dist2 = dist2(d, p0, px, py);
  best_pi = -1;

  while (best_pi != pi) {

    for (i = 0; i < d->e[pi].n; i ++) {
      d2 = dist2(d, d->e[pi].edge[i], px, py);
      if (d2 < best_dist2) {
	best_dist2 = d2;
	best_pi = d->e[pi].edge[i];
      }
    }

    if (best_pi < 0) {
      best_pi = pi;
    } else {
      pi = best_pi;
      best_pi = -1;
    }
  }

  if (include_corners == 0 &&
      pi < 4) {
    /*
     * We can't use this point, recycle through neighbours to find the nearest.
     */
    best_pi = -1;
    best_dist2 = d->d2max;

    for (i = 0; i < d->e[pi].n; i ++) {

      if (d->e[pi].edge[i] > 3) {
	d2 = dist2(d, d->e[pi].edge[i], px, py);
	if (d2 < best_dist2) {
	  best_dist2 = d2;
	  best_pi = d->e[pi].edge[i];
	}
      }
    }
    
    if (best_pi < 0) {
      fprintf(stderr, "delaunay2d_nearest_from: unable to find nearest point\n");
      delaunay2d_validate_edges(d);
      delaunay2d_validate_delaunay(d);
      /* fprintf(stderr, "points\n"); */
      /* delaunay2d_print_points(d); */
      /* fprintf(stderr, "triangles\n"); */
      /* delaunay2d_print_triangles(d); */
    }
    pi = best_pi;
  }


  return pi;
}
      

int 
delaunay2d_point_of_index(const delaunay2d_t *d,
			  int i,
			  double *px,
			  double *py)
{
  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_point_of_index: NULL delaunay\n");
  RJMCMC_CONDITIONCHECKINT(i >= d->np, "delaunay2d_point_of_index: invalid index\n");
  
  *px = d->p[i].x;
  *py = d->p[i].y;

  return 0;
}

int delaunay2d_index_of_point(const delaunay2d_t *d,
			      double x,
			      double y)
{
  int i;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_point_of_index: NULL delaunay\n");

  for (i = 0; i < d->np; i ++) {
    if ((x == d->p[i].x) && (y == d->p[i].y)) {
      return i;
    }
  }

  return -1;
}

int
delaunay2d_npoints(const delaunay2d_t *d)
{
  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_npoints: NULL delaunay\n");

  return d->np;
}

int
delaunay2d_polygon_bound(const delaunay2d_t *d,
			 int pi,
			 bbox2d_t *bound)
{
  int i;
  int ti;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_polygon_bound: NULL delaunay\n");

  /*
   * Find the first triangle that has the point as one of its vertices
   */
  ti = -1;
  for (i = 0; i < d->nt; i ++) {

    if ((d->t[i].v[0] == pi) ||
	(d->t[i].v[1] == pi) ||
	(d->t[i].v[2] == pi)) {
      ti = i;
      break;
    }
  }

  RJMCMC_CONDITIONCHECKINT(ti < 0, "delaunay2d_polygon_bound: unable to find first triangle\n");
  
  return 0;
}

int
delaunay2d_ct_update(delaunay2d_t *d,
		     const double *point_z,
		     int n)
{
  int i;
  int j;
  double dx;
  double dy;
  double dz;

  double pdx;
  double pdy;

  double x1, y1, z1, z1x, z1y;
  double x2, y2, z2, z2x, z2y;
  double x3, y3, z3, z3x, z3y;

  double x4, y4;
  double theta1, theta2, theta3;

  if (d->np != n) {
    fprintf(stderr, "error: points mismatch %d != %d\n", d->np, n);
    return -1;
  }

  /* Go through points and set z */
  for (i = 0; i < n; i ++) {
    d->p[i].z = point_z[i];
  }

  /* Go through points and compute dx and dy */
  for (i = 0; i < n; i ++) {
    dx = 0.0;
    dy = 0.0;

    for (j = 0; j < d->e[i].n; j ++) {
      pdx = d->p[j].x - d->p[i].x;
      pdy = d->p[j].y - d->p[i].y;
      dz = d->p[j].z - d->p[i].z;

      if (pdx != 0.0) {
	dx += dz/pdx;
      }
      if (pdy != 0.0) {
	dy += dz/pdy;
      }
    }
  }
      
  /* Go through triangles and compute ct parameters */
  for (i = 0; i < d->nt; i ++) {
    
    x1 = d->p[d->t[i].v[0]].x;
    x2 = d->p[d->t[i].v[1]].x;
    x3 = d->p[d->t[i].v[2]].x;

    x4 = (x1 + x2 + x3)/3.0;
    y4 = (y1 + y2 + y3)/3.0;

    y1 = d->p[d->t[i].v[0]].y;
    y2 = d->p[d->t[i].v[1]].y;
    y3 = d->p[d->t[i].v[2]].y;

    z1 = d->p[d->t[i].v[0]].z;
    z2 = d->p[d->t[i].v[1]].z;
    z3 = d->p[d->t[i].v[2]].z;

    z1x = d->p[d->t[i].v[0]].dx;
    z2x = d->p[d->t[i].v[1]].dx;
    z3x = d->p[d->t[i].v[2]].dx;
    
    z1y = d->p[d->t[i].v[0]].dy;
    z2y = d->p[d->t[i].v[1]].dy;
    z3y = d->p[d->t[i].v[2]].dy;

    d->t[i].ct[1] = z1;
    d->t[i].ct[2] = z2;
    d->t[i].ct[3] = z3;

    d->t[i].ct[4] = ((x2 - x1) * z1x + (y2 - y1) * z1y)/3.0 + z1;
    d->t[i].ct[5] = ((x4 - x1) * z1x + (y4 - y1) * z1y)/3.0 + z1;
    d->t[i].ct[6] = ((x3 - x1) * z1x + (y3 - y1) * z1y)/3.0 + z1;

    d->t[i].ct[7] = ((x3 - x2) * z2x + (y3 - y2) * z2y)/3.0 + z2;
    d->t[i].ct[8] = ((x4 - x2) * z2x + (y4 - y2) * z2y)/3.0 + z2;
    d->t[i].ct[9] = ((x1 - x2) * z2x + (y1 - y2) * z2y)/3.0 + z2;

    d->t[i].ct[10] = ((x1 - x2) * z3x + (y1 - y2) * z3y)/3.0 + z3;
    d->t[i].ct[11] = ((x4 - x2) * z3x + (y4 - y2) * z3y)/3.0 + z3;
    d->t[i].ct[12] = ((x2 - x2) * z3x + (y2 - y2) * z3y)/3.0 + z3;

    pdx = x2 - x1;
    pdy = y2 - y1;
    theta1 = ((x4 - x1)*pdx + (y4 - y1)*pdy)/
      (pdx*pdx + pdy*pdy);

    pdx = x3 - x2;
    pdy = y3 - y2;
    theta2 = ((x4 - x2)*pdx + (y4 - y2)*pdy)/
      (pdx*pdx + pdy*pdy);

    pdx = x1 - x3;
    pdy = y1 - y3;
    theta3 = ((x4 - x3)*pdx + (y4 - y3)*pdy)/
      (pdx*pdx + pdy*pdy);
    

    d->t[i].ct[13] = (d->t[i].ct[5] + 
		      d->t[i].ct[8] + 
		      (theta1 - 1.0) *d->t[i].ct[1] + 
		      (2.0 - 3.0*theta1) * d->t[i].ct[4] + 
		      (3.0 * theta1 - 1.0) * d->t[i].ct[9] - 
		      theta1 * d->t[i].ct[2])/2.0;
    d->t[i].ct[14] = (d->t[i].ct[8] + 
		      d->t[i].ct[11] + 
		      (theta2 - 1.0) *d->t[i].ct[2] + 
		      (2.0 - 3.0*theta2) * d->t[i].ct[7] + 
		      (3.0 * theta2 - 1.0) * d->t[i].ct[12] - 
		      theta2 * d->t[i].ct[3])/2.0;
    d->t[i].ct[15] = (d->t[i].ct[11] + 
		      d->t[i].ct[5] + 
		      (theta3 - 1.0) *d->t[i].ct[3] + 
		      (2.0 - 3.0*theta3) * d->t[i].ct[10] + 
		      (3.0 * theta3 - 1.0) * d->t[i].ct[6] - 
		      theta3 * d->t[i].ct[1])/2.0;

    d->t[i].ct[16] = (d->t[i].ct[15] + 
		      d->t[i].ct[5] + 
		      d->t[i].ct[13])/3.0;
    d->t[i].ct[17] = (d->t[i].ct[13] + 
		      d->t[i].ct[8] + 
		      d->t[i].ct[14])/3.0;
    d->t[i].ct[18] = (d->t[i].ct[14] + 
		      d->t[i].ct[11] + 
		      d->t[i].ct[15])/3.0;

    d->t[i].ct[19] = (d->t[i].ct[18] + 
		      d->t[i].ct[16] + 
		      d->t[i].ct[17])/3.0;
  }

  return 0;
    
}

int 
delaunay2d_ct_value_at(delaunay2d_t *d,
		       int t0,
		       double px,
		       double py,
		       double *z)
{
   
  int t;
  int pa, pb, pc;
  double ba, bb, bc;

  int subtriangle;

  double x1, y1;
  double x2, y2;
  double x3, y3;
  double x4, y4;

  double detT;
  double dx, dy;
  double u, v, w;

  t = delaunay2d_find_enclosing_triangle(d,
					 t0,
					 px,
					 py,
					 &pa,
					 &pb, 
					 &pc,
					 &ba,
					 &bb,
					 &bc);

  subtriangle = -1;
  if (ba > bb) {
    if (ba > bc) {
      subtriangle = 3;
    } else {
      subtriangle = 1;
    }
  } else {
    if (ba > bc) {
      subtriangle = 1;
    } else {
      subtriangle = 2;
    }
  }

  x1 = d->p[d->t[t].v[0]].x;
  y1 = d->p[d->t[t].v[0]].y;
  
  x2 = d->p[d->t[t].v[1]].x;
  y2 = d->p[d->t[t].v[1]].y;
  
  x3 = d->p[d->t[t].v[2]].x;
  y3 = d->p[d->t[t].v[2]].y;

  x4 = (x1 + x2 + x3)/3.0;
  y4 = (y1 + y2 + y3)/3.0;

  switch(subtriangle) {
  case 1:
    /* 1 2 4 */
    detT = (x1 - x4) * (y2 - y4) - (x2 - x4) * (y1 - y4);

    dx = (px - x4);
    dy = (py - y4);

    u = ((y2 - y4)*dx + (x4 - x2)*dy)/detT;
    v = ((y4 - y1)*dx + (x1 - x4)*dy)/detT;
    w = 1.0 - (u + v);

    *z = 
      d->t[t].ct[19] *w*w*w + 
      3.0 * d->t[t].ct[16] *w*w*u + 
      3.0 * d->t[t].ct[17] *w*w*v +
      3.0 * d->t[t].ct[5] *w*u*u + 
      6.0 * d->t[t].ct[13] *w*u*v + 
      3.0 * d->t[t].ct[8] * w*v*v + 
      d->t[t].ct[1] * u*u*u +
      3.0 * d->t[t].ct[4] * u*u*v + 
      3.0 * d->t[t].ct[9] * u*v*v + 
      d->t[t].ct[2] * v*v*v;
    break;

  case 2:
    /* 2 3 4 */
    detT = (x2 - x4) * (y3 - y4) - (x3 - x4) * (y2 - y4);

    dx = (px - x4);
    dy = (py - y4);

    u = ((y3 - y4)*dx + (x4 - x3)*dy)/detT;
    v = ((y4 - y2)*dx + (x2 - x4)*dy)/detT;
    w = 1.0 - (u + v);

    *z = 
      d->t[t].ct[19] *w*w*w + 
      3.0 * d->t[t].ct[17] *w*w*u + 
      3.0 * d->t[t].ct[18] *w*w*v +
      3.0 * d->t[t].ct[8] *w*u*u + 
      6.0 * d->t[t].ct[14] *w*u*v + 
      3.0 * d->t[t].ct[11] * w*v*v + 
      d->t[t].ct[2] * u*u*u +
      3.0 * d->t[t].ct[7] * u*u*v + 
      3.0 * d->t[t].ct[12] * u*v*v + 
      d->t[t].ct[3] * v*v*v;
    break;

  case 3:
    /* 3 1 4 */
    detT = (x3 - x4) * (y1 - y4) - (x1 - x4) * (y3 - y4);

    dx = (px - x4);
    dy = (py - y4);

    u = ((y1 - y4)*dx + (x4 - x1)*dy)/detT;
    v = ((y4 - y3)*dx + (x3 - x4)*dy)/detT;
    w = 1.0 - (u + v);

    *z = 
      d->t[t].ct[19] *w*w*w + 
      3.0 * d->t[t].ct[18] *w*w*u + 
      3.0 * d->t[t].ct[16] *w*w*v +
      3.0 * d->t[t].ct[11] *w*u*u + 
      6.0 * d->t[t].ct[15] *w*u*v + 
      3.0 * d->t[t].ct[5] * w*v*v + 
      d->t[t].ct[3] * u*u*u +
      3.0 * d->t[t].ct[10] * u*u*v + 
      3.0 * d->t[t].ct[6] * u*v*v + 
      d->t[t].ct[1] * v*v*v;
    break;
    
  default:
    fprintf(stderr, "error: invalid subtriangle %d\n", subtriangle);
    return -1;
  }

  return 0;
  
}



static int add_point(delaunay2d_t *d,
		     double x,
		     double y)
{
  RJMCMC_CONDITIONCHECKINT(d == NULL, "add_point: NULL delaunay\n");
  RJMCMC_CONDITIONCHECKINT(d->np >= d->maxpoints, "add_point: max points reached\n");
 
  /*
   * Add actual point
   */
  d->p[d->np].x = x;
  d->p[d->np].y = y;

  /*
   * Reset edge list
   */
  d->e[d->np].n = 0;

  d->np ++;

  return 0;
}

static int add_edge(delaunay2d_t *d,
		    int p,
		    int np)
{
  int i;
  int existing;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "add_edge: NULL delaunay\n");
  RJMCMC_CONDITIONCHECKINT(p >= d->np, "add_edge: invalid point\n");
  RJMCMC_CONDITIONCHECKINT(np >= d->np, "add_edge: invalidate neighbour point\n");

  /*
   * Add n -> np 
   */
  existing = 0;
  for (i = 0; i < d->e[p].n; i ++) {
    if (d->e[p].edge[i] == np) {
      existing = -1;
      break;
    }
  }
  if (!existing) {
    d->e[p].edge[d->e[p].n] = np;
    d->e[p].n ++;
  }

  /*
   * Add np -> n
   */
  existing = 0;
  for (i = 0; i < d->e[np].n; i ++) {
    if (d->e[np].edge[i] == p) {
      existing = -1;
      break;
    }
  }
  if (!existing) {
    d->e[np].edge[d->e[np].n] = p;
    d->e[np].n ++;
  }

  return 0;
}

static int add_triangle(delaunay2d_t *d,
			int a, 
			int b,
			int c,
			int nab, 
			int nbc, 
			int nca)
{
  int new_st;
  triangle_t *new_t;
  int i;
  double detT;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "add_triangle: NULL delaunay\n");
  RJMCMC_CONDITIONCHECKINT(a >= d->np, "add_triangle: invalid a index\n");
  RJMCMC_CONDITIONCHECKINT(b >= d->np, "add_triangle: invalid b index\n");
  RJMCMC_CONDITIONCHECKINT(c >= d->np, "add_triangle: invalid c index\n");

  detT = compute_detT(d, a, b, c);
  if (detT == 0.0) {
    rjmcmc_error("add_triangle: colinear points: (%f %f) (%f %f) (%f %f)\n",
		 d->p[a].x, d->p[a].y,
		 d->p[b].x, d->p[b].y,
		 d->p[c].x, d->p[c].y);
    return -1;
  }

  if (d->nt == d->st) {
    new_st = d->st + TRIANGLE_INCR;
    new_t = malloc(sizeof(triangle_t) * new_st);
    RJMCMC_NULLCHECKINT(new_t, "add_triangle: failed to resize triangles\n");
    
    for (i = 0; i < d->nt; i ++) {
      triangle_copy(&(d->t[i]), &(new_t[i]));
    }

    free(d->t);
    d->t = new_t;
    d->st = new_st;
  }


  i = d->nt;
  d->t[i].v[0] = a;
  d->t[i].v[1] = b;
  d->t[i].v[2] = c;
    
  d->t[i].n[0] = nab;
  d->t[i].n[1] = nbc;
  d->t[i].n[2] = nca;

  d->t[i].detT = detT;

  RJMCMC_CONDITIONCHECKINT(circumcircle_of_triangle(d,
				       i,
				       &(d->t[i].cx),
				       &(d->t[i].cy),
				       &(d->t[i].cr2)) < 0,
	      "add_triangle: failed to compute circumcircle\n");

  /*
   * Add edges to the edge list
   */
  RJMCMC_CONDITIONCHECKINT(add_edge(d, a, b) < 0,
			   "add_triangle: failed to add edge ab\n");
  RJMCMC_CONDITIONCHECKINT(add_edge(d, b, c) < 0,
			   "add_triangle: failed to add edge bc\n");
  RJMCMC_CONDITIONCHECKINT(add_edge(d, c, a) < 0,
			   "add_triangle: failed to add edge ca\n");
    
  d->nt ++;

  return 0;
}

static void triangle_copy(const triangle_t *src,
			  triangle_t *dst)
{
  int i;

  for (i = 0; i < 3; i ++) {
    dst->v[i] = src->v[i];
    dst->n[i] = src->n[i];
  }

  dst->detT = src->detT;

  dst->cx = src->cx;
  dst->cy = src->cy;
  dst->cr2 = src->cr2;
}

static int triangle_replace_neighbour(delaunay2d_t *d,
				      int ti,
				      int oldn,
				      int newn)
{
  int i;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "triangle_replace_neighbour: null delaunay\n");
  RJMCMC_CONDITIONCHECKINT(ti >= d->nt, "triangle_replace_neighbour: triangle out of range\n");
  RJMCMC_CONDITIONCHECKINT(oldn >= d->nt, "triangle_replace_neighbour: old neighbour out of range\n");
  RJMCMC_CONDITIONCHECKINT(newn >= d->nt, "triangle_replace_neighbour: new neighbour out of range\n");

  for (i = 0; i < 3; i ++) {
    if (d->t[ti].n[i] == oldn) {
      d->t[ti].n[i] = newn;
      return 0;
    }
  }

  fprintf(stderr, "triangle_replace_neighbour: no neighbour %d in triangle %d (%d %d %d)\n",
	  oldn,
	  ti,
	  d->t[ti].n[0],
	  d->t[ti].n[1],
	  d->t[ti].n[2]);
  return -1;
}

static int triangle_verify_neighbour(delaunay2d_t *d,
				     int ti,
				     int ei)
{
  int nti;
  int p1;
  int p2;

  RJMCMC_CONDITIONCHECKINT(((ei < 0) || (ei > 2)), "triangle_verify_neighbour: invalid edge index\n");
  nti = d->t[ti].n[ei];
  p1 = d->t[ti].v[ei];
  p2 = d->t[ti].v[(ei + 1) % 3];

  if (d->t[nti].v[0] == p2) {
    if (d->t[nti].v[1] != p1) {
      fprintf(stderr, "triangle_verify_neighbour: edge mismatch (%d %d) != (%d %d)\n", p1, p2, 
	      d->t[nti].v[1], d->t[nti].v[0]);
      return -1;
    }
    d->t[nti].n[0] = ti;
  } else if (d->t[nti].v[1] == p2) {
    if (d->t[nti].v[2] != p1) {
      fprintf(stderr, "triangle_verify_neighbour: edge mismatch (%d %d) != (%d %d)\n", p1, p2, 
	      d->t[nti].v[2], d->t[nti].v[1]);
      return -1;
    }
    d->t[nti].n[1] = ti;
    
  } else if (d->t[nti].v[2] == p2) {
    if (d->t[nti].v[0] != p1) {
      fprintf(stderr, "triangle_verify_neighbour: edge mismatch (%d %d) != (%d %d)\n", p1, p2, 
	      d->t[nti].v[0], d->t[nti].v[2]);
      return -1;
    }
    d->t[nti].n[2] = ti;

  } else {

    fprintf(stderr, "triangle_verify_neighbour: point %d not in triangle %d\n", p2, nti);
    return -1;
  }
  

  return 0;
}


static int triangle_update(delaunay2d_t *d,
			   int ti)
{
  RJMCMC_CONDITIONCHECKINT(d == NULL, "triangle_update: NULL delaunay\n");
  RJMCMC_CONDITIONCHECKINT(ti >= d->nt, "triangle_update: index out of range\n");

  RJMCMC_CONDITIONCHECKINT(circumcircle_of_triangle(d, 
				       ti,
				       &(d->t[ti].cx),
				       &(d->t[ti].cy),
				       &(d->t[ti].cr2)) < 0, 
	      "triangle_update: failed to compute circumcircle\n");

  d->t[ti].detT = compute_detT(d, 
			       d->t[ti].v[0],
			       d->t[ti].v[1],
			       d->t[ti].v[2]);
  RJMCMC_CONDITIONCHECKINT(d->t[ti].detT == 0.0,
	      "triangle_update: detT == 0.0\n");

  return 0;
}

static int validate_triangle_pair(delaunay2d_t *d,
				  int ti,
				  int edge, 
				  int bidirectional)
{
  /*
   * This function does the edge flipping required to maintain the Delaunay
   * criterion.
   */

  int n;
  int fp;
  int oedge;

  double dx;
  double dy;
  double r2;

  int a1, b1, c1;
  int nab1, nbc1, nca1;

  int a2, b2, c2;
  int nab2, nbc2, nca2;

  int n1;
  int n2;
  int f1;
  int f2;

  int e1a, e1b;
  int e2a, e2b;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "validate_triangle_pair: NULL delaunay\n");
  RJMCMC_CONDITIONCHECKINT(ti >= d->nt, "validate_triangle_pair: index out of range\n");
  
  n = d->t[ti].n[edge];
  if (n >= 0) {
    
    if (d->t[n].n[0] == ti) {
      fp = d->t[n].v[2];
      oedge = 0;
    } else if (d->t[n].n[1] == ti) {
      fp = d->t[n].v[0];
      oedge = 1;
    } else if (d->t[n].n[2] == ti) {
      fp = d->t[n].v[1];
      oedge = 2;
    } else {
      fprintf(stderr, "validate_triangle_pair: failed to find neighbour\n");
      fprintf(stderr, "                      : %d %d: %d %d %d\n",
	      ti, n, d->t[n].n[0], d->t[n].n[1], d->t[n].n[2]);
      return -1;
    }
    
    dx = d->t[ti].cx - d->p[fp].x;
    dy = d->t[ti].cy - d->p[fp].y;
    
    r2 = dx*dx + dy*dy;
    
    if (r2 < d->t[ti].cr2) {
      a1 = d->t[ti].v[edge];
      b1 = d->t[n].v[(oedge + 2) % 3];
      c1 = d->t[ti].v[(edge + 2) % 3];
      
      nab1 = d->t[n].n[(oedge + 1) % 3];
      nbc1 = n;
      nca1 = d->t[ti].n[(edge + 2) % 3];
      
      n1 = d->t[ti].n[(edge + 1) % 3];
      f1 = 0;
      
      a2 = d->t[n].v[oedge];
      b2 = d->t[ti].v[(edge + 2) % 3];
      c2 = d->t[n].v[(oedge + 2) % 3];
	
      nab2 = d->t[ti].n[(edge + 1) % 3];
      nbc2 = ti;
      nca2 = d->t[n].n[(oedge + 2) % 3];
      
      n2 = d->t[n].n[(oedge + 1) % 3];
      f2 = 2;
      
      e1a = d->t[ti].v[edge];
      e1b = d->t[n].v[oedge];
      e2a = d->t[ti].v[(edge + 2) % 3];
      e2b = d->t[n].v[(oedge + 2) % 3];
#if defined(DEBUG_FLIPS)
      fprintf(stderr, "flip: 0 0 %d (%d %d %d) %d (%d %d %d)\n", ti, a1, b1, c1, n, a2, b2, c2);
#endif
      
      /*
       * Update flipped triangles
       */
      d->t[ti].v[0] = a1;
      d->t[ti].v[1] = b1;
      d->t[ti].v[2] = c1;
      
      d->t[ti].n[0] = nab1;
      d->t[ti].n[1] = nbc1;
      d->t[ti].n[2] = nca1;
      
      RJMCMC_CONDITIONCHECKINT(triangle_update(d, ti) < 0,
			       "validate_triangle_pair: failed to update flipped triangle\n");
      
      if (n1 >= 0) {
	RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, n1, ti, n) < 0,
				 "validate_triangle_pair: failed to update flipped triangle 1's neighbour (0)\n");
      }
      
      d->t[n].v[0] = a2;
      d->t[n].v[1] = b2;
      d->t[n].v[2] = c2;
      
      d->t[n].n[0] = nab2;
      d->t[n].n[1] = nbc2;
      d->t[n].n[2] = nca2;
      
      RJMCMC_CONDITIONCHECKINT(triangle_update(d, n) < 0,
			       "validate_triangle_pair: failed to update flipped triangle\n");
      
      if (n2 >= 0) {
	RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, n2, n, ti) < 0,
				 "validate_triangle_pair: failed to update flipped triangle 2's neighbour (0)\n");
      }
      
      /* Update flipped edges */
      RJMCMC_CONDITIONCHECKINT(flip_edge(d, e1a, e1b, e2a, e2b) < 0,
			       "validate_triangle_pair: failed to update edges\n");
      
      /* Recursively flip */
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, ti, f1, bidirectional) < 0,
			       "validate_triangle_pair: recursive 1 failed\n");
      RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, n, f2, bidirectional) < 0,
			       "validate_triangle_pair: recursive 2 failed\n");

      if (bidirectional) {
	RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, ti, (f1 + 2) % 3, bidirectional) < 0,
				 "validate_triangle_pair: recursive 1 failed bidirectional\n");
	RJMCMC_CONDITIONCHECKINT(validate_triangle_pair(d, n, (f2 + 1) % 3, bidirectional) < 0,
				 "validate_triangle_pair: recursive 2 failed bidirectional\n");
      }
      
    }
  }
  
  return 0;
}

static double compute_detT(const delaunay2d_t *d,
			   int a,
			   int b,
			   int c)
{
  return 
    ((d->p[a].x - d->p[c].x) * (d->p[b].y - d->p[c].y)) -
    ((d->p[b].x - d->p[c].x) * (d->p[a].y - d->p[c].y));
}

static int barycentre_of_point(const delaunay2d_t *d,
			       int ti,
			       double x,
			       double y,
			       double *ba,
			       double *bb,
			       double *bc)
{
  int a, b, c;

  double dx3;
  double dy3;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "barycentre_of_point: NULL delaunay\n");
  RJMCMC_CONDITIONCHECKINT(ti >= d->nt, "barycentre_of_point: invalid triangle index\n");

  a = d->t[ti].v[0];
  b = d->t[ti].v[1];
  c = d->t[ti].v[2];

  dx3 = (x - d->p[c].x);
  dy3 = (y - d->p[c].y);

  *ba = (((d->p[b].y - d->p[c].y) * dx3) + ((d->p[c].x - d->p[b].x)*dy3))/d->t[ti].detT;
  *bb = (((d->p[c].y - d->p[a].y) * dx3) + ((d->p[a].x - d->p[c].x)*dy3))/d->t[ti].detT;
  *bc = 1.0 - (*ba) - (*bb);

  return 0;
}

int
barycentre(double px, double py,
	   double x1, double y1,
	   double x2, double y2,
	   double x3, double y3,
	   double *ba, 
	   double *bb,
	   double *bc)
{
  double detT;
  double dx3;
  double dy3;

  detT = ((x1 - x3) * (y2 - y3) -
	  ((x2 - x3) * (y1 - y3)));
  if (detT == 0.0) {
    return -1;
  }
	  
  dx3 = (px - x3);
  dy3 = (py - y3);

  *ba = ((y2 - y3)*dx3 + (x3 - x2)*dy3)/detT;
  *bb = ((y3 - y1)*dx3 + (x1 - x3)*dy3)/detT;
  *bc = 1.0 - (*ba) - (*bb);
  return 0;
}
   

int point_in_triangle(double px,
		      double py,
		      double x1, double y1,
		      double x2, double y2,
		      double x3, double y3)
{
  double ba, bb, bc;

  if (barycentre(px, py, x1, y1, x2, y2, x3, y3, &ba, &bb, &bc) < 0) {
    return 0;
  }

  if (ba < 0.0 || ba > 1.0) {
    return 0;
  }
  if (bb < 0.0 || bb > 1.0) {
    return 0;
  }
  if (bc < 0.0 || bc > 1.0) {
    return 0;
  }

  return -1;
}


static int circumcircle_of_triangle(const delaunay2d_t *d,
				    int ti,
				    double *cx,
				    double *cy,
				    double *r2)
{
  RJMCMC_CONDITIONCHECKINT(d == NULL, "circumcircle_of_triangle: NULL delaunay\n");

  return triangle_circumcircle(d->p[d->t[ti].v[0]].x, d->p[d->t[ti].v[0]].y,
			       d->p[d->t[ti].v[1]].x, d->p[d->t[ti].v[1]].y,
			       d->p[d->t[ti].v[2]].x, d->p[d->t[ti].v[2]].y,
			       cx, cy, r2);
}

static int neighbour_index(const delaunay2d_t *d,
			   int neighbour,
			   int tri)
{
  int i;

  /*
   * This is a bit counter-intuitive. This function returns the edge index of triangle
   * neighbour in triangle tri.
   */
  if (neighbour < 0 || tri < 0) {
    return -1;
  }

  for (i = 0; i < 3; i ++) {
    if (d->t[tri].n[i] == neighbour) {
      return i;
    }
  }

  return -1;
}

int 
delaunay2d_validate_circumcircles(const delaunay2d_t *d)
{
  int i;
  int j;
  double dx;
  double dy;
  double r2;
  
  double err;

  int error_count;
  
  RJMCMC_CONDITIONCHECKINT(d == NULL, "delaunay2d_validate_circumcircles: NULL delaunay\n");

  error_count = 0;

  for (i = 0; i < d->nt; i ++) {
    
    for (j = 0; j < 3; j ++) {
      dx = d->p[d->t[i].v[j]].x - d->t[i].cx;
      dy = d->p[d->t[i].v[j]].y - d->t[i].cy;
      
      r2 = dx*dx + dy*dy;
    
      err = fabs(r2 - d->t[i].cr2);
      if (err/d->t[i].cr2 > 1.0e-6) {
	fprintf(stderr, 
	      "delaunay2d_validate_circumcircles: triangle %d: %d: %g %g %g\n", 
		i, j, err, r2, d->t[i].cr2);
	fprintf(stderr,
		"                                 : %f %f -> %f %f\n", 
		d->p[d->t[i].v[0]].x, d->p[d->t[i].v[0]].y,
		d->t[i].cx, d->t[i].cy);
	error_count ++;
      }
    }
  }

  if (error_count > 0) {
    delaunay2d_print_points(d);
    delaunay2d_print_triangles(d);
    return -1;
  }

  return 0;
}

int
delaunay2d_validate_edges(const delaunay2d_t *d)
{
  int i;
  int j;
  int k;
  int p;
  int np;

  int found;
  int errorcount;

  errorcount = 0;

  /*
   * Edges mirrored
   */
  for (i = 0; i < d->np; i ++) {
    if (d->e[i].n == 0) {
      fprintf(stderr, "delaunay2d_validate_edges: point %d has no edges\n", i);
      errorcount ++;
    }

    for (j = 0; j < d->e[i].n; j ++) {
      p = d->e[i].edge[j];
      
      if (p >= d->np) {
	fprintf(stderr, "delaunay2d_validate_edges: point %d has invalid neighbour %d\n", i, p);
	errorcount ++;
      } else {
	found = 0;
	for (k = 0; k < d->e[p].n; k ++) {
	  if (d->e[p].edge[k] == i) {
	    found = -1;
	    break;
	  }
	}
	if (found == 0) {
	  fprintf(stderr, "delaunay2d_validate_edges: point %d has neighbour %d unmirrored\n", i, p);
	  errorcount ++;
	}
      }
    }
  }

  /*
   * Edges not duplicated
   */
  for (i = 0; i < d->np; i ++) {
    for (j = 0; j < (d->e[i].n - 1); j ++) {
      for (k = j + 1; k < d->e[i].n; k ++) {
	if (d->e[i].edge[j] == d->e[i].edge[k]) {
	  fprintf(stderr, "delaunay2d_validate_edges: point %d has duplicate neighbour %d\n",
		  i, d->e[i].edge[j]);
	  errorcount ++;
	}
      }
    }
  }

  /*
   * Edges for each triangle correct
   */
  for (i = 0; i < d->nt; i ++) {
    found = 0;
    p = d->t[i].v[0];
    np = d->t[i].v[1];
    for (j = 0; j < d->e[p].n; j ++) {
      if (d->e[p].edge[j] == np) {
	found = -1;
	break;
      }
    }
    if (!found) {
      fprintf(stderr, "delaunay2d_validate_edges: triangle %d ab edge %d -> %d missing\n", i, p, np);
      errorcount ++;
    }

    found = 0;
    p = d->t[i].v[1];
    np = d->t[i].v[2];
    for (j = 0; j < d->e[p].n; j ++) {
      if (d->e[p].edge[j] == np) {
	found = -1;
	break;
      }
    }
    if (!found) {
      fprintf(stderr, "delaunay2d_validate_edges: triangle %d bc edge %d -> %d missing\n", i, p, np);
      errorcount ++;
    }

    found = 0;
    p = d->t[i].v[2];
    np = d->t[i].v[0];
    for (j = 0; j < d->e[p].n; j ++) {
      if (d->e[p].edge[j] == np) {
	found = -1;
	break;
      }
    }
    if (!found) {
      fprintf(stderr, "delaunay2d_validate_edges: triangle %d ca edge %d -> %d missing\n", i, p, np);
      errorcount ++;
    }
  }

  /*
   * Every edge is part of a triangle
   */
  for (i = 0; i < d->np; i ++) {
    for (j = 0; j < d->e[i].n; j ++) {

      p = i;
      np = d->e[i].edge[j];

      found = 0;
      for (k = 0; k < d->nt; k ++) {
	if ((d->t[k].v[0] == p && d->t[k].v[1] == np) ||
	    (d->t[k].v[1] == p && d->t[k].v[2] == np) ||
	    (d->t[k].v[2] == p && d->t[k].v[0] == np) ||
	    (d->t[k].v[0] == np && d->t[k].v[1] == p) ||
	    (d->t[k].v[1] == np && d->t[k].v[2] == p) ||
	    (d->t[k].v[2] == np && d->t[k].v[0] == p)) {
	  found = -1;
	  break;
	}
      }

      if (!found) {
	fprintf(stderr, "delaunay2d_validate_edges: edge %d %d not in any triangle\n", p, np);
	errorcount ++;
      }
    }
  }

  if (errorcount > 0) {
    return -1;
  }

  return 0;
}

int
delaunay2d_validate_delaunay(const delaunay2d_t *d)
{
  int i;
  int j;

  double dx;
  double dy;
  double r2;
  
  int error_count;
  
  error_count = 0;

  for (i = 0; i < d->nt; i ++) {

    for (j = 0; j < d->np; j ++) {
      if ((j == d->t[i].v[0]) ||
	  (j == d->t[i].v[1]) ||
	  (j == d->t[i].v[2])) {
	continue;
      }


      dx = d->p[j].x - d->t[i].cx;
      dy = d->p[j].y - d->t[i].cy;

      r2 = dx*dx + dy*dy;
    
      if ((r2 + 1.0e-6) < d->t[i].cr2) {
	fprintf(stderr, "point %d (%f, %f) is in cc of triangle %d: (%f, %f) %f %f %.10f\n", 
		j, d->p[j].x, d->p[j].y,
		i, d->t[i].cx, d->t[i].cy, d->t[i].cr2, r2, fabs(d->t[i].cr2 - r2));
	error_count ++;
      }
    }
  }

  if (error_count > 0) {
    return -1;
  }

  return 0;
}

int
delaunay2d_validate_neighbours(const delaunay2d_t *d)
{
  int i;
  int error_count;
  int ni;

  error_count = 0;
  for (i = 0; i < d->nt; i ++) {

    /*
     * Neighbour AB
     */
    if (d->t[i].n[0] >= 0) {
      ni = neighbour_index(d, i, d->t[i].n[0]);
      switch(ni) {
      case 0:
	if (d->t[i].v[0] != d->t[d->t[i].n[0]].v[1] ||
	    d->t[i].v[1] != d->t[d->t[i].n[0]].v[0]) {
	  fprintf(stderr, "neighbour mismatch: %d %d != %d %d\n",
		  d->t[i].v[0], d->t[i].v[1],
		  d->t[d->t[i].n[0]].v[1], d->t[d->t[i].n[0]].v[0]);
	  error_count ++;
	}
	break;

      case 1:
	if (d->t[i].v[0] != d->t[d->t[i].n[0]].v[2] ||
	    d->t[i].v[1] != d->t[d->t[i].n[0]].v[1]) {
	  fprintf(stderr, "neighbour mismatch: %d %d != %d %d\n",
		  d->t[i].v[0], d->t[i].v[1],
		  d->t[d->t[i].n[0]].v[2], d->t[d->t[i].n[0]].v[1]);
	  error_count ++;
	}
	break;

      case 2:
	if (d->t[i].v[0] != d->t[d->t[i].n[0]].v[0] ||
	    d->t[i].v[1] != d->t[d->t[i].n[0]].v[2]) {
	  fprintf(stderr, "neighbour mismatch: %d %d != %d %d\n",
		  d->t[i].v[0], d->t[i].v[1],
		  d->t[d->t[i].n[0]].v[0], d->t[d->t[i].n[0]].v[2]);
	  error_count ++;
	}
	break;

      default:
	fprintf(stderr, "neighbour mismatch: %d edge 0 %d\n", i, d->t[i].n[0]);
	error_count ++;
	break;
      }
    }

    /*
     * Neighbour BC
     */
    if (d->t[i].n[1] >= 0) {
      ni = neighbour_index(d, i, d->t[i].n[1]);
      switch(ni) {
      case 0:
	if (d->t[i].v[1] != d->t[d->t[i].n[1]].v[1] ||
	    d->t[i].v[2] != d->t[d->t[i].n[1]].v[0]) {
	  fprintf(stderr, "neighbour mismatch: %d %d != %d %d\n",
		  d->t[i].v[1], d->t[i].v[2],
		  d->t[d->t[i].n[1]].v[1], d->t[d->t[i].n[1]].v[0]);
	  error_count ++;
	}
	break;

      case 1:
	if (d->t[i].v[1] != d->t[d->t[i].n[1]].v[2] ||
	    d->t[i].v[2] != d->t[d->t[i].n[1]].v[1]) {
	  fprintf(stderr, "neighbour mismatch: %d %d != %d %d\n",
		  d->t[i].v[1], d->t[i].v[2],
		  d->t[d->t[i].n[1]].v[2], d->t[d->t[i].n[1]].v[1]);
	  error_count ++;
	}
	break;

      case 2:
	if (d->t[i].v[1] != d->t[d->t[i].n[1]].v[0] ||
	    d->t[i].v[2] != d->t[d->t[i].n[1]].v[2]) {
	  fprintf(stderr, "neighbour mismatch: %d %d != %d %d\n",
		  d->t[i].v[1], d->t[i].v[2],
		  d->t[d->t[i].n[1]].v[0], d->t[d->t[i].n[1]].v[2]);
	  error_count ++;
	}
	break;

      default:
	fprintf(stderr, "neighbour mismatch: %d edge 0 %d\n", i, d->t[i].n[1]);
	error_count ++;
	break;
      }
    }

    if (d->t[i].n[2] >= 0) {
      ni = neighbour_index(d, i, d->t[i].n[2]);
      switch(ni) {
      case 0:
	if (d->t[i].v[2] != d->t[d->t[i].n[2]].v[1] ||
	    d->t[i].v[0] != d->t[d->t[i].n[2]].v[0]) {
	  fprintf(stderr, "neighbour mismatch: triangle %d edge 2: %d %d != %d %d\n",
		  i,
		  d->t[i].v[2], d->t[i].v[0],
		  d->t[d->t[i].n[2]].v[1], d->t[d->t[i].n[2]].v[0]);
	  error_count ++;
	}
	break;

      case 1:
	if (d->t[i].v[2] != d->t[d->t[i].n[2]].v[2] ||
	    d->t[i].v[0] != d->t[d->t[i].n[2]].v[1]) {
	  fprintf(stderr, "neighbour mismatch: triangle %d edge 2: %d %d != %d %d\n",
		  i,
		  d->t[i].v[2], d->t[i].v[0],
		  d->t[d->t[i].n[2]].v[2], d->t[d->t[i].n[2]].v[1]);
	  error_count ++;
	}
	break;

      case 2:
	if (d->t[i].v[2] != d->t[d->t[i].n[2]].v[0] ||
	    d->t[i].v[0] != d->t[d->t[i].n[2]].v[2]) {
	  fprintf(stderr, "neighbour mismatch: triangle %d edge 2: %d %d != %d %d\n",
		  i,
		  d->t[i].v[2], d->t[i].v[0],
		  d->t[d->t[i].n[2]].v[0], d->t[d->t[i].n[2]].v[2]);
	  error_count ++;
	}
	break;

      default:
	fprintf(stderr, "neighbour mismatch: %d edge 0 %d\n", i, d->t[i].n[2]);
	error_count ++;
	break;
      }
    }
  }

  if (error_count > 0) {
    fprintf(stderr, "\n");
    delaunay2d_print_points(d);
    delaunay2d_print_triangles(d);
    return -1;
  }

  return 0;
    
}

int
delaunay2d_validate_nonintersecting(const delaunay2d_t *d)
{
  int i;
  int j;

  int error_count;

  double ba, bb, bc;

  error_count = 0;
  for (i = 0; i < d->nt; i ++) {

    for (j = 0; j < d->np; j ++) {

      if (j == d->t[i].v[0] ||
	  j == d->t[i].v[1] ||
	  j == d->t[i].v[2]) {
	continue;
      }

      if (barycentre_of_point(d,
			      i,
			      d->p[j].x,
			      d->p[j].y,
			      &ba,
			      &bb,
			      &bc) < 0) {
	fprintf(stderr, "delaunay2d_validate_nonintersecting: failed to compute barycentre\n");
	error_count ++;
	continue;
      }

      if (ba == 0.0 && (bb <= 1.0 && bc <= 1.0)) {
	fprintf(stderr, "delaunay2d_validate_nonintersecting: point %d in edge of triangle %d\n",
		j, i);
	error_count ++;
      } else if (bb == 0.0 && (ba <= 1.0 && bc <= 1.0)) {
	fprintf(stderr, "delaunay2d_validate_nonintersecting: point %d in edge of triangle %d\n",
		j, i);
	error_count ++;
      } else if (bc == 0.0 && (ba <= 1.0 && bb <= 1.0)) {
	fprintf(stderr, "delaunay2d_validate_nonintersecting: point %d in edge of triangle %d\n",
		j, i);
	error_count ++;
      } else if ((ba > 0.0 && bb > 0.0 && bc > 0.0) &&
		 (ba <= 1.0 && bb <= 1.0 && bc <= 1.0)) {
	fprintf(stderr, "delaunay2d_validate_nonintersecting: point %d inside triangle %d\n",
		j, i);
	error_count ++;
      }
    }
  }

  if (error_count == 0) {
    return 0;
  }

  return -1;
}

void
delaunay2d_print_points(const delaunay2d_t *d)
{
  int i;

  for (i = 0; i < d->np; i ++) {
    fprintf(stderr, "%4d: %10.6f %10.6f\n",
	    i, 
	    d->p[i].x,
	    d->p[i].y);
  }
}

void
delaunay2d_print_edges(const delaunay2d_t *d)
{
  int i;
  int j;

  for (i = 0; i < d->np; i ++) {
    fprintf(stderr, "%4d: %4d: ", i, d->e[i].n);

    for (j = 0; j < d->e[i].n; j ++) {
      fprintf(stderr, "%4d ", d->e[i].edge[j]);
    }
    fprintf(stderr, "\n");
  }
}


void
delaunay2d_print_triangles(const delaunay2d_t *d)
{
  int i;

  for (i = 0; i < d->nt; i ++) {
    fprintf(stderr, "%4d: %4d %4d %4d : %4d %4d %4d\n",
	    i,
	    d->t[i].v[0],
	    d->t[i].v[1],
	    d->t[i].v[2],
	    d->t[i].n[0],
	    d->t[i].n[1],
	    d->t[i].n[2]);
  }
}

static int *
resize_int_array(int *a, int size, int new_size)
{
  int *r;
  int i;

  r = malloc(sizeof(int) * new_size);
  RJMCMC_NULLCHECKPTR(r, "resize_int_array: failed to allocate memory\n");
  
  for (i = 0; i < size; i ++) {
    r[i] = a[i];
  }

  free(a);
  return r;
}

static int delete_triangle(delaunay2d_t *d,
			   int ti,
			   int *new_ti)
{
  int i;
  int nti;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "delete_triangle: null delaunay\n");
  RJMCMC_CONDITIONCHECKINT(ti < 0, "delete_triangle: invalid index\n");
  RJMCMC_CONDITIONCHECKINT(ti >= d->nt, "delete_triangle: invalid index (num triangles)\n");

  /*
   * Find the last triangle that isn't deleted.
   */
  for (i = d->nt - 1; i > ti; i --) {
    if (d->t[i].v[0] >= 0) {
      break;
    }
  }

  if (i > ti) {

    nti = i;
    *new_ti = nti;

    /*
     * Overwrite the deleted triangle 
     */
    d->t[ti].v[0] = d->t[nti].v[0];
    d->t[ti].v[1] = d->t[nti].v[1];
    d->t[ti].v[2] = d->t[nti].v[2];

    d->t[ti].n[0] = d->t[nti].n[0];
    d->t[ti].n[1] = d->t[nti].n[1];
    d->t[ti].n[2] = d->t[nti].n[2];

    d->t[ti].detT = d->t[nti].detT;
    d->t[ti].cx = d->t[nti].cx;
    d->t[ti].cy = d->t[nti].cy;
    d->t[ti].cr2 = d->t[nti].cr2;

    /*
     * Update references
     */
    if (d->t[ti].n[0] >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, d->t[ti].n[0], nti, ti) < 0,
		  "delete_triangle: failed to replace neighbour nab\n");
    }
    if (d->t[ti].n[1] >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, d->t[ti].n[1], nti, ti) < 0,
		  "delete_triangle: failed to replace neighbour nbc\n");
    }
    if (d->t[ti].n[2] >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_replace_neighbour(d, d->t[ti].n[2], nti, ti) < 0,
		  "delete_triangle: failed to replace neighbour nca\n");
    }
  } else {
    *new_ti = -1;
  }
  
  return 0;
}

static int delete_point(delaunay2d_t *d,
			int pi)
{
  int npi;
  int i;

  npi = d->np - 1;

  if (pi == npi) {
    /* Don't need to do anything */
  } else {
    /* Shift points down, updating triangles and edges */

    for (i = pi + 1; i < d->np; i ++) {
      d->p[i - 1].x = d->p[i].x;
      d->p[i - 1].y = d->p[i].y;
    }

    for (i = 0; i < d->nt; i ++) {

      if (d->t[i].v[0] >= pi) {
	d->t[i].v[0] --;
      }
      if (d->t[i].v[1] >= pi) {
	d->t[i].v[1] --;
      }
      if (d->t[i].v[2] >= pi) {
	d->t[i].v[2] --;
      }

    }

    if (shift_down_edge(d, pi) < 0) {
      return -1;
    }
  }

  return 0;
}

static int delete_edge(delaunay2d_t *d,
		       int p,
		       int np)
{
  int i;
  int ni;

  /*
   * Remove p -> np
   */
  ni = -1;

  for (i = 0; i < d->e[p].n; i ++) {
    if (d->e[p].edge[i] == np) {
      ni = i;
      break;
    }
  }

  if (ni >= 0) {
    for (i = ni + 1; i < d->e[p].n; i ++) {
      d->e[p].edge[i - 1] = d->e[p].edge[i];
    }
    d->e[p].n --;
  }
    
  /*
   * Remove np -> p
   */
  ni = -1;

  for (i = 0; i < d->e[np].n; i ++) {
    if (d->e[np].edge[i] == p) {
      ni = i;
      break;
    }
  }

  if (ni >= 0) {
    for (i = ni + 1; i < d->e[np].n; i ++) {
      d->e[np].edge[i - 1] = d->e[np].edge[i];
    }
    d->e[np].n --;
  }

  return 0;
}

static int flip_edge(delaunay2d_t *d,
		     int e1a, 
		     int e1b,
		     int e2a,
		     int e2b)
{
  if (delete_edge(d, e1a, e1b) < 0) {
    return -1;
  }

  if (add_edge(d, e2a, e2b) < 0) {
    return -1;
  }

  return 0;
}

static int shift_down_edge(delaunay2d_t *d,
			  int p)
{
  int i;
  int j;
  int np;

  RJMCMC_CONDITIONCHECKINT(d == NULL, "shiftdown_edge: null delaunay\n");
  RJMCMC_CONDITIONCHECKINT(p >= d->np, "shiftdown_edge: invalid index\n");

  for (i = p + 1; i < d->np; i ++) {

    d->e[i - 1].n = d->e[i].n;
    for (j = 0; j < d->e[i].n; j ++) {

      np = d->e[i].edge[j];
      if (np >= p) {
	np --;
      }
      d->e[i - 1].edge[j] = np;
    }
  }

  for (i = 0; i < p; i ++) {
    for (j = 0; j < d->e[i].n; j ++) {
      if (d->e[i].edge[j] >= p) {
	d->e[i].edge[j] --;
      }
    }
  }
    

  return 0;
}


static void rotate_int_array(int *a, int size, int steps)
{
  int i;
  int j;
  int first;

  for (j = 0; j < steps; j ++) {
    first = a[0];
    
    for (i = 1; i < size; i ++) {
      a[i - 1] = a[i];
    }
    
    a[size - 1] = first;
  }
}

#if 0
int
triangle_circumcircle(double x1, double y1,
		      double x2, double y2,
		      double x3, double y3,
		      double *cx, double *cy,
		      double *r2)
{
  double dx21;
  double dy21;
  double dx32;
  double dy32;
  double denom;
  double t;
  double dx;
  double dy;

  dx32 = x3 - x2;
  dx21 = x2 - x1;
  dy21 = y2 - y1;
  dy32 = y3 - y2;

  if (dx32 != 0.0) {
    denom = dy21*dx32 - dy32*dx21;
    if (denom == 0.0) {
      fprintf(stderr, "circumcircle_of_triangle: 0 denominator for dx32: (%f %f), (%f %f), (%f %f)\n",
	      x1, y1, 
	      x2, y2, 
	      x3, y3);
    }
    RJMCMC_CONDITIONCHECKINT(denom == 0.0, "circumcircle_of_triangle: 0 denominator for dx32 != 0\n");

    t = ((x3 - x1)*dx32 + (y3 - y1)*dy32)/(2.0 * denom);

    *cx = (x1 + x2)/2.0 + dy21*t;
    *cy = (y2 + y1)/2.0 - dx21*t;

    dx = (*cx - x1);
    dy = (*cy - y1);
    *r2 = dx*dx + dy*dy;

  } else if (dx21 != 0.0) {

    *cy = (y3 + y2)/2.0;
    t = (*cy + (y1 + y2)/2.0)/dx21;

    *cx = 0.5*(x1 + x2) + dy21*t;

    dx = (*cx - x1);
    dy = (*cy - y1);
    *r2 = dx*dx + dy*dy;

  } else {
    /* First 2 points of triangle coincident. */
    fprintf(stderr, "circumcircle_of_triangle: colinear points %f %f %f %f %f\n",
	    x1, x2, x3, dx21, dx32);
    return -1;
  }

  return 0;
}
#endif

int
triangle_circumcircle(double x1, double y1,
		      double x2, double y2,
		      double x3, double y3,
		      double *cx, double *cy,
		      double *r2)
{
  /*
   * From a discussion here: http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
   * Algorithm by Dave Watson
   */
  double D;
  double dx;
  double dy;

  D = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
  if (D == 0.0) {
    fprintf(stderr, "circumcircle_of_triangle: 0 determinant : (%f %f) (%f %f) (%f %f)\n",
	    x1, y1, 
	    x2, y2, 
	    x3, y3);
    /* fprintf(stderr, "  %f x %f - %f x %f\n", x1 - x3, y2 - y3, x2 - x3, y1 - y3); */
    return -1;
  }
  
  *cx = (((x1 - x3) * (x1 + x3) + (y1 - y3) * (y1 + y3)) / 2 * (y2 - y3) 
	 -  ((x2 - x3) * (x2 + x3) + (y2 - y3) * (y2 + y3)) / 2 * (y1 - y3)) 
    / D;

  *cy = (((x2 - x3) * (x2 + x3) + (y2 - y3) * (y2 + y3)) / 2 * (x1 - x3)
	 -  ((x1 - x3) * (x1 + x3) + (y1 - y3) * (y1 + y3)) / 2 * (x2 - x3))
    / D;

  dx = (*cx - x1);
  dy = (*cy - y1);
  *r2 = dx*dx + dy*dy;
  
  return 0;
}    


int
delaunay2d_all_visible(const delaunay2d_t *d,
		       int *vlist,
		       int nvertices)
{
  int i;

  double nx;
  double ny;

  double tx;
  double ty;
  
  for (i = 1; i < (nvertices - 1); i ++) {
    
    nx = d->p[vlist[i]].y - d->p[vlist[0]].y;
    ny = -(d->p[vlist[i]].x - d->p[vlist[0]].x);

    tx = d->p[vlist[i + 1]].x - d->p[vlist[i]].x;
    ty = d->p[vlist[i + 1]].y - d->p[vlist[i]].y;

    if ((nx*tx + ny*ty) < 0.0) {
      return 0;
    }
  }

  return -1;
}

int
delaunay2d_save_geo(const delaunay2d_t *d, const char *filename)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "delaunay2d_save_geo: failed to open file\n");
    return -1;
  }

  fprintf(fp, "PGEOMETRY V5\n");
  fprintf(fp, "NPoints %d NPrims %d\n", d->np, d->nt);
  fprintf(fp, "NPointGroups 0 NPrimGroups 0\n");
  fprintf(fp, "NPointAttrib 0 NVertexAttrib 0 NPrimAttrib 0 NAttrib 0\n");
  
  for (i = 0; i < d->np; i ++) {
    fprintf(fp, "%f %f 0.0 1.0\n", d->p[i].x, d->p[i].y);
  }

  fprintf(fp, "Run %d Poly\n", d->nt);
  for (i = 0; i < d->nt; i ++) {
    fprintf(fp, " 3 < %d %d %d\n", d->t[i].v[0], d->t[i].v[1], d->t[i].v[2]);
  }

  fprintf(fp, "beginExtra\n");
  fprintf(fp, "endExtra\n");

  fclose(fp);

  return 0;
}

int
delaunay2d_save_cc_geo(const delaunay2d_t *d, const char *filename)
{
  FILE *fp;
  int i;
  double r;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "delaunay2d_save_cc_geo: failed to open file\n");
    return -1;
  }

  fprintf(fp, "PGEOMETRY V5\n");
  fprintf(fp, "NPoints %d NPrims %d\n", d->nt, d->nt);
  fprintf(fp, "NPointGroups 0 NPrimGroups 0\n");
  fprintf(fp, "NPointAttrib 0 NVertexAttrib 0 NPrimAttrib 0 NAttrib 0\n");
  
  for (i = 0; i < d->nt; i ++) {
    fprintf(fp, "%f %f 0.0 1.0\n", d->t[i].cx, d->t[i].cy);
  }

  for (i = 0; i < d->nt; i ++) {
    r = sqrt(d->t[i].cr2);
    fprintf(fp, "Circle %d %f 0 0 0 %f 0 0 0 1\n",
	    i,
	    r,
	    r);
  }

  fprintf(fp, "beginExtra\n");
  fprintf(fp, "endExtra\n");

  fclose(fp);

  return 0;
}


static int delaunay2d_can_close_first_triangle(const delaunay2d_t *d, 
					       int nd)
{
  int i;
  double nx;
  double ny;

  double bcx;
  double bcy;

  /*
   * First test if angle abc is oriented correctly
   */
  nx = d->p[d->vlist[1]].y - d->p[d->vlist[0]].y;
  ny = -(d->p[d->vlist[1]].x - d->p[d->vlist[0]].x);

  bcx = d->p[d->vlist[2]].x - d->p[d->vlist[1]].x;
  bcy = d->p[d->vlist[2]].y - d->p[d->vlist[1]].y;

  if ((nx*bcx + ny*bcy) <= 0.0) {
    /* Angle abc is reflex/straight so can't form a valid triangle */
    return 0;
  }

  /*
   * Now check if other points are within the proposed triangle
   */
  for (i = 3; i < nd; i ++) {
    if (point_in_triangle(d->p[d->vlist[i]].x, d->p[d->vlist[i]].y,
			  d->p[d->vlist[0]].x, d->p[d->vlist[0]].y,
			  d->p[d->vlist[1]].x, d->p[d->vlist[1]].y,
			  d->p[d->vlist[2]].x, d->p[d->vlist[2]].y)) {
      return 0;
    }
  }

  return -1;
}

static int delaunay2d_fill_hole(delaunay2d_t *d, 
				int nd)
{
  int i;
  int ti;
  int ni;

  if (nd == 3) {
    /* Replace triangle 0 with triangle made from vertices 0, 1, 2 */

    ti = d->tlist[0];
    d->t[ti].v[0] = d->vlist[0];
    d->t[ti].v[1] = d->vlist[1];
    d->t[ti].v[2] = d->vlist[2];


    d->t[ti].n[0] = d->nlist[0];
    if (d->t[ti].n[0] >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_verify_neighbour(d, ti, 0) < 0, "delaunay2d_fill_hole: failed to verify neighbour ab\n");
    }
    d->t[ti].n[1] = d->nlist[1];
    if (d->t[ti].n[1] >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_verify_neighbour(d, ti, 1) < 0, "delaunay2d_fill_hole: failed to verify neighbour bc\n");
    }
    d->t[ti].n[2] = d->nlist[2];
    if (d->t[ti].n[2] >= 0) {
      RJMCMC_CONDITIONCHECKINT(triangle_verify_neighbour(d, ti, 2) < 0, "delaunay2d_fill_hole: failed to verify neighbour ca\n");
    }


    /* printf("fill: filling %d with %d %d %d (%d %d %d)\n", ti,  */
    /* 	   d->t[ti].v[0], d->t[ti].v[1], d->t[ti].v[2], */
    /* 	   d->t[ti].n[0], d->t[ti].n[1], d->t[ti].n[2]); */

    for (i = 1; i < nd; i ++) {
      d->tlist[i - 1] = d->tlist[i];
    }

    RJMCMC_CONDITIONCHECKINT(triangle_update(d, ti) < 0, "delaunay2d_fill_hole: failed to update filled triangle\n");

    RJMCMC_CONDITIONCHECKINT(queue_flip(d, ti, 0) < 0, "delaunay2d_fill_hole: failed to queue flip fill 0\n");
    RJMCMC_CONDITIONCHECKINT(queue_flip(d, ti, 1) < 0, "delaunay2d_fill_hole: failed to queue flip fill 1\n");
    RJMCMC_CONDITIONCHECKINT(queue_flip(d, ti, 2) < 0, "delaunay2d_fill_hole: failed to queue flip fill 2\n");

    ni = neighbour_index(d, d->t[ti].n[0], ti);
    if (ni >= 0) {
      RJMCMC_CONDITIONCHECKINT(queue_flip(d, d->t[ti].n[0], ni) < 0,
		  "delaunay2d_fill_hole: failed to queue flip fill inverse nab\n");
    }
    
    ni = neighbour_index(d, d->t[ti].n[1], ti);
    if (ni >= 0) {
      RJMCMC_CONDITIONCHECKINT(queue_flip(d, d->t[ti].n[1], ni) < 0,
		  "delaunay2d_fill_hole: failed to queue flip fill inverse nbc\n");
    }

    ni = neighbour_index(d, d->t[ti].n[2], ti);
    if (ni >= 0) {
      RJMCMC_CONDITIONCHECKINT(queue_flip(d, d->t[ti].n[2], ni) < 0,
		  "delaunay2d_fill_hole: failed to queue flip fill inverse nca\n");
    }

    return 0;

  } else {
    if (delaunay2d_can_close_first_triangle(d, nd)) {

      ti = d->tlist[0];
      d->t[ti].v[0] = d->vlist[0];
      d->t[ti].v[1] = d->vlist[1];
      d->t[ti].v[2] = d->vlist[2];

      /*
       * Edge 0->1 and 1->2 already exist, but need to add edge 2->0
       */
      if (add_edge(d, d->t[ti].v[2], d->t[ti].v[0]) < 0) {
	fprintf(stderr, "delaunay2d_fill_hole: failed to add edge 2 0\n");
	return -1;
      }
      
      d->t[ti].n[0] = d->nlist[0];
      if (d->t[ti].n[0] >= 0) {
	RJMCMC_CONDITIONCHECKINT(triangle_verify_neighbour(d, ti, 0) < 0, 
				 "delaunay2d_fill_hole: failed to verify neighbour ab\n");
      }
      d->t[ti].n[1] = d->nlist[1];
      if (d->t[ti].n[1] >= 0) {
	RJMCMC_CONDITIONCHECKINT(triangle_verify_neighbour(d, ti, 1) < 0, 
				 "delaunay2d_fill_hole: failed to verify neighbour bc\n");
      }
      d->t[ti].n[2] = -1;

      /* printf("fill: replacing %d with %d %d %d (%d %d %d)\n", ti,  */
      /* 	     d->t[ti].v[0], d->t[ti].v[1], d->t[ti].v[2], */
      /* 	     d->t[ti].n[0], d->t[ti].n[1], d->t[ti].n[2]); */

      /*
       * Remove vertices and neighbours that have been used up 
       */
      d->nlist[0] = ti;

      for (i = 1; i < nd; i ++) {
	d->tlist[i - 1] = d->tlist[i];
      }

      for (i = 2; i < nd; i ++) {
	d->vlist[i - 1] = d->vlist[i];
	d->nlist[i - 1] = d->nlist[i];
      }

      RJMCMC_CONDITIONCHECKINT(triangle_update(d, ti) < 0, "delaunay2d_fill_hole: failed to update replaced triangle\n");

      RJMCMC_CONDITIONCHECKINT(queue_flip(d, ti, 0) < 0, "delaunay2d_fill_hole: failed to queue flip 0\n");
      RJMCMC_CONDITIONCHECKINT(queue_flip(d, ti, 1) < 0, "delaunay2d_fill_hole: failed to queue flip 1\n");
      RJMCMC_CONDITIONCHECKINT(queue_flip(d, ti, 2) < 0, "delaunay2d_fill_hole: failed to queue flip 2\n");

      return delaunay2d_fill_hole(d, nd - 1);

    } else {
      
      rotate_int_array(d->vlist, nd, 1);
      rotate_int_array(d->nlist, nd, 1);
      rotate_int_array(d->tlist, nd, 1);

      return delaunay2d_fill_hole(d, nd);
    }

  }
  return -1;
}

static int clear_flip(delaunay2d_t *d)
{
  d->nf = 0;
  return 0;
}

static int queue_flip(delaunay2d_t *d, 
		      int ti,
		      int ei)
{
  if (d->nf == d->sf) {
    int newsize = d->sf + DLIST_INCR;
    int *newarrayt = malloc(sizeof(int) * newsize);
    int *newarraye = malloc(sizeof(int) * newsize);
    int i;

    if (newarrayt == NULL ||
	newarraye == NULL) {
      return -1;
    }

    for (i = 0; i < d->nf; i ++) {
      newarrayt[i] = d->ftlist[i];
      newarraye[i] = d->felist[i];
    }
    
    free(d->ftlist);
    free(d->felist);

    d->ftlist = newarrayt;
    d->felist = newarraye;

    d->sf = newsize;
  }

  d->ftlist[d->nf] = ti;
  d->felist[d->nf] = ei;
  d->nf ++;

  return 0;
}

static int update_flip(delaunay2d_t *d,
		       int old_ti,
		       int new_ti)
{
  int i;

  for (i = 0; i < d->nf; i ++) {
    if (d->ftlist[i] == old_ti) {
      d->ftlist[i] = new_ti;
    }
  }

  return 0;
}

static int fill_fan_about_point(delaunay2d_t *d,
				int pi,
				int ti,
				bbox2d_t *bound)
{
  int i;
  int ni;
  
  bbox2d_initialize(bound, d->p[pi].x, d->p[pi].y);
  
  /*
   * Determine the other triangles in the polygon that surrounds the point
   */
  ni = ti;
  i = 0;
  do {

    /* Add ni to the list of triangles to remove */
    d->tlist[i] = ni;

    /* Determine next triangle in clockwise direction */
    if (d->t[ni].v[0] == pi) {

      /* add point to ring d->t[ni].v[1]*/
      d->vlist[i] = d->t[ni].v[1];

      /* add neighbour to ring d->t[ni].n[1] */
      d->nlist[i] = d->t[ni].n[1];

      /* Next triangle is the ca neighbor */
      ni = d->t[ni].n[2];

    } else if (d->t[ni].v[1] == pi) {
      /* add point to ring d->t[ni].v[2]*/
      d->vlist[i] = d->t[ni].v[2];
      
      /* add neighbour to ring d->t[ni].n[2] */
      d->nlist[i] = d->t[ni].n[2];

      /* Next triangle is the ab neighbor */
      ni = d->t[ni].n[0];

    } else if (d->t[ni].v[2] == pi) {
      /* add point to ring d->t[ni].v[0]*/
      d->vlist[i] = d->t[ni].v[0];

      /* add neighbour to ring d->t[ni].n[0] */
      d->nlist[i] = d->t[ni].n[0];

      /* Next triangle is the bc neighbor */
      ni = d->t[ni].n[1];
    } else {
      fprintf(stderr, "delaunay2d_delete: neighbour doesn't contain point %d (%d %d %d) %d\n",
	      ni, 
	      d->t[ni].v[0], d->t[ni].v[1], d->t[ni].v[2], 
	      pi);
      return -1;
    }

    bbox2d_expand(bound, d->p[d->vlist[i]].x, d->p[d->vlist[i]].y);

    i ++;
    if (i >= d->listsize) {
      /* Resize the lists */
      d->vlist = resize_int_array(d->vlist, d->listsize, d->listsize + DLIST_INCR);
      RJMCMC_NULLCHECKINT(d->vlist, "delaunay2d_delete: null vlist");
      d->nlist = resize_int_array(d->nlist, d->listsize, d->listsize + DLIST_INCR);
      RJMCMC_NULLCHECKINT(d->vlist, "delaunay2d_delete: null nlist");
      d->tlist = resize_int_array(d->tlist, d->listsize, d->listsize + DLIST_INCR);
      RJMCMC_NULLCHECKINT(d->vlist, "delaunay2d_delete: null tlist");

      d->listsize += DLIST_INCR;
    }
  } while (ni != ti);

  return i;
}

static double dist2(delaunay2d_t *d,
		    int pi,
		    double px,
		    double py)
{
  double dx;
  double dy;

  if (d->ci != d->p[pi].ci) {
    dx = px - d->p[pi].x;
    dy = py - d->p[pi].y;
    
    d->p[pi].pd2 = dx*dx + dy*dy;
    d->p[pi].ci = d->ci;
  }

  return d->p[pi].pd2;
}
