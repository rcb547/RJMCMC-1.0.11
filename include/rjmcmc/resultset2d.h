#ifndef resultset2d_h
#define resultset2d_h

typedef enum {
  RESULTSET2D_MEAN       = 0x01,
  RESULTSET2D_MEDIAN     = 0x02,
  RESULTSET2D_MODE       = 0x04,
  RESULTSET2D_CREDIBLE   = 0x08,

  RESULTSET2D_GRADIENT   = 0x10
} resultset2d_result_t;

typedef struct _resultset2d resultset2d_t;

resultset2d_t *
resultset2d_create(int burnin,
		   int total,
		   int xsamples,
		   int ysamples,
		   int zsamples,
		   int nprocesses,
		   int maxpartitions,
		   double xmin,
		   double xmax,
		   double ymin,
		   double ymax,
		   double zmin,
		   double zmax,
		   double credible_interval,
		   int results);

void 
resultset2d_destroy(resultset2d_t *r);


/*
 * Sampling methods
 */

void
resultset2d_propose(resultset2d_t *r,
		    int p);

void
resultset2d_accept(resultset2d_t *r,
		   int p);

void
resultset2d_sample(resultset2d_t *r,
		   int i,
		   const double **v);

void 
resultset2d_sample_misfit(resultset2d_t *r,
			  int i,
			  double misfit);

void
resultset2d_sample_centre(resultset2d_t *r,
			  double x,
			  double y);

void 
resultset2d_sample_sigma(resultset2d_t *r,
			 int i,
			 double sigma);

void
resultset2d_sample_npartitions(resultset2d_t *r,
			       int i,
			       int npartitions);

/*
 * Assemble after simulation (needed for histogram derived values, eg median)
 */

void
resultset2d_assemble_results(resultset2d_t *r);

/*
 * Getting results
 */

const int *
resultset2d_get_propose(resultset2d_t *r,
			int *nprocesses);

const int *
resultset2d_get_accept(resultset2d_t *r,
		       int *nprocesses);

const double *
resultset2d_get_misfit(resultset2d_t *r);

const double *
resultset2d_get_lambda(resultset2d_t *r);

const int *
resultset2d_get_partitions(resultset2d_t *r);

const double **
resultset2d_get_mean(resultset2d_t *r);

const int **
resultset2d_get_centres(resultset2d_t *r);

const double **
resultset2d_get_median(resultset2d_t *r);

const double **
resultset2d_get_mode(resultset2d_t *r);

const double **
resultset2d_get_credible_min(resultset2d_t *r);
		     
const double **
resultset2d_get_credible_max(resultset2d_t *r);

const double *
resultset2d_get_user_ensemble(resultset2d_t *r,
			      int ei);

void
resultset2d_fill_xcoord_vector(resultset2d_t *r,
			       double *x,
			       int *l);

void
resultset2d_fill_ycoord_vector(resultset2d_t *r,
			       double *y,
			       int *l);

void
resultset2d_fill_zcoord_vector(resultset2d_t *r,
			       double *z,
			       int *l);

#endif /* resultset_h */
