
#include <stdio.h>
#include <stdlib.h>

#include "rjmcmc/resultset2dfm.h"

#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_defines.h"

struct _resultset2dfm {

  int results;

  int burnin;
  int total;
  int thin;

  int xsamples;
  int ysamples;
  int zsamples;

  int nprocesses;

  int max_partitions;

  double xmin;
  double xmax;
  double ymin;
  double ymax;

  int nglobalparameters;
  const forwardmodelparameter_t *global_parameters;
  
  int nlocalparameters;
  const forwardmodelparameter_t *local_parameters;

  int *propose;
  int *accept;

  /*
   * Misfit [total]
   */
  double *misfit;

  /*
   * N Partitions [total]
   */
  int *partitions;

  /*
   * Global Parameter History [nglobalparameters][total];
   */
  double **global;

  /*
   * Hierarchical Parameter History [nhierarchicalparameters][total];
   */
  int nhierarchical;
  double **hierarchical;

  /*
   * Centres [xsamples][ysamples]
   */
  int **centres;

  /*
   * Histogram [nlocalparameters][xsamples][ysamples][zsamples]
   */
  int ****hist;

  /*
   * Mean [nlocalparameters][xsamples][ysamples]
   */
  double ***local_mean;
  double ***local_variance;
  int local_samples;

  /*
   * Mode [nlocalparameters][xsamples][ysamples]
   */
  double ***local_mode;

  /*
   * Median [nlocalparameters][xsamples][ysamples]
   */
  double ***local_median;

  /*
   * Credible Intervals [nlocalparameters][xsamples][ysamples]
   */
  double credible_interval;
  double ***local_cred_min;
  double ***local_cred_max;


};

typedef struct {
  int (*write_header)(FILE *fp);
  int (*write_variable_header)(FILE *fp, const char *name);
  int (*write_variable_footer)(FILE *fp);
  int (*write_int)(FILE *fp, int i);
  int (*write_int_array)(FILE *fp, int *a, int size);
  int (*write_int_array_2d)(FILE *fp, int **a, int width, int height);
  int (*write_int_array_3d)(FILE *fp, int ***a, int width, int height, int depth);
  int (*write_int_array_4d)(FILE *fp, int ****a, int width, int height, int depth, int length);
  int (*write_double)(FILE *fp, double d);
  int (*write_double_array)(FILE *fp, double *a, int size);
  int (*write_double_array_2d)(FILE *fp, double **a, int width, int height);
  int (*write_double_array_3d)(FILE *fp, double ***a, int width, int height, int depth);
  int (*write_footer)(FILE *fp);
} writer_t;

typedef struct {
  int (*read_header)(FILE *fp);
  int (*read_variable_header)(FILE *fp, const char *name);
  int (*read_variable_footer)(FILE *fp);
  int (*read_int)(FILE *fp, int *i);
  int (*read_int_array)(FILE *fp, int *a, int size);
  int (*read_int_array_2d)(FILE *fp, int **a, int width, int height);
  int (*read_int_array_3d)(FILE *fp, int ***a, int width, int height, int depth);
  int (*read_int_array_4d)(FILE *fp, int ****a, int width, int height, int depth, int length);
  int (*read_double)(FILE *fp, double *d);
  int (*read_double_array)(FILE *fp, double *a, int size);
  int (*read_double_array_2d)(FILE *fp, double **a, int width, int height);
  int (*read_double_array_3d)(FILE *fp, double ***a, int width, int height, int depth);
  int (*read_footer)(FILE *fp);
} reader_t;

static int binary_write_header(FILE *fp);
static int binary_write_variable_header(FILE *fp, const char *name);
static int binary_write_variable_footer(FILE *fp);
static int binary_write_int(FILE *fp, int i);
static int binary_write_int_array(FILE *fp, int *a, int size);
static int binary_write_int_array_2d(FILE *fp, int **a, int width, int height);
static int binary_write_int_array_3d(FILE *fp, int ***a, int width, int height, int depth);
static int binary_write_int_array_4d(FILE *fp, int ****a, int width, int height, int depth, int length);
static int binary_write_double(FILE *fp, double d);
static int binary_write_double_array(FILE *fp, double *a, int size);
static int binary_write_double_array_2d(FILE *fp, double **a, int width, int height);
static int binary_write_double_array_3d(FILE *fp, double ***a, int width, int height, int depth);
static int binary_write_footer(FILE *fp);

static int binary_read_header(FILE *fp);
static int binary_read_variable_header(FILE *fp, const char *name);
static int binary_read_variable_footer(FILE *fp);
static int binary_read_int(FILE *fp, int *i);
static int binary_read_int_array(FILE *fp, int *a, int size);
static int binary_read_int_array_2d(FILE *fp, int **a, int width, int height);
static int binary_read_int_array_3d(FILE *fp, int ***a, int width, int height, int depth);
static int binary_read_int_array_4d(FILE *fp, int ****a, int width, int height, int depth, int length);
static int binary_read_double(FILE *fp, double *d);
static int binary_read_double_array(FILE *fp, double *a, int size);
static int binary_read_double_array_2d(FILE *fp, double **a, int width, int height);
static int binary_read_double_array_3d(FILE *fp, double ***a, int width, int height, int depth);
static int binary_read_footer(FILE *fp);

static writer_t binary_writer = {
  binary_write_header,
  binary_write_variable_header,
  binary_write_variable_footer,
  binary_write_int,
  binary_write_int_array,
  binary_write_int_array_2d,
  binary_write_int_array_3d,
  binary_write_int_array_4d,
  binary_write_double,
  binary_write_double_array,
  binary_write_double_array_2d,
  binary_write_double_array_3d,
  binary_write_footer
};

static reader_t binary_reader = {
  binary_read_header,
  binary_read_variable_header,
  binary_read_variable_footer,
  binary_read_int,
  binary_read_int_array,
  binary_read_int_array_2d,
  binary_read_int_array_3d,
  binary_read_int_array_4d,
  binary_read_double,
  binary_read_double_array,
  binary_read_double_array_2d,
  binary_read_double_array_3d,
  binary_read_footer
};

resultset2dfm_t *
resultset2dfm_create(int burnin,
		     int total,
		     int thin,
		     int nglobalparameters,
		     const forwardmodelparameter_t *global_parameters,
		     int nlocalparameters,
		     const forwardmodelparameter_t *local_parameters,
		     int nhierarchicalparameters,
		     int xsamples,
		     int ysamples,
		     int zsamples,
		     int maxpartitions,
		     double xmin,
		     double xmax,
		     double ymin,
		     double ymax,
		     int nprocesses,
		     double credible_interval,
		     int results)
{
  resultset2dfm_t *r;
  int i;

  r = (resultset2dfm_t *)malloc(sizeof(resultset2dfm_t));
  if (r == NULL) {
    rjmcmc_fatal("resultset2dfm_create: failed to allocate memory\n");
    return NULL;
  }

  r->results = results;
  r->burnin = burnin;
  r->total = total;
  r->thin = thin;

  r->xsamples = xsamples;
  r->ysamples = ysamples;
  r->zsamples = zsamples;

  r->nprocesses = nprocesses;

  r->max_partitions = maxpartitions;

  r->xmin = xmin;
  r->xmax = xmax;
  r->ymin = ymin;
  r->ymax = ymax;

  r->nglobalparameters = nglobalparameters;
  r->global_parameters = global_parameters;

  r->nlocalparameters = nlocalparameters;
  r->local_parameters = local_parameters;

  r->credible_interval = credible_interval;
  
  r->propose = rjmcmc_create_int_array_1d(nprocesses);
  if (r->propose == NULL) {
    return NULL;
  }
  
  r->accept = rjmcmc_create_int_array_1d(nprocesses);
  if (r->accept == NULL) {
    return NULL;
  }
  
  r->misfit = rjmcmc_create_array_1d(total);
  if (r->misfit == NULL) {
    return NULL;
  }

  r->partitions = rjmcmc_create_int_array_1d(total);
  if (r->partitions == NULL) {
    return NULL;
  }
  
  r->global = NULL;
  if (nglobalparameters > 0) {
    r->global = rjmcmc_create_array_2d(nglobalparameters, total);
    if (r->global == NULL) {
      return NULL;
    }
  }

  r->hierarchical = NULL;
  r->nhierarchical = nhierarchicalparameters;
  if (nhierarchicalparameters > 0) {
    r->hierarchical = rjmcmc_create_array_2d(nhierarchicalparameters, total);
    if (r->hierarchical == NULL) {
      return NULL;
    }
  }

  r->local_mean = rjmcmc_create_array_3d(nlocalparameters, xsamples, ysamples);
  if (r->local_mean == NULL) {
    return NULL;
  }
  r->local_variance = rjmcmc_create_array_3d(nlocalparameters, xsamples, ysamples);
  if (r->local_variance == NULL) {
    return NULL;
  }
  r->local_samples = 0;
   
  r->centres = rjmcmc_create_int_array_2d(xsamples, ysamples);
  if (r->centres == NULL) {
    return NULL;
  }

  /*
   * Histogram [nlocalparameters][xsamples][ysamples][zsamples]
   */
  r->hist = NULL;
  if ((results & RESULTSET2DFM_MODE) ||
      (results & RESULTSET2DFM_MEDIAN) ||
      (results & RESULTSET2DFM_CREDIBLE)) {
    r->hist = rjmcmc_create_int_array_4d(nlocalparameters,
					 xsamples,
					 ysamples,
					 zsamples);
    if (r->hist == NULL) {
      return NULL;
    }
  }

  r->local_mode = NULL;
  if ((results & RESULTSET2DFM_MODE)) {
    r->local_mode = rjmcmc_create_array_3d(nlocalparameters,
					   xsamples,
					   ysamples);
    if (r->local_mode == NULL) {
      return NULL;
    }
  }

  r->local_median = NULL;
  if ((results & RESULTSET2DFM_MEDIAN)) {
    r->local_median = rjmcmc_create_array_3d(nlocalparameters,
					     xsamples,
					     ysamples);
    if (r->local_median == NULL) {
      return NULL;
    }
  }

  r->local_cred_min = NULL;
  r->local_cred_max = NULL;
  if ((results & RESULTSET2DFM_CREDIBLE)) {
    r->local_cred_min = rjmcmc_create_array_3d(nlocalparameters,
					     xsamples,
					     ysamples);
    if (r->local_cred_min == NULL) {
      return NULL;
    }
    r->local_cred_max = rjmcmc_create_array_3d(nlocalparameters,
					       xsamples,
					       ysamples);
    if (r->local_cred_max == NULL) {
      return NULL;
    }
  }

  return r;
}

#define WRITE_INT(name, var) \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_int(fp, var), \
		     "resultset2dfm_save: failed to write int\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");

#define WRITE_INT_ARRAY(name, var, size) \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_int_array(fp, var, size),			\
		     "resultset2dfm_save: failed to write int array\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");

#define WRITE_INT_ARRAY_2D(name, var, width, height)		 \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_int_array_2d(fp, var, width, height),		\
		     "resultset2dfm_save: failed to write int array\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");

#define WRITE_INT_ARRAY_3D(name, var, width, height, depth)	 \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_int_array_3d(fp, var, width, height, depth),	\
		     "resultset2dfm_save: failed to write int array\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");

#define WRITE_INT_ARRAY_4D(name, var, width, height, depth, length)	\
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_int_array_4d(fp, var, width, height, depth, length), \
		     "resultset2dfm_save: failed to write int array\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");

#define WRITE_DOUBLE(name, var) \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_double(fp, var), \
		     "resultset2dfm_save: failed to write int\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");
  
#define WRITE_DOUBLE_ARRAY(name, var, size) \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_double_array(fp, var, size),			\
		     "resultset2dfm_save: failed to write double array\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");

#define WRITE_DOUBLE_ARRAY_2D(name, var, width, height)		 \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_double_array_2d(fp, var, width, height),		\
		     "resultset2dfm_save: failed to write double array\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");

#define WRITE_DOUBLE_ARRAY_3D(name, var, width, height, depth) \
  RJMCMC_INTCHECKINT(w->write_variable_header(fp, name), \
		     "resultset2dfm_save: failed to write variable header\n"); \
  RJMCMC_INTCHECKINT(w->write_double_array_3d(fp, var, width, height, depth),	\
		     "resultset2dfm_save: failed to write double array\n"); \
  RJMCMC_INTCHECKINT(w->write_variable_footer(fp), \
		     "resultset2dfm_save: failed to write variable footer\n");


int resultset2dfm_save(const resultset2dfm_t *r,
		       const char *filename,
		       int fmt)
{
  writer_t *w;
  FILE *fp;
  int i;

  switch(fmt) {
  case RESULTSET2DFM_BINARY:
    w = &binary_writer;
    break;

  default:
    rjmcmc_error("resultset2dfm_save: unimplemented format.\n");
    return -1;
  }

  fp = fopen(filename, "w");
  RJMCMC_NULLCHECKINT(fp, "resultset2dfm_save: failed to create file\n");

  RJMCMC_INTCHECKINT(w->write_header(fp), 
		     "resultset2dfm_save: failed to write header\n");

  WRITE_INT("results", r->results);
  WRITE_INT("burnin", r->burnin);
  WRITE_INT("total", r->total);
  WRITE_INT("thin", r->thin);
  WRITE_INT("xsamples", r->xsamples);
  WRITE_INT("ysamples", r->ysamples);
  WRITE_INT("zsamples", r->zsamples);
  WRITE_INT("nprocesses", r->nprocesses);
  WRITE_INT("maxpartitions", r->max_partitions);
  WRITE_DOUBLE("xmin", r->xmin);
  WRITE_DOUBLE("xmax", r->xmax);
  WRITE_DOUBLE("ymin", r->ymin);
  WRITE_DOUBLE("ymax", r->ymax);
  WRITE_DOUBLE("credible_interval", r->credible_interval);
  WRITE_INT("nhierarchical", r->nhierarchical);

  WRITE_INT("nglobalparameters", r->nglobalparameters);
  for (i = 0; i < r->nglobalparameters; i ++) {
    WRITE_DOUBLE("min", r->global_parameters[i].fmin);
    WRITE_DOUBLE("max", r->global_parameters[i].fmax);
    WRITE_DOUBLE("std", r->global_parameters[i].fstd_value);
    WRITE_DOUBLE("std_bd", r->global_parameters[i].fstd_bd);
  }

  WRITE_INT("nlocalparameters", r->nlocalparameters);
  for (i = 0; i < r->nlocalparameters; i ++) {
    WRITE_DOUBLE("min", r->local_parameters[i].fmin);
    WRITE_DOUBLE("max", r->local_parameters[i].fmax);
    WRITE_DOUBLE("std", r->local_parameters[i].fstd_value);
    WRITE_DOUBLE("std_bd", r->local_parameters[i].fstd_bd);
  }

  WRITE_INT_ARRAY("propose", r->propose, r->nprocesses);
  WRITE_INT_ARRAY("accept", r->accept, r->nprocesses);

  WRITE_DOUBLE_ARRAY("misfit", r->misfit, r->total);
    
  WRITE_INT_ARRAY("partitions", r->partitions, r->total);

  WRITE_DOUBLE_ARRAY_2D("global", r->global, r->nglobalparameters, r->total);

  WRITE_DOUBLE_ARRAY_2D("hierarchical", r->hierarchical, r->nhierarchical, r->total);

  WRITE_INT_ARRAY_2D("centres", r->centres, r->xsamples, r->ysamples);

  if ((r->results & (RESULTSET2DFM_MODE | RESULTSET2DFM_CREDIBLE | RESULTSET2DFM_MEDIAN)) > 0) {
    WRITE_INT_ARRAY_4D("hist", r->hist, r->nlocalparameters, r->xsamples, r->ysamples, r->zsamples);
  }

  if ((r->results & RESULTSET2DFM_MEAN) > 0) {
    WRITE_DOUBLE_ARRAY_3D("local_mean", r->local_mean, r->nlocalparameters, r->xsamples, r->ysamples);
    WRITE_DOUBLE_ARRAY_3D("local_delta", r->local_variance, r->nlocalparameters, r->xsamples, r->ysamples);
  }
  WRITE_INT("local_samples", r->local_samples);

  if ((r->results & RESULTSET2DFM_MODE) > 0) {
    WRITE_DOUBLE_ARRAY_3D("local_mode", r->local_mode, r->nlocalparameters, r->xsamples, r->ysamples);
  }

  if ((r->results & RESULTSET2DFM_MEDIAN) > 0) {
    WRITE_DOUBLE_ARRAY_3D("local_median", r->local_median, r->nlocalparameters, r->xsamples, r->ysamples);
  }
  
  if ((r->results & RESULTSET2DFM_CREDIBLE) > 0) {
    WRITE_DOUBLE_ARRAY_3D("local_cred_min", r->local_cred_min, r->nlocalparameters, r->xsamples, r->ysamples);
    WRITE_DOUBLE_ARRAY_3D("local_cred_max", r->local_cred_max, r->nlocalparameters, r->xsamples, r->ysamples);
  }

  RJMCMC_INTCHECKINT(w->write_footer(fp), 
		     "resultset2dfm_save: failed to write footer\n");
  fclose(fp);
  return 0;
}

#define READ_INT(name, var) \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_int(fp, var), \
		     "resultset2dfm_save: failed to read int\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");

#define READ_INT_ARRAY(name, var, size) \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_int_array(fp, var, size),			\
		     "resultset2dfm_save: failed to read int array\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");

#define READ_INT_ARRAY_2D(name, var, width, height)		 \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_int_array_2d(fp, var, width, height),		\
		     "resultset2dfm_save: failed to read int array\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");

#define READ_INT_ARRAY_3D(name, var, width, height, depth)	 \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_int_array_3d(fp, var, width, height, depth),	\
		     "resultset2dfm_save: failed to read int array\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");

#define READ_INT_ARRAY_4D(name, var, width, height, depth, length)	\
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_int_array_4d(fp, var, width, height, depth, length), \
		     "resultset2dfm_save: failed to read int array\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");

#define READ_DOUBLE(name, var) \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_double(fp, var), \
		     "resultset2dfm_save: failed to read int\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");
  
#define READ_DOUBLE_ARRAY(name, var, size) \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_double_array(fp, var, size),			\
		     "resultset2dfm_save: failed to read double array\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");

#define READ_DOUBLE_ARRAY_2D(name, var, width, height)		 \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_double_array_2d(fp, var, width, height),		\
		     "resultset2dfm_save: failed to read double array\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");

#define READ_DOUBLE_ARRAY_3D(name, var, width, height, depth) \
  RJMCMC_INTCHECKPTR(r->read_variable_header(fp, name), \
		     "resultset2dfm_save: failed to read variable header\n"); \
  RJMCMC_INTCHECKPTR(r->read_double_array_3d(fp, var, width, height, depth),	\
		     "resultset2dfm_save: failed to read double array\n"); \
  RJMCMC_INTCHECKPTR(r->read_variable_footer(fp), \
		     "resultset2dfm_save: failed to read variable footer\n");
 

resultset2dfm_t *resultset2dfm_load(const char *filename, int fmt,
				    int *results,
				    int *burnin,
				    int *total,
				    int *thin,
				    int *xsamples,
				    int *ysamples,
				    int *zsamples,
				    int *nprocesses,
				    int *maxpartitions,
				    double *xmin,
				    double *xmax,
				    double *ymin,
				    double *ymax,
				    double *credibleinterval,
				    int *nhierarchical,
				    int *nglobal,
				    int *nlocal,
				    forwardmodelparameter_t **global_parameters,
				    forwardmodelparameter_t **local_parameters)
{
  FILE *fp;
  reader_t *r;

  /*
   * Boot strap variables.
   */
  forwardmodelparameter_t *global;
  forwardmodelparameter_t *local;

  resultset2dfm_t *rs;
  int i;
  
  
  fp = fopen(filename, "r");
  RJMCMC_NULLCHECKPTR(fp, "resultset2dfm_load: failed to open file\n");
  
  switch(fmt) {
  case RESULTSET2DFM_BINARY:
    r = &binary_reader;
    break;

  default:
    rjmcmc_error("resultset2dfm_load: unimplemented format.\n");
    return NULL;
  }

  RJMCMC_INTCHECKPTR(r->read_header(fp),
		     "resultset2dfm_load: failed to read header\n");

  READ_INT("result", results);
  READ_INT("burnin", burnin);
  READ_INT("total", total);
  READ_INT("thin", thin);
  READ_INT("xsamples", xsamples);
  READ_INT("ysamples", ysamples);
  READ_INT("zsamples", zsamples);
  READ_INT("nprocesses", nprocesses);
  READ_INT("maxpartitions", maxpartitions);
  READ_DOUBLE("xmin", xmin);
  READ_DOUBLE("xmax", xmax);
  READ_DOUBLE("ymin", ymin);
  READ_DOUBLE("ymax", ymax);
  READ_DOUBLE("credibleinterval", credibleinterval);
  READ_INT("nheirarchical", nhierarchical);
  
  READ_INT("nglobalparameters", nglobal);
  global = NULL;
  if ((*nglobal) >= 0) {
    global = forwardmodelparameter_create(*nglobal);
    for (i = 0; i < (*nglobal); i ++) {
      READ_DOUBLE("min", &(global[i].fmin));
      READ_DOUBLE("max", &(global[i].fmax));
      READ_DOUBLE("std", &(global[i].fstd_value));
      READ_DOUBLE("std_bd", &(global[i].fstd_bd));
    }
  }
  *global_parameters = global;


  READ_INT("nlocalparameters", nlocal);
  local = NULL;
  if ((*nlocal) >= 0) {
    local = forwardmodelparameter_create(*nlocal);
    for (i = 0; i < (*nlocal); i ++) {
      READ_DOUBLE("min", &(local[i].fmin));
      READ_DOUBLE("max", &(local[i].fmax));
      READ_DOUBLE("std", &(local[i].fstd_value));
      READ_DOUBLE("std_bd", &(local[i].fstd_bd));
    }
  } else {
    rjmcmc_error("resultset2dfm_load: invalid number of local parameters.\n");
    return NULL;
  }
  *local_parameters = local;

  /*
   * We now have enough information to create a resultset.
   */
  rs = resultset2dfm_create(*burnin,
			    *total,
			    *thin,
			    *nglobal,
			    global,
			    *nlocal,
			    local,
			    *nhierarchical,
			    *xsamples,
			    *ysamples,
			    *zsamples,
			    *maxpartitions,
			    *xmin,
			    *xmax,
			    *ymin,
			    *ymax,
			    *nprocesses,
			    *credibleinterval,
			    *results);
  if (rs == NULL) {
    return NULL;
  }

  /*
   * Now fill in the currently stored results.
   */
  READ_INT_ARRAY("propose", rs->propose, rs->nprocesses);
  READ_INT_ARRAY("accept", rs->accept, rs->nprocesses);

  READ_DOUBLE_ARRAY("misfit", rs->misfit, rs->total);

  READ_INT_ARRAY("partitions", rs->partitions, rs->total);
  
  READ_DOUBLE_ARRAY_2D("global", rs->global, rs->nglobalparameters, rs->total);
  READ_DOUBLE_ARRAY_2D("hierarchical", rs->hierarchical, rs->nhierarchical, rs->total);

  READ_INT_ARRAY_2D("centres", rs->centres, rs->xsamples, rs->ysamples);

  if ((rs->results & (RESULTSET2DFM_MODE | RESULTSET2DFM_CREDIBLE | RESULTSET2DFM_MEDIAN)) > 0) {
    READ_INT_ARRAY_4D("hist", rs->hist, rs->nlocalparameters, rs->xsamples, rs->ysamples, rs->zsamples);
  }

  if ((rs->results & RESULTSET2DFM_MEAN) > 0) {
    READ_DOUBLE_ARRAY_3D("local_mean", rs->local_mean, rs->nlocalparameters, rs->xsamples, rs->ysamples);
    READ_DOUBLE_ARRAY_3D("local_delta", rs->local_variance, rs->nlocalparameters, rs->xsamples, rs->ysamples);
  }
  READ_INT("local_samples", &(rs->local_samples));

  if ((rs->results & RESULTSET2DFM_MODE) > 0) {
    READ_DOUBLE_ARRAY_3D("local_mode", rs->local_mode, rs->nlocalparameters, rs->xsamples, rs->ysamples);
  }

  if ((rs->results & RESULTSET2DFM_MEDIAN) > 0) {
    READ_DOUBLE_ARRAY_3D("local_median", rs->local_median, rs->nlocalparameters, rs->xsamples, rs->ysamples);
  }
  
  if ((rs->results & RESULTSET2DFM_CREDIBLE) > 0) {
    READ_DOUBLE_ARRAY_3D("local_cred_min", rs->local_cred_min, rs->nlocalparameters, rs->xsamples, rs->ysamples);
    READ_DOUBLE_ARRAY_3D("local_cred_max", rs->local_cred_max, rs->nlocalparameters, rs->xsamples, rs->ysamples);
  }

  RJMCMC_INTCHECKPTR(r->read_footer(fp), 
		     "resultset2dfm_read: failed to read footer\n");

  fclose(fp);

  return rs;
}

void 
resultset2dfm_destroy(resultset2dfm_t *r)
{
  if (r != NULL) {

    rjmcmc_destroy_int_array_1d(r->propose);
    rjmcmc_destroy_int_array_1d(r->accept);

    rjmcmc_destroy_array_1d(r->misfit);
    rjmcmc_destroy_int_array_1d(r->partitions);

    rjmcmc_destroy_array_2d(r->nglobalparameters, r->global);
    rjmcmc_destroy_array_3d(r->nlocalparameters, r->xsamples, r->local_mean);
    rjmcmc_destroy_array_3d(r->nlocalparameters, r->xsamples, r->local_variance);
    rjmcmc_destroy_array_3d(r->nlocalparameters, r->xsamples, r->local_median);
    rjmcmc_destroy_array_3d(r->nlocalparameters, r->xsamples, r->local_mode);
    rjmcmc_destroy_array_3d(r->nlocalparameters, r->xsamples, r->local_cred_min);
    rjmcmc_destroy_array_3d(r->nlocalparameters, r->xsamples, r->local_cred_max);

    rjmcmc_destroy_array_2d(r->nhierarchical, r->hierarchical);

    rjmcmc_destroy_int_array_2d(r->xsamples, r->centres);

    rjmcmc_destroy_int_array_4d(r->nlocalparameters,
				r->xsamples,
				r->ysamples,
				r->hist);

    free(r);
  }
}

void resultset2dfm_propose(resultset2dfm_t *r,
			   int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses, 
			"resultset2dfm_propose: invalid process\n");

  r->propose[p] ++;
}

void
resultset2dfm_accept(resultset2dfm_t *r,
		   int p)
{
  RJMCMC_INDEXCHECKVOID(p, 
			r->nprocesses, 
			"resultset2dfm_accept: invalid process\n");

  r->accept[p] ++;
}

void
resultset2dfm_sample_global_parameter(resultset2dfm_t *r,
				      int i,
				      int gi,
				      double v)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset2dfm_sample_global_parameter: invalid process\n");
  RJMCMC_INDEXCHECKVOID(gi, 
			r->nglobalparameters,
			"resulset2dfm_sample_global_paramaeter: invalid index\n");

  r->global[gi][i] = v;
}

void
resultset2dfm_sample_hierarchical_parameter(resultset2dfm_t *r,
					    int i,
					    int hi,
					    double v)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset2dfm_sample_hierarchical_parameter: invalid process\n");
  RJMCMC_INDEXCHECKVOID(hi, 
			r->nhierarchical,
			"resulset2dfm_sample_hierarchical_paramaeter: invalid index\n");

  r->hierarchical[hi][i] = v;
}

void
resultset2dfm_sample_local_parameter(resultset2dfm_t *r,
				     int i,
				     int li,
				     double **v)
{
  int ix;
  int iy;

  double delta;

  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset2dfm_sample_local_parameter: invalid index\n");
  RJMCMC_INDEXCHECKVOID(li, 
			r->nlocalparameters, 
			"resultset2dfm_sample_local_parameter: invalid index\n");

  if (i >= r->burnin) {

    if (r->thin == 0 ||
	(i % r->thin) == 0) {
      r->local_samples ++;
      for (ix = 0; ix < r->xsamples; ix ++) {
	for (iy = 0; iy < r->ysamples; iy ++) {
	  delta = v[ix][iy] - r->local_mean[li][ix][iy];
	  r->local_mean[li][ix][iy] += delta/(double)r->local_samples;
	  r->local_variance[li][ix][iy] += delta*(v[ix][iy] - r->local_mean[li][ix][iy]);
	}
      }
    }

    if (r->hist != NULL) {
      for (ix = 0; ix < r->xsamples; ix ++) {
	for (iy = 0; iy < r->ysamples; iy ++) {
	  r->hist[li][ix][iy][rjmcmc_map_to_index(v[ix][iy],
						  r->local_parameters[li].fmin,
						  r->local_parameters[li].fmax,
						  r->zsamples)] ++;
	}
      }
    }
  }
}

void 
resultset2dfm_sample_misfit(resultset2dfm_t *r,
			    int i,
			    double misfit)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset2dfm_sample_misfit: invalid index\n");

  r->misfit[i] = misfit;
}

void
resultset2dfm_sample_centre(resultset2dfm_t *r,
			    double x,
			    double y)
{
  int ix;
  int iy;

  ix = rjmcmc_map_to_index(x, r->xmin, r->xmax, r->xsamples);
  iy = rjmcmc_map_to_index(y, r->ymin, r->ymax, r->ysamples);
  r->centres[ix][iy] ++;
}

void
resultset2dfm_sample_npartitions(resultset2dfm_t *r,
			       int i,
			       int npartitions)
{
  RJMCMC_INDEXCHECKVOID(i, 
			r->total, 
			"resultset2dfm_sample_npartitions: invalid index\n");

  r->partitions[i] = npartitions;
}

void
resultset2dfm_assemble_results(resultset2dfm_t *r)
{
  int i;
  int j;
  int k;
  int li;
  double denom;
  int cred_samples;


  if (r->local_mean) {
    
    if (r->local_samples > 0) {
      denom = (double)r->local_samples;
    } else {
      denom = 1.0;
    }

    for (li = 0; li < r->nlocalparameters; li ++) {
      for (j = 0; j < r->xsamples; j ++) {
	for (k = 0; k < r->ysamples; k ++) {
#if 0 // Mean is now computed running so we don't need this.
	  r->local_mean[li][j][k] /= denom;
#endif 
	  r->local_variance[li][j][k] /= (denom - 1.0);
	}
      }
    }
  }

  if (r->local_median) {

    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) {
	for (k = 0; k < r->ysamples; k ++) {

	  r->local_median[i][j][k] = 
	    rjmcmc_median_from_histogram(r->hist[i][j][k],
					 r->local_parameters[i].fmin,
					 r->local_parameters[i].fmax,
					 r->zsamples);
	  
	}
      }
    }
  }

  if (r->local_mode) {
    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) {
	for (k = 0; k < r->ysamples; k ++) {

	  r->local_mode[i][j][k] = 
	    rjmcmc_mode_from_histogram(r->hist[i][j][k],
				       r->local_parameters[i].fmin,
				       r->local_parameters[i].fmax,
				       r->zsamples);
	}
      }
    }
  }

  if (r->local_cred_min &&
      r->local_cred_max) {

    cred_samples = (double)(r->total - r->burnin) * (1.0 - r->credible_interval)/2.0;
    
    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) {
	for (k = 0; k < r->ysamples; k ++) {

	  r->local_cred_min[i][j][k] = 
	    rjmcmc_head_from_histogram(r->hist[i][j][k],
				       r->local_parameters[i].fmin,
				       r->local_parameters[i].fmax,
				       r->zsamples,
				       cred_samples);
							     
	  r->local_cred_max[i][j][k] = 
	    rjmcmc_tail_from_histogram(r->hist[i][j][k],
				       r->local_parameters[i].fmin,
				       r->local_parameters[i].fmax,
				       r->zsamples,
				       cred_samples);

	}
      }
    }
  }
}

#if defined(HAVE_MPI_H)
void MPI_resultset2dfm_assemble_results(resultset2dfm_t *r,
					int mpisize,
					int mpirank,
					int root,
					MPI_Comm comm)
{
  int ivsize;
  int *iv;

  int vsize;
  double *v;

  int cred_samples;

  int i;
  int j;
  int k;

  double denom;

  vsize = r->xsamples;
  v = rjmcmc_create_array_1d(vsize);
  if (v == NULL) {
    rjmcmc_error("MPI_resultset1dfm_assemble_results: "
		 "failed to create temporary array\n");
    return;
  }

  ivsize = r->nprocesses;
  if (r->hist && ivsize < r->zsamples) {
    ivsize = r->zsamples;
  }
  iv = rjmcmc_create_int_array_1d(ivsize);
  if (iv == NULL) {
    rjmcmc_error("MPI_resultset1dfm_assemble_results: "
		 "failed to create temporary int array\n");
    return;
  }

  /*
   * Aggregate the total acceptance/proposal counts
   */
  MPI_Reduce(r->propose, iv, r->nprocesses, MPI_INT, MPI_SUM, root, comm);
  if (mpirank == root) {
    for (j = 0; j < r->nprocesses; j ++) {
      r->propose[j] = iv[j];
    }
  }

  MPI_Reduce(r->accept, iv, r->nprocesses, MPI_INT, MPI_SUM, root, comm);
  if (mpirank == root) {
    for (j = 0; j < r->nprocesses; j ++) {
      r->accept[j] = iv[j];
    }
  }

  /*
   * Aggregate the means. First compute locally then distribute.
   */
  if (r->local_mean) {

    if (r->local_samples > 0) {
      denom = (double)r->local_samples;
    } else {
      denom = 1.0;
    }

    for (i = 0; i < r->nlocalparameters; i ++) {
      for (j = 0; j < r->xsamples; j ++) { 

	for (k = 0; k < r->ysamples; k ++) {

	  r->local_variance[i][j][k] /= (denom - 1.0);

	  /*
	   * Add the mean squared so that we have just the square of the values
	   */
	  r->local_variance[i][j][k] += r->local_mean[i][j][k] * r->local_mean[i][j][k];

#if 0
	  r->local_mean[i][j][k] /= denom;
#endif 
	}

	MPI_Reduce(r->local_mean[i][j], v, r->ysamples, 
		   MPI_DOUBLE, MPI_SUM, root, comm);

	if (mpirank == root) {
	  for (k = 0; k < r->ysamples; k ++) {
	    r->local_mean[i][j][k] = v[k]/(double)mpisize;
	  }
	}

	/*
	 * Sum the variances
	 */
	MPI_Reduce(r->local_variance[i][j], v, r->ysamples,
		   MPI_DOUBLE, MPI_SUM, root, comm);

	if (mpirank == root) {
	  for (k = 0; k < r->ysamples; k ++) {
	    r->local_variance[i][j][k] = v[k]/(double)mpisize - r->local_mean[i][j][k]*r->local_mean[i][j][k];
	  }
	}

      }
    }
  }

  /*
   * Aggregate the histograms
   */
  if (r->hist) {

    cred_samples = ((double)(r->total - r->burnin) * (1.0 - r->credible_interval)/2.0) * mpisize;

    for (i = 0; i < r->nlocalparameters; i ++) {
      
      for (j = 0; j < r->xsamples; j ++) {
	for (k = 0; k < r->ysamples; k ++) {
	
	  
	  MPI_Reduce(r->hist[i][j][k], iv, r->zsamples, 
		     MPI_INT, MPI_SUM, root, comm);

	
	  
	  if (mpirank == root) {
	    
	    if (r->local_median) {
	      r->local_median[i][j][k] = 
		rjmcmc_median_from_histogram(iv,
					     r->local_parameters[i].fmin,
					     r->local_parameters[i].fmax,
					     r->zsamples);
	    }
	    
	    if (r->local_mode) {
	      r->local_median[i][j][k] = 
		rjmcmc_mode_from_histogram(iv,
					   r->local_parameters[i].fmin,
					   r->local_parameters[i].fmax,
					   r->zsamples);
	    }
	    
	    if (r->local_cred_min &&
		r->local_cred_max) {
	      
	      r->local_cred_min[i][j][k] = 
		rjmcmc_head_from_histogram(iv,
					   r->local_parameters[i].fmin,
					   r->local_parameters[i].fmax,
					   r->zsamples,
					   cred_samples);
	      
	      r->local_cred_max[i][j][k] = 
		rjmcmc_tail_from_histogram(iv,
					   r->local_parameters[i].fmin,
					   r->local_parameters[i].fmax,
					   r->zsamples,
					   cred_samples);
	      
	    }
	  }
	}
      }
    }
  }      
  
  rjmcmc_destroy_int_array_1d(iv);
  rjmcmc_destroy_array_1d(v);
}

#endif /* HAVE_MPI_H */

/*
 * Getting results
 */

int
resultset2dfm_get_nparameters(resultset2dfm_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->nprocesses;
}

int 
resultset2dfm_get_total(resultset2dfm_t *r)
{
  if (r == NULL) {
    return -1;
  }

  return r->total;
}

const int *
resultset2dfm_get_propose(resultset2dfm_t *r,
			  int *nprocesses)
{
  if (r == NULL) {
    return NULL;
  }

  if (nprocesses != NULL) {
    *nprocesses = r->nprocesses;
  }

  return r->propose;
}

int
resultset2dfm_get_propose_f(resultset2dfm_t *r,
			    int maxsize,
			    int *propose)
{
  int i;
  int n;

  if (r == NULL) {
    return -1;
  }

  n = resultset2dfm_get_nparameters(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    propose[i] = r->propose[i];
  }

  return n;
}

const int *
resultset2dfm_get_accept(resultset2dfm_t *r,
			 int *nprocesses)
{
  if (r == NULL) {
    return NULL;
  }

  if (nprocesses != NULL) {
    *nprocesses = r->nprocesses;
  }

  return r->accept;
}

int
resultset2dfm_get_accept_f(resultset2dfm_t *r,
			    int maxsize,
			    int *accept)
{
  int i;
  int n;

  if (r == NULL) {
    return -1;
  }

  n = resultset2dfm_get_nparameters(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    accept[i] = r->accept[i];
  }

  return n;
}

const double *
resultset2dfm_get_misfit(resultset2dfm_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->misfit;
}

int
resultset2dfm_get_misfit_f(resultset2dfm_t *r,
			   int maxsize,
			   double *misfit)
{
  int i;
  int n;

  if (r == NULL) {
    return -1;
  }

  n = resultset2dfm_get_total(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    misfit[i] = r->misfit[i];
  }

  return n;
}

const int *
resultset2dfm_get_partitions(resultset2dfm_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return r->partitions;
}

int
resultset2dfm_get_partitions_f(resultset2dfm_t *r,
			       int maxsize,
			       int *npartitions)
{
  int i;
  int n;

  if (r == NULL) {
    return -1;
  }

  n = resultset2dfm_get_total(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    npartitions[i] = r->partitions[i];
  }

  return n;
}

const double *
resultset2dfm_get_global_parameter(resultset2dfm_t *r, int gi)
{
  if (r == NULL) {
    return NULL;
  }

  return r->global[gi];
}

int
resultset2dfm_get_global_parameter_f(resultset2dfm_t *r,
				     int gi,
				     int maxsize,
				     double *global_parameter)
{
  int i;
  int n;

  if (r == NULL) {
    return -1;
  }

  if (gi < 0 || gi >= r->nglobalparameters) {
    return -1;
  }

  n = resultset2dfm_get_total(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    global_parameter[i] = r->global[gi][i];
  }

  return n;
}


const double *
resultset2dfm_get_hierarchical_parameter(resultset2dfm_t *r, int gi)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset2dfm_get_hierarchical_parameters: " 
		      "null results\n");
  RJMCMC_INDEXCHECKPTR(gi, r->nhierarchical,
		      "resultset2dfm_get_hierarchical_parameters: " 
		      "invalid index\n");
		      
  return r->hierarchical[gi];
}

int
resultset2dfm_get_hierarchical_parameter_f(resultset2dfm_t *r,
					   int hi,
					   int maxsize,
					   double *hierarchical_parameter)
{
  int i;
  int n;

  if (r == NULL) {
    return -1;
  }

  if (hi < 0 || hi >= r->nhierarchical) {
    return -1;
  }

  n = resultset2dfm_get_total(r);
  if (maxsize < n) {
    n = maxsize;
  }

  for (i = 0; i < n; i ++) {
    hierarchical_parameter[i] = r->hierarchical[hi][i];
  }

  return n;
}

const double **
resultset2dfm_get_local_parameter_mean(resultset2dfm_t *r, int li)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset2dfm_get_local_parameter_mean: "
		      "null results\n");
  RJMCMC_INDEXCHECKPTR(li,
		       r->nlocalparameters,
		      "resultset2dfm_get_local_parameter_mean: "
		      "invalid index\n");
  RJMCMC_NULLCHECKPTR(r->local_mean, 
		      "resultset2dfm_get_local_parameter_mean: "
		      "null mean\n");
		       
  return (const double **)r->local_mean[li];
}

int
resultset2dfm_get_local_parameter_mean_f(resultset2dfm_t *r,
					 int li,
					 int xsamples,
					 int ysamples,
					 double *mean)
{
  const double **lmean;
  int x;
  int y;

  lmean = resultset2dfm_get_local_parameter_mean(r, li);
  if (lmean == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset2dfm_get_local_parameter_mean_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_local_parameter_mean_f: "
			   "invalid ysamples\n");

  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      mean[x + y * xsamples] = lmean[x][y];

    }
  }
  
  return 0;
}

const double **
resultset2dfm_get_local_parameter_variance(resultset2dfm_t *r, int li)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset2dfm_get_local_parameter_variance: "
		      "null results\n");
  RJMCMC_INDEXCHECKPTR(li,
		       r->nlocalparameters,
		      "resultset2dfm_get_local_parameter_variance: "
		      "invalid index\n");
  RJMCMC_NULLCHECKPTR(r->local_variance, 
		      "resultset2dfm_get_local_parameter_variance: "
		      "null variance\n");
		       
  return (const double **)r->local_variance[li];
}


int
resultset2dfm_get_local_parameter_variance_f(resultset2dfm_t *r,
					     int li,
					     int xsamples,
					     int ysamples,
					     double *variance)
{
  const double **lvariance;
  int x;
  int y;

  lvariance = resultset2dfm_get_local_parameter_variance(r, li);
  if (lvariance == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset2dfm_get_local_parameter_variance_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_local_parameter_variance_f: "
			   "invalid ysamples\n");

  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      variance[x + y * xsamples] = lvariance[x][y];

    }
  }
  
  return 0;
}


const double **
resultset2dfm_get_local_parameter_mode(resultset2dfm_t *r, int li)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset2dfm_get_local_parameter_mode: "
		      "null results\n");
  RJMCMC_INDEXCHECKPTR(li,
		       r->nlocalparameters,
		      "resultset2dfm_get_local_parameter_mode: "
		      "invalid index\n");
  RJMCMC_NULLCHECKPTR(r->local_mode, 
		      "resultset2dfm_get_local_parameter_mode: "
		      "null mode\n");
		       
  return (const double **)r->local_mode[li];
}

int
resultset2dfm_get_local_parameter_mode_f(resultset2dfm_t *r,
					 int li,
					 int xsamples,
					 int ysamples,
					 double *mode)
{
  const double **lmode;
  int x;
  int y;

  lmode = resultset2dfm_get_local_parameter_mode(r, li);
  if (lmode == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset2dfm_get_local_parameter_mode_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_local_parameter_mode_f: "
			   "invalid ysamples\n");

  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      mode[x + y * xsamples] = lmode[x][y];

    }
  }
  
  return 0;
}

const double **
resultset2dfm_get_local_parameter_median(resultset2dfm_t *r, int li)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset2dfm_get_local_parameter_median: "
		      "null results\n");
  RJMCMC_INDEXCHECKPTR(li,
		       r->nlocalparameters,
		      "resultset2dfm_get_local_parameter_median: "
		      "invalid index\n");
  RJMCMC_NULLCHECKPTR(r->local_median, 
		      "resultset2dfm_get_local_parameter_median: "
		      "null median\n");
  
		       
  return (const double **)r->local_median[li];
}

int
resultset2dfm_get_local_parameter_median_f(resultset2dfm_t *r,
					 int li,
					 int xsamples,
					 int ysamples,
					 double *median)
{
  const double **lmedian;
  int x;
  int y;

  lmedian = resultset2dfm_get_local_parameter_median(r, li);
  if (lmedian == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset2dfm_get_local_parameter_median_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_local_parameter_median_f: "
			   "invalid ysamples\n");

  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      median[x + y * xsamples] = lmedian[x][y];

    }
  }
  
  return 0;
}

const double **
resultset2dfm_get_local_parameter_credible_min(resultset2dfm_t *r, int li)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset2dfm_get_local_parameter_credible_min: "
		      "null results\n");
  RJMCMC_INDEXCHECKPTR(li,
		       r->nlocalparameters,
		      "resultset2dfm_get_local_parameter_credible_min: "
		      "invalid index\n");
  RJMCMC_NULLCHECKPTR(r->local_cred_min,
		      "resultset2dfm_get_local_parameter_credible_min: "
		      "null credible min\n");
		       
  return (const double **)r->local_cred_min[li];
}

int
resultset2dfm_get_local_parameter_credible_min_f(resultset2dfm_t *r,
						 int li,
						 int xsamples,
						 int ysamples,
						 double *credible_min)
{
  const double **lcredible_min;
  int x;
  int y;

  lcredible_min = resultset2dfm_get_local_parameter_credible_min(r, li);
  if (lcredible_min == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset2dfm_get_local_parameter_credible_min_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_local_parameter_credible_min_f: "
			   "invalid ysamples\n");

  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      credible_min[x + y * xsamples] = lcredible_min[x][y];

    }
  }
  
  return 0;
}


const double ** resultset2dfm_get_local_parameter_credible_max(resultset2dfm_t *r, int li)
{
  RJMCMC_NULLCHECKPTR(r, 
		      "resultset2dfm_get_local_parameter_credible_max: "
		      "null results\n");
  RJMCMC_INDEXCHECKPTR(li,
		       r->nlocalparameters,
		      "resultset2dfm_get_local_parameter_credible_max: "
		      "invalid index\n");
  RJMCMC_NULLCHECKPTR(r->local_cred_max,
		      "resultset2dfm_get_local_parameter_credible_min: "
		      "null credible max\n");
		       
  return (const double **)r->local_cred_max[li];
}

int
resultset2dfm_get_local_parameter_credible_max_f(resultset2dfm_t *r,
						 int li,
						 int xsamples,
						 int ysamples,
						 double *credible_max)
{
  const double **lcredible_max;
  int x;
  int y;

  lcredible_max = resultset2dfm_get_local_parameter_credible_max(r, li);
  if (lcredible_max == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset2dfm_get_local_parameter_credible_max_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_local_parameter_credible_max_f: "
			   "invalid ysamples\n");

  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      credible_max[x + y * xsamples] = lcredible_max[x][y];

    }
  }
  
  return 0;
}

const int **
resultset2dfm_get_centres(resultset2dfm_t *r)
{
  if (r == NULL) {
    return NULL;
  }

  return (const int**)r->centres;
}

int
resultset2dfm_get_centres_f(resultset2dfm_t *r,
			    int li,
			    int xsamples,
			    int ysamples,
			    int *centres)
{
  const int **lcentres;
  int x;
  int y;

  lcentres = resultset2dfm_get_centres(r);
  if (lcentres == NULL) {
    return -1;
  }

  RJMCMC_CONDITIONCHECKINT(xsamples != r->xsamples,
			   "resultset2dfm_get_centres_f: "
			   "invalid xsamples\n");
  RJMCMC_CONDITIONCHECKINT(ysamples != r->ysamples,
			   "resultset2dfm_get_centres_f: "
			   "invalid ysamples\n");

  for (y = 0; y < ysamples; y ++) {
    for (x = 0; x < xsamples; x ++) {

      centres[x * xsamples + y] = lcentres[x][y];

    }
  }
  
  return 0;
}

void
resultset2dfm_fill_xcoord_vector(resultset2dfm_t *r,
				 double *x,
				 int *l)
{
  int i;
  int e;

  if (r == NULL) {
    return;
  }

  e = *l;
  if (e > r->xsamples) {
    e = r->xsamples;
    *l = e;
  }

  rjmcmc_fill_coord_vector(r->xmin, r->xmax, e, x);
}

int
resultset2dfm_fill_xcoord_vector_f(resultset2dfm_t *r,
				   int maxsize,
				   double *x)
{
  int s;

  s = maxsize;
  resultset2dfm_fill_xcoord_vector(r, x, &s);

  return s;
}

void
resultset2dfm_fill_ycoord_vector(resultset2dfm_t *r,
				 double *y,
				 int *l)
{
  int i;
  int e;

  if (r == NULL) {
    return;
  }

  e = *l;
  if (e > r->ysamples) {
    e = r->ysamples;
    *l = e;
  }

  rjmcmc_fill_coord_vector(r->ymin, r->ymax, e, y);
}

int
resultset2dfm_fill_ycoord_vector_f(resultset2dfm_t *r,
				   int maxsize,
				   double *y)
{
  int s;

  s = maxsize;
  resultset2dfm_fill_ycoord_vector(r, y, &s);

  return s;
}

static int binary_write_header(FILE *fp)
{
  return 0;
}

static int binary_write_variable_header(FILE *fp, const char *name)
{
  return 0;
}

static int binary_write_variable_footer(FILE *fp)
{
  return 0;
}

static int binary_write_int(FILE *fp, int i)
{
  return binary_write_int_array(fp, &i, 1);
}

static int binary_write_int_array(FILE *fp, int *a, int size)
{
  if (fwrite(a, sizeof(int), size, fp) != size) {
    return -1;
  }

  return 0;
}

static int binary_write_int_array_2d(FILE *fp, int **a, int width, int height)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_write_int_array(fp, a[i], height) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_write_int_array_3d(FILE *fp, int ***a, int width, int height, int depth)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_write_int_array_2d(fp, a[i], height, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_write_int_array_4d(FILE *fp, int ****a, int width, int height, int depth, int length)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_write_int_array_3d(fp, a[i], height, depth, length) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_write_double(FILE *fp, double d)
{
  return binary_write_double_array(fp, &d, 1);
}

static int binary_write_double_array(FILE *fp, double *a, int size)
{
  if (fwrite(a, sizeof(double), size, fp) != size) {
    return -1;
  }

  return 0;
}

static int binary_write_double_array_2d(FILE *fp, double **a, int width, int height)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_write_double_array(fp, a[i], height) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_write_double_array_3d(FILE *fp, double ***a, int width, int height, int depth)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_write_double_array_2d(fp, a[i], height, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_write_footer(FILE *fp)
{
  return 0;
}

static int binary_read_header(FILE *fp)
{
  return 0;
}

static int binary_read_variable_header(FILE *fp, const char *name)
{
  return 0;
}

static int binary_read_variable_footer(FILE *fp)
{ 
  return 0;
}

static int binary_read_int(FILE *fp, int *i)
{
  return binary_read_int_array(fp, i, 1);
}

static int binary_read_int_array(FILE *fp, int *a, int size)
{
  if (fread(a, sizeof(int), size, fp) != size) {
    return -1;
  }

  return 0;
}

static int binary_read_int_array_2d(FILE *fp, int **a, int width, int height)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_read_int_array(fp, a[i], height) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_read_int_array_3d(FILE *fp, int ***a, int width, int height, int depth)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_read_int_array_2d(fp, a[i], height, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_read_int_array_4d(FILE *fp, int ****a, int width, int height, int depth, int length)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_read_int_array_3d(fp, a[i], height, depth, length) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_read_double(FILE *fp, double *d)
{
  return binary_read_double_array(fp, d, 1);
}

static int binary_read_double_array(FILE *fp, double *a, int size)
{
  if (fread(a, sizeof(double), size, fp) != size) {
    return -1;
  }

  return 0;
}

static int binary_read_double_array_2d(FILE *fp, double **a, int width, int height)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_read_double_array(fp, a[i], height) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_read_double_array_3d(FILE *fp, double ***a, int width, int height, int depth)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (binary_read_double_array_2d(fp, a[i], height, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

static int binary_read_footer(FILE *fp)
{
  return 0;
}
