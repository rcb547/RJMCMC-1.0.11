
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultsetfm.h>

#include <rjmcmc/rjmcmc_util.h>

static double unif(void)
{
  return rjmcmc_uniform();
}

static double norm(void)
{
  return rjmcmc_normal();
}

typedef enum {
  P_XC = 0,
  P_YC,
  P_TOTAL
} parameter_t;

static double gaussian_misfit(void *user_arg,
			      int nvalues,
			      const double *values);

FILE *path = NULL;

double *path_x;
double *path_y;
int path_i;

#define TOTAL 20000

#define HIST_XSAMPLES 50
#define HIST_YSAMPLES 50
#define HIST_XMIN (-10.0)
#define HIST_XMAX (10.0)
#define HIST_YMIN (-10.0)
#define HIST_YMAX (10.0)

int hist2d[HIST_XSAMPLES][HIST_YSAMPLES];
int hsamples;
int hthin;

static void histogram_init(void);
static void histogram_addpoint(double x, double y);
static void histogram_save(const char *filename);

static void histogram_save_as_geo(const char *filename, int samples);
static void path_save_as_geo(const char *filename, int total);

int main(int argc, char *argv[]) 
{
  int burnin = 1000;
  int total = TOTAL;
  int samples = 100;
  double confidence_interval = 0.95;
  int requested_results = RESULTSETFM_MEAN;

  forwardmodelparameter_t parameters[P_TOTAL];

  int nproc;
  const double *v;
  const int *iv;
  int i;

  double *xcoords;
  double *y;
  int xcl;

  resultsetfm_t *results;

  histogram_init();
  hsamples = 0;
  hthin = 10;

  path = fopen("path.txt", "w");
  if (path == NULL) {
    fprintf(stderr, "error: failed to create path file\n");
    return -1;
  }

  path_x = rjmcmc_create_array_1d(total);
  path_y = rjmcmc_create_array_1d(total);
  if (path_x == NULL || path_y == NULL) {
    fprintf(stderr, "error: failed to create path vectors\n");
    return -1;
  }
  path_i = 0;

  /*
   * Initialize the search space for the parameters
   */
  
  parameters[P_XC].fmin = -10.0;
  parameters[P_XC].fmax = 10.0;
  parameters[P_XC].fstd_value = 1.0;
  parameters[P_XC].fstd_bd = 0.0;

  parameters[P_YC].fmin = -10.0;
  parameters[P_YC].fmax = 10.0;
  parameters[P_YC].fstd_value = 1.0;
  parameters[P_YC].fstd_bd = 0.0;

  /*
   * Run the forward model
   */
  results = single_forwardmodel(burnin,
				total,
				unif, /*rjmcmc_uniform,*/
				norm, /*rjmcmc_normal,*/
				P_TOTAL,
				parameters,
				gaussian_misfit,
				NULL,
				samples,
				confidence_interval,
				requested_results);

  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run functionfit\n");
    return -1;
  }

  v = resultsetfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("functionfit.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultsetfm_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultsetfm_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  printf("x mean: %f\n", 
	 resultsetfm_get_parameter_mean(results, P_XC));
  printf("y mean: %f\n", 
	 resultsetfm_get_parameter_mean(results, P_YC));

  resultsetfm_destroy(results);

  fclose(path);

  rjmcmc_destroy_array_1d(path_x);
  rjmcmc_destroy_array_1d(path_y);

  histogram_save("final_histogram.txt");

  histogram_save_as_geo("final_histogram.geo", total);

  return 0;
}

static double gaussian_misfit(void *user_arg,
			      int nvalues,
			      const double *values)
{
  double xc = values[P_XC];
  double yc = values[P_YC];

  histogram_addpoint(xc, yc);
  hsamples ++;

  fprintf(path, "%f %f\n", xc, yc);
  path_x[path_i] = xc;
  path_y[path_i] = yc;
  path_i ++;

  if (hsamples % hthin == 0) {
    char hfilename[256];

    sprintf(hfilename, "model/hist2d_%d.geo", hsamples/hthin);
    histogram_save_as_geo(hfilename, hsamples);

    sprintf(hfilename, "model/path_%d.geo", hsamples/hthin);
    path_save_as_geo(hfilename, TOTAL);
  }
  

  return (xc*xc + yc*yc)/2.0;
}


static void histogram_init(void)
{
  int i;
  int j;
  
  for (i = 0; i < HIST_XSAMPLES; i ++) {
    for (j = 0; j < HIST_YSAMPLES; j ++) {
      hist2d[i][j] = 0;
    }
  }
}

static void histogram_addpoint(double x, double y)
{
  int i;
  int j;
  
  i = (x - HIST_XMIN)/(HIST_XMAX - HIST_XMIN) * (double)(HIST_XSAMPLES);
  if (i < 0) { 
    i = 0;
  }
  if (i >= HIST_XSAMPLES) {
    i = HIST_XSAMPLES - 1;
  }

  
  j = (y - HIST_YMIN)/(HIST_YMAX - HIST_YMIN) * (double)(HIST_YSAMPLES);
  if (j < 0) { 
    j = 0;
  }
  if (j >= HIST_YSAMPLES) {
    j = HIST_YSAMPLES - 1;
  }

  hist2d[i][j] ++;
}

static void histogram_save(const char *filename)
{
  FILE *fp;
  int i;
  int j;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "histogram_save: failed to create file\n");
    return;
  }

  for (j = 0; j < HIST_YSAMPLES; j ++) {
    for (i = 0; i < HIST_XSAMPLES; i ++) {
      fprintf(fp, "%d ", hist2d[i][j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

static void histogram_save_as_geo(const char *filename, int samples)
{

  FILE *fp;

  int xi;
  int yi;

  double xw;
  double yw;

  double xc;
  double yc;
  double zh;

  double *px;
  double *py;
  double *pz;

  int nboxes;
  int npoints;
  int nfaces;

  int offset;

  int i;

  xw = 0.8 * (HIST_XMAX - HIST_XMIN)/(2.0 * (double)(HIST_XSAMPLES));
  yw = 0.8 * (HIST_YMAX - HIST_YMIN)/(2.0 * (double)(HIST_YSAMPLES));
    
  if (samples <= 0) {
    samples = 1;
  }

  nboxes = HIST_XSAMPLES * HIST_YSAMPLES;
  npoints = nboxes * 8;
  nfaces = nboxes * 6;
  px = malloc(sizeof(double) * npoints);
  py = malloc(sizeof(double) * npoints);
  pz = malloc(sizeof(double) * npoints);
  if (px == NULL ||
      py == NULL ||
      pz == NULL) {
    fprintf(stderr, "histogram_save_as_geo: failed to allocate memory\n");
    return;
  }

  for (i = 0; i < npoints; i ++) {
    px[i] = -100.0;
    py[i] = -100.0;
    pz[i] = -100.0;
  }
    

  offset = 0;
  for (xi = 0; xi < HIST_XSAMPLES; xi ++) {
  
    xc = ((double)(xi) + 0.5)/(double)(HIST_XSAMPLES) * (HIST_XMAX - HIST_XMIN) + HIST_XMIN;

    for (yi = 0; yi < HIST_YSAMPLES; yi ++, offset += 8) {
      
      yc = ((double)(yi) + 0.5)/(double)(HIST_YSAMPLES) * (HIST_YMAX - HIST_YMIN) + HIST_YMIN;
      zh = (double)(hist2d[xi][yi])/(double)(samples);
      
      px[offset + 0] = xc - xw;
      py[offset + 0] = yc - yw;
      pz[offset + 0] = 0.0;

      px[offset + 1] = xc + xw;
      py[offset + 1] = yc - yw;
      pz[offset + 1] = 0.0;

      px[offset + 2] = xc + xw;
      py[offset + 2] = yc + yw;
      pz[offset + 2] = 0.0;
      
      px[offset + 3] = xc - xw;
      py[offset + 3] = yc + yw;
      pz[offset + 3] = 0.0;
      
      px[offset + 4] = xc - xw;
      py[offset + 4] = yc - yw;
      pz[offset + 4] = zh;

      px[offset + 5] = xc + xw;
      py[offset + 5] = yc - yw;
      pz[offset + 5] = zh;

      px[offset + 6] = xc + xw;
      py[offset + 6] = yc + yw;
      pz[offset + 6] = zh;
      
      px[offset + 7] = xc - xw;
      py[offset + 7] = yc + yw;
      pz[offset + 7] = zh;
      
    }
  }

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "histogram_save_as_geo: failed to create file\n");
    return;
  }

  fprintf(fp, "PGEOMETRY V5\n");
  fprintf(fp, "NPoints %d NPrims %d\n", npoints, nfaces);
  fprintf(fp, "NPointGroups 0 NPrimGroups 0\n");
  fprintf(fp, "NPointAttrib 0 NVertexAttrib 0 NPrimAttrib 0 NAttrib 0\n");
  
  for (i = 0; i < npoints; i ++) {
    fprintf(fp, "%f %f %f 1.0\n", px[i], py[i], pz[i]);
  }

  offset = 0;
  for (i = 0; i < nboxes; i ++, offset += 8) {
    fprintf(fp, "Run 6 Poly\n");
    fprintf(fp, " 4 < %d %d %d %d\n", offset + 0, offset + 1, offset + 2, offset + 3);
    fprintf(fp, " 4 < %d %d %d %d\n", offset + 0, offset + 1, offset + 5, offset + 4);
    fprintf(fp, " 4 < %d %d %d %d\n", offset + 1, offset + 2, offset + 6, offset + 5);
    fprintf(fp, " 4 < %d %d %d %d\n", offset + 2, offset + 3, offset + 7, offset + 6);
    fprintf(fp, " 4 < %d %d %d %d\n", offset + 3, offset + 0, offset + 4, offset + 7);
    fprintf(fp, " 4 < %d %d %d %d\n", offset + 4, offset + 5, offset + 6, offset + 7);

  }

  fprintf(fp, "beginExtra\n");
  fprintf(fp, "endExtra\n");
 
  free(px);
  free(py);
  free(pz);

  fclose(fp);
}

static void path_save_as_geo(const char *filename, int total)
{
  FILE *fp;

  int npoints;

  int i;

  npoints = path_i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "path_save_as_geo: failed to create file\n");
    return;
  }

  fprintf(fp, "PGEOMETRY V5\n");
  fprintf(fp, "NPoints %d NPrims %d\n", npoints, npoints - 1);
  fprintf(fp, "NPointGroups 0 NPrimGroups 0\n");
  fprintf(fp, "NPointAttrib 0 NVertexAttrib 0 NPrimAttrib 0 NAttrib 0\n");

  for (i = 0; i < npoints; i ++) {
    fprintf(fp, "%f %f %f 1.0\n", path_x[i], path_y[i], (double)(npoints - 1 - i)/(double)total * 10.0);
  }

  for (i = 1; i < npoints; i ++) {
    fprintf(fp, "Poly 2 : %d %d\n", i - 1, i);
  }

  fprintf(fp, "beginExtra\n");
  fprintf(fp, "endExtra\n");

  fclose(fp);
}
