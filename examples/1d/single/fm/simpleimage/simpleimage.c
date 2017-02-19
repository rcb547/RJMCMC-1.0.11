/* simpleimage

   This example segments/locates a gaussian shape within a noise image. It
   does this by using a model for the image based on several parameters
   and perturbing these using a MCMC scheme and convergence is achieved
   by implementing a misfit function that compares the model to the 
   actual image (a synthetic test image created by mkimage.c).

   The model is constructed using the formula:

   I(x,y) = B + n + A*exp(ga*(x - x0)^2 + 
                          2.0*gb*(px - x)*(py - y) + 
                          gc*(y - y0)^2)

   Where B = The background intensity
         n = Noise (gaussian)
         A = Shape peak intensity
         x0 = Horizontal centre (centre of image is 0,0)
         y0 = Vertical centre
         ga, gb, gc = Coefficients of gaussian derived from the remaining
	 model parameters
         sigma_x = x std. deviation of the gaussian
         sigma_y = y std. deviation of the gaussian
         theta = angle of the gaussian

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_util.h>

#include "image.h"

/* We use an enum to index our model parameters */
enum {
  P_A = 0,
  P_sigmax,
  P_sigmay,
  P_theta,
  P_x0,
  P_y0,
  P_B,
  P_TOTAL
};

/* The corresponding names for our parameters (for nicely printing out the
   accept/reject ratios */
static const char *P_NAMES[] = {
  "amplitude",
  "sigmax",
  "sigmay",
  "theta",
  "x0",
  "y0", 
  "background",
  "noise"
};

/* This is the data structure used to hold the image data */
struct image {
  int width;
  int height;
  double sigma;
  double *data;
};

/* Loads the PGM image with a fixed sigma. When running the hierarchical
   version, the sigma value is stored but not used in the MCMC process.
*/
struct image *image_load(const char *filename, double sigma);

/* The misfit function for the normal forward model */
double image_likelihood(void *user_arg, 
			int n,
			const double *values);

/* The misfit function for the hierarchical forward model */
double image_likelihood_hierarchical(void *user_arg, 
				     int n,
				     const double *values,
				     int hierarchical,
				     int nhierarchical,
				     const double *hierarchical_values,
				     double *logdetce);

double sqr(double x)
{
  return x*x;
}

int main(int argc, char *argv[])
{
  int burnin;
  int total;
  int thin;

  int nparameters;
  forwardmodelparameter_t *parameters;

  int nhierarchicalparameters;
  forwardmodelparameter_t *hierarchicalparameters;
  
  struct image *img;

  resultsetfm_t *results;

  double sigma;

  const double *v;
  const int *propose;
  const int *accept;
  int nprocess;

  int i;
  char filename[256];

  int hierarchical;

  int samples;
  double confidence_interval;
  int requested_results;

  double mean[P_TOTAL];

  /* Set the number of iterations etc */
  burnin = 2000;
  total = 20000;

  /* Create the forward model parameters array */
  nparameters = P_TOTAL;
  parameters = forwardmodelparameter_create(nparameters);

  /* The Amplitude parameter */
  parameters[P_A].fmin = 0.0;
  parameters[P_A].fmax = 1.0;
  parameters[P_A].fstd_value = 0.05;
  parameters[P_A].fstd_bd = 0.0;

  /* The sigma_x parameter, essentially the radius in the x direction */
  parameters[P_sigmax].fmin = 0.001;
  parameters[P_sigmax].fmax = 0.5;
  parameters[P_sigmax].fstd_value = 0.01;
  parameters[P_sigmax].fstd_bd = 0.0;
  
  /* The sigma_y parameter, essentially the radius in the y direction */
  parameters[P_sigmay].fmin = 0.001;
  parameters[P_sigmay].fmax = 0.5;
  parameters[P_sigmay].fstd_value = 0.01;
  parameters[P_sigmay].fstd_bd = 0.0;

  /* The theta parameter. Note that due to symmetry we only need to check
     the range 0 - 90 degrees. */
  parameters[P_theta].fmin = 0.0;
  parameters[P_theta].fmax = M_PI/2.0;
  parameters[P_theta].fstd_value = M_PI/36.0;
  parameters[P_theta].fstd_bd = 0.0;

  /* The x centre of the gaussian */
  parameters[P_x0].fmin = -0.5;
  parameters[P_x0].fmax = 0.5;
  parameters[P_x0].fstd_value = 0.003;
  parameters[P_x0].fstd_bd = 0.0;
  
  /* The y centre of the gaussian */
  parameters[P_y0].fmin = -0.5;
  parameters[P_y0].fmax = 0.5;
  parameters[P_y0].fstd_value = 0.003;
  parameters[P_y0].fstd_bd = 0.0;

  /* The background level */
  parameters[P_B].fmin = 0.0;
  parameters[P_B].fmax = 0.5;
  parameters[P_B].fstd_value = 0.003;
  parameters[P_B].fstd_bd = 0.0;

  /* Set this value to 0 for a simple forward model or 1 to estimate the
     noise level */
  hierarchical = 1;

  /* These values represent the prior on the noise for hierarchical 
     simulations */
  nhierarchicalparameters = 1;
  hierarchicalparameters = forwardmodelparameter_create(nhierarchicalparameters);
  hierarchicalparameters[0].fmin = 0.01;
  hierarchicalparameters[0].fmax = 0.3;
  hierarchicalparameters[0].fstd_value = 0.001;
  hierarchicalparameters[0].fstd_bd = 0.0;

  /* The fixed sigma value for simple forward model simulations */
  sigma = 0.15;

  img = image_load("image.pgm", sigma);
  if (img == NULL) {
    return -1;
  }

  /* Set the type/quality of the results */
  samples = 100;
  confidence_interval = 0.95;
  requested_results = RESULTSETFM_MEAN;

  if (hierarchical) {
    results = single_forwardmodel_hierarchical(burnin,
					       total,
					       rjmcmc_uniform,
					       rjmcmc_normal,
					       nparameters,
					       parameters,
					       nhierarchicalparameters,
					       hierarchicalparameters,
					       image_likelihood_hierarchical,
					       img,
					       samples,
					       confidence_interval,
					       requested_results);
  } else {
    results = single_forwardmodel(burnin,
				  total,
				  rjmcmc_uniform,
				  rjmcmc_normal,
				  nparameters,
				  parameters,
				  image_likelihood,
				  img,
				  samples,
				  confidence_interval,
				  requested_results);
  }
  if (results == NULL) {
    fprintf(stderr, "error: failed to run simulation\n");
    return -1;
  }

  v = resultsetfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to misfit history\n");
    return -1;
  }
  if (rjmcmc_save_vector("simpleimage.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit history\n");
    return -1;
  }

  for (i = 0; i < nparameters; i ++) {
    v = resultsetfm_get_parameter_history(results, i);
    if (v == NULL) {
      fprintf(stderr, "error: failed to get parameter %d results\n", i);
      return -1;
    }

    sprintf(filename, "simpleimage.%s", P_NAMES[i]);
    if (rjmcmc_save_vector(filename, v, total) < 0) {
      fprintf(stderr, "error: failed to save parameter %d\n", i);
      return -1;
    }

    mean[i] = rjmcmc_mean_skip(v, burnin, total);
  }

  if (hierarchical) {
    v = resultsetfm_get_hierarchical_parameter_history(results, 0);
    if (v == NULL) {
      fprintf(stderr, "error: failed to get sigma history\n");
      return -1;
    }
    
    if (rjmcmc_save_vector("simpleimage.sigma", v, total) < 0) {
      fprintf(stderr, "error: failed to save sigma history\n");
      return -1;
    }
  }

  if (hierarchical) {
    for (i = 0; i < (P_TOTAL + 1); i ++) {
      printf("%12s ", P_NAMES[i]);
    }
  } else {
    for (i = 0; i < P_TOTAL; i ++) {
      printf("%12s ", P_NAMES[i + 1]);
    }
  }    
  printf("\n");

  propose = resultsetfm_get_propose(results, &nprocess);
  if (propose == NULL) {
    fprintf(stderr, "error: failed to get the proposed counts\n");
    return -1;
  }
  for (i = 0; i < nprocess; i ++) {
    printf("%12d ", propose[i]);
  }
  printf("\n");

  accept = resultsetfm_get_accept(results, &nprocess);
  if (accept == NULL) {
    fprintf(stderr, "error: failed to get the accepted counts\n");
    return -1;
  }
  for (i = 0; i < nprocess; i ++) {
    printf("%12d ", accept[i]);
  }
  printf("\n");
  for (i = 0; i < nprocess; i ++) {
    printf("%11.3f%% ", (float)accept[i]/(float)propose[i] * 100.0);
  }
  printf("\n");

  printf("\n");
  for (i = 0; i < P_TOTAL; i ++) {
    printf("%12s (mean): %f\n", P_NAMES[i], mean[i]);
  }

  if (image_createandwritepgm("image_mean.pgm",
			      img->width,
			      img->height,
			      mean[P_A],
			      mean[P_sigmax],
			      mean[P_sigmay],
			      mean[P_theta],
			      mean[P_x0],
			      mean[P_y0],
			      mean[P_B],
			      0.0) < 0) {
    fprintf(stderr, "error: failed to save mean image\n");
    return -1;
  }

  resultsetfm_destroy(results);
  forwardmodelparameter_destroy(parameters);
  free(img->data);
  free(img);

  return 0;
}

struct image *image_load(const char *filename, double sigma)
{
  FILE *fp;
  char buffer[256];
  int width;
  int height;
  int depth;
  int i;
  int j;
  int c;
  double *p;

  struct image *r;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, 
	    "error: failed to load %s, have you run mkimage?\n", 
	    filename);
    return NULL;
  }

  if (fgets(buffer, sizeof(buffer), fp) == NULL) {
    fprintf(stderr, 
	    "error: failed to read header\n");
    return NULL;
  }

  if (strncmp("P5", buffer, 2) != 0) {
    fprintf(stderr, 
	    "error: invalid header\n");
    return NULL;
  }

  if (fscanf(fp, "%d %d\n%d\n", &width, &height, &depth) != 3) {
    fprintf(stderr,
	    "error: failed to read image dimensions\n");
    return NULL;
  }

  if (depth != 255) {
    fprintf(stderr,
	    "error: invalid depth\n");
    return NULL;
  }

  r = (struct image*)malloc(sizeof(struct image));
  if (r == NULL) {
    fprintf(stderr, "error: failed to allocate memory for image\n");
    return NULL;
  }
  r->data = (double*)malloc((sizeof(double) * width * height));
  if (r->data == NULL) {
    fprintf(stderr, "error: failed to allocate memory for image data\n");
    return NULL;
  }

  r->width = width;
  r->height = height;
  r->sigma = sigma;
  p = r->data;
  for (j = 0; j < height; j ++) {
    for (i = 0; i < width; i ++, p ++) {
      c = fgetc(fp);
      if (c == EOF) {
	fprintf(stderr, "error: premature end of file\n");
	return NULL;
      }

      *p = (double)c/255.0;
    }
  }

  fclose(fp);
  return r;
}

double image_likelihood(void *user_arg, 
			int n,
			const double *values)
{
  struct image *img = (struct image*)user_arg;

  return image_likelihood_hierarchical(user_arg,
				       n,
				       values,
				       0,
				       1,
				       &(img->sigma),
				       NULL);
}

double image_likelihood_hierarchical(void *user_arg, 
				     int n,
				     const double *values,
				     int hierarchical,
				     int nhierarchicalparameters,
				     const double *hierarchicalvalues,
				     double *logdetce)
{
  struct image *img = (struct image*)user_arg;

  double A = values[P_A];
  double sigmax = values[P_sigmax];
  double sigmay = values[P_sigmay];
  double theta = values[P_theta];
  double x0 = values[P_x0];
  double y0 = values[P_y0];
  double B = values[P_B];
  
  double c2t = sqr(cos(theta));
  double s2t = sqr(sin(theta));
  double sx2 = sqr(sigmax);
  double sy2 = sqr(sigmay);

  double ga = c2t/(2.0 * sx2) + s2t/(2.0 * sy2);
  double gb = sin(2.0*theta) * (1.0/(4.0*sy2) - 1.0/(4.0*sx2));
  double gc = s2t/(2.0 * sx2) + c2t/(2.0 * sy2);

  double sum = 0.0;

  double sigma = hierarchicalvalues[0];

  int i;
  int j;

  double x;
  double y;

  double v;
  double dv;

  for (j = 0; j < img->height; j ++) {
    y = 1.0 - ((double)j + 0.5)/(double)img->height * 2.0;
    for (i = 0; i < img->width; i ++) {
      x = ((double)i + 0.5)/(double)img->width * 2.0 - 1.0;

      v = B + A * exp(-(ga * sqr(x - x0) + 
			2.0*gb*(x - x0)*(y - y0) + 
			gc*sqr(y - y0)));
      dv = v - img->data[j * img->width + i];

      sum += sqr(dv)/sqr(sigma);
	
    }
  }

  if (hierarchical) {

    *logdetce = 2.0 * (double)img->width * (double)img->height * log(sigma);

    /* printf("  %f\n", (*logdetce)); */
  }

  return sum/2.0;
}
