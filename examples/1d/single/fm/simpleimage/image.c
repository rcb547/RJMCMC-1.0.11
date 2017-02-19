
#include <stdio.h>
#include <math.h>

#include <rjmcmc/rjmcmc_random.h>

#include "image.h"

static unsigned char ftoc(double f)
{
  if (f < 0.0) {
    return 0;
  }
  if (f > 1.0) {
    return 255;
  }
  
  return (unsigned char)(255.0 * f);
}

static double sqr(double x)
{
  return x*x;
}


int image_createandwritepgm(const char *filename,
			    int width,
			    int height,
			    double A,
			    double sigmax,
			    double sigmay,
			    double theta,
			    double x,
			    double y,
			    double B,
			    double n)
{
  FILE *fp;

  double ga;
  double gb;
  double gc;
  double c2t;
  double s2t;
  double sx2;
  double sy2;

  int i;
  int j;
  
  double px;
  double py;

  c2t = sqr(cos(theta));
  s2t = sqr(sin(theta));
  sx2 = sqr(sigmax);
  sy2 = sqr(sigmay);

  ga = c2t/(2.0 * sx2) + s2t/(2.0 * sy2);
  gb = sin(2.0*theta) * (1.0/(4.0*sy2) - 1.0/(4.0*sx2));
  gc = s2t/(2.0 * sx2) + c2t/(2.0 * sy2);

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create image file\n");
    return -1;
  }

  fprintf(fp, "P5\n%d %d\n255\n", width, height);

  for (j = 0; j < height; j ++) {
    py = 1.0 - ((double)j + 0.5)/(double)height * 2.0;

    for (i = 0; i < width; i ++) {

      px = ((double)i + 0.5)/(double)width * 2.0 - 1.0;

      fputc(ftoc(A * exp(-(ga * sqr(px - x) + 
			   2.0*gb*(px - x)*(py - y) + 
			   gc*sqr(py - y))) + 
		 B +  
		 (rjmcmc_normal() * n)),
	    fp);
    }
  }

  fclose(fp);

  return 0;
}
