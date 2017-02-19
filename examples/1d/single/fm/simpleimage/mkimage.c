
/* mkimage
 
   A simple program to generate a synthetic image to use in the simpleimage
   example program.

   The image is a grey scale image which has a 2D gaussian shape within it
   that the simple image example infers the parameters of.

   The image is constructed using the formula:

   I(x,y) = B + n + A*exp(ga*(x - x0)^2 + 
                          2.0*gb*(px - x)*(py - y) + 
                          gc*(y - y0)^2)

   Where B = The background intensity
         n = Noise (gaussian)
         A = Shape peak intensity
         x0 = Horizontal centre (centre of image is 0,0)
         y0 = Vertical centre
         ga, gb, gc = Coefficients of gaussian derived from sigma_x, sigma_y
         and theta (the gaussian shape parameters).

   The image is saved as a PGM file to test.pgm.
*/


#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include <rjmcmc/rjmcmc_random.h>

int main(int argc, char *argv[])
{
  FILE *fp;
  int i;
  int j;
  double px;
  double py;

  int width;
  int height;
  
  double A;
  double sigmax;
  double sigmay;
  double theta;
  double x;
  double y;
  double B;
  double n;

  width = 256;
  height = 256;

  A = 0.5;
  sigmax = 0.1;
  sigmay = 0.2;
  theta = M_PI/4.0;
  x = -0.1;
  y = 0.2;
  B = 0.3;
  n = 0.15;

  if (image_createandwritepgm("image.pgm",
			      width, height,
			      A,
			      sigmax,
			      sigmay,
			      theta,
			      x,
			      y,
			      B,
			      n) < 0) {
    fprintf(stderr, "error: failed to write image.pgm\n");
    return -1;
  }
  
  /*
   * Write the actual image without noise for comparison.
   */
  if (image_createandwritepgm("image_real.pgm",
			      width, height,
			      A,
			      sigmax,
			      sigmay,
			      theta,
			      x,
			      y,
			      B,
			      0.0) < 0) {
    fprintf(stderr, "error: failed to write image.pgm\n");
    return -1;
  }

  return 0;
}
