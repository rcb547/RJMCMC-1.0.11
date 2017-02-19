#ifndef IMAGE_H
#define IMAGE_H

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
			    double n);

#endif /* IMAGE_H */
