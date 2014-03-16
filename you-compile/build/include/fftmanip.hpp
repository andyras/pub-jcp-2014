#ifndef __FFTMANIP_H__
#define __FFTMANIP_H__

#include <fftw3.h>

/* This function implements something similar to MATLAB's 'fftshift':
 * it takes the FFT'd vector invec and puts a shifted copy in outvec.*/
void fftshift(fftw_complex * invec, fftw_complex * outvec, int n);

/* Same as fftshift, but for doubles! */
void fftshift_double(double * invec, double * outvec, int n);

/* This function writes the FT of the input time-sampled vector to a file.
 * It assumes evenly spaced time data (with at least two points).
 */
void writeFT(const char * fileName, double * invec, double * times, int n);

#endif
