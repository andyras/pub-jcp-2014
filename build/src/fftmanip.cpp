#include <math.h>

#include "output.hpp"
#include "fftmanip.hpp"

/* This function implements something similar to MATLAB's 'fftshift':
 * it takes the FFT'd vector invec and puts a shifted copy in outvec.*/
void fftshift(fftw_complex * invec, fftw_complex * outvec, int n) {
  // number of elements to shift output array
  int shift = n/2;
  // shifted counter
  int ii;

  // do the shift!
  for (int i = 0; i < n; i++) {
    ii = (i+shift) % n;
    outvec[ii][0] = invec[i][0];
    outvec[ii][1] = invec[i][1];
  }

  return;
}

void fftshift_double(double * invec, double * outvec, int n) {
  // number of elements to shift output array
  int shift = n/2;
  // shifted counter
  int ii;

  // do the shift!
  for (int i = 0; i < n; i++) {
    ii = (i+shift) % n;
    outvec[ii] = invec[i];
    outvec[ii] = invec[i];
  }

  return;
}

/* This function writes the FT of the input time-sampled vector to a file.
 * It assumes evenly spaced time data (with at least two points).
 */
void writeFT(const char * fileName, double * invec, double * times, int n) {
  // compute time spacing
  double dt = times[1] - times[0];
  // compute energy spacing
  double dE = 2.0*M_PI/(dt*n);
  // build energies array
  double * energies = new double [n];
  for (int i = 0; i < n; i++) {
    energies[i] = (i-n/2)*dE;
  }

  // FFTW arrays
  fftw_complex * in;
  fftw_complex * out;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

  // make FFT plan
  fftw_plan fftw_fwd;
  fftw_fwd = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  // set up in array
  for (int ii = 0; ii < n; ii++) {
    in[ii][0] = invec[ii];
    in[ii][1] = 0;
  }

  // do the FFT
  fftw_execute(fftw_fwd);

  // print the result
  outputFFTWVectorShift(fileName, out, energies, n);

  // clean up
  delete [] energies;
  return;
}
