#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <map>
#include <cmath>
#include <string>
#include <fftw3.h>
#include <nvector/nvector_serial.h>

#include "fftmanip.h"

/* TYPE DEFINITIONS */

/* structure to hold complex numbers for linear algebra routines */
typedef struct {
 double re;
 double im;
} complex16;

/* FUNCTIONS */

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputFFTWVector(const char * fileName, fftw_complex * vec, double * x, int len);

/* Wrapper to outputFFTWVector, which fftshifts the output.
 */
void outputFFTWVectorShift(const char * fileName, fftw_complex * vec, double * x, int len);

/* prints out initial wave function.  Inputs are the wave function array and
 * the number of equations.
 */
void outputYData(realtype * ydata, int n, std::map<std::string, bool> &outs);

/* prints a vector W of length N */
void outputVector(realtype * W, int N, char * fileName);

/* prints the elementwise square of a complex vector */
void outputPsiSquare(complex16 * v, realtype * evals,  int N, char * fileName);

/* prints a complex vector v of length N */
void outputCVector(complex16 * v, int N, char * fileName);

/* prints a complex vector v of length N with M time steps */
void outputCVectorTime(complex16 * v, int N, int M, char * fileName);

/* prints a square matrix M of dimension N */
void outputSquareMatrix(realtype * M, int N, char * fileName);

/* prints a square matrix stored as a 2D array */
void output2DSquareMatrix(realtype ** M, int N, char * fileName);

#endif
