#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <map>
#include <cmath>
#include <string>
#include <fftw3.h>
#include <nvector/nvector_serial.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "conversions.hpp"
#include "numerical.hpp"
#include "fftmanip.hpp"
#include "params.hpp"
#include "rhs.hpp"

/* TYPE DEFINITIONS */

/* structure to hold complex numbers for linear algebra routines */
typedef struct {
  double re;
  double im;
} complex16;

/* FUNCTIONS */

bool isOutput(std::map<const std::string, bool> &myMap, const std::string myStr);

std::string outputFileName(char * fileName, PARAMETERS * p);

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
void outputWavefunction(realtype * psi, int n);

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

void outputIntegratedDM(const std::string fileName, const int start, const int end,
    const realtype * dmt, struct PARAMETERS * p);

void outputIntegratedWfn(const std::string fileName, const int start, const int end,
    const realtype * wfnt, struct PARAMETERS * p);

void outputIntegralDM(const std::string fileName, const int start, const int end,
    const realtype * dmt, struct PARAMETERS * p);

void outputIntegralWfn(const std::string fileName, const int start, const int end,
    const realtype * wfnt, struct PARAMETERS * p);

void outputXProbsWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p);

void outputXProbs(char * fileName, int start, int end, realtype * dmt,
    struct PARAMETERS * p);

void outputtXprobWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p);

void outputtXprob(char * fileName, int start, int end, realtype * dmt,
    struct PARAMETERS * p);

void outputDMZ(realtype * dmt, struct PARAMETERS * p);

void outputDMCoherences(realtype * dmt, struct PARAMETERS * p);

void outputDMIm(char * fileName, realtype * dmt, struct PARAMETERS * p);

void outputDMRe(char * fileName, realtype * dmt, struct PARAMETERS * p);

void outputEnergy(char * fileName, struct PARAMETERS * p);

void outputEnergyExpWfn(const char * fileName, struct PARAMETERS * p);

void outputTimes(char * fileName, struct PARAMETERS * p);

void outputEnergyExp(char * fileName, realtype * dmt,
    struct PARAMETERS * p);

void outputDynamicMu(char * fileName, realtype * dmt, int bandFlag, struct PARAMETERS * p);

void outputMuFromPops(char * fileName, realtype * dmt,
    struct PARAMETERS * p);

void outputTorsion(struct PARAMETERS * p, char * fileName);

void outputRTA(std::map<const std::string, bool> &outs,
    realtype * dmt, struct PARAMETERS * p);

void outputCouplings(struct PARAMETERS * p, char * fileName);

void findPeaksWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p);

void outputDeriv(std::string outFile, int n, realtype * deriv, struct PARAMETERS * p);

void outputDerivsWfn(std::map<const std::string, bool> &outs, realtype * wfnt,
    struct PARAMETERS * p);

void outputDerivsDM(std::map<const std::string, bool> &outs, realtype * dmt,
    struct PARAMETERS * p);

void computeGeneralOutputs(std::map<const std::string, bool> &outs,
    struct PARAMETERS * p);

void computeWfnOutput(realtype * wfnt, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p);

void computeDMOutput(realtype * dmt, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p);
#endif
