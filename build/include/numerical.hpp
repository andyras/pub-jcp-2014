#ifndef __NUMERICAL_H__
#define __NUMERICAL_H__

#include <map>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <nvector/nvector_serial.h>
#include <mkl.h>

#include "params.hpp"
#include "constants.hpp"

void updateHamiltonian(PARAMETERS * p, realtype t);

/* returns the number of numbers in a file.  This way, it doesn't matter if
 * they are one per line or multiple per line.
 */
int numberOfValuesInFile(const char * nameOfFile);

/* reads in the values from file; returns an array the length of the number of 
 * numbers in the file
 */
void readArrayFromFile(realtype * array, const char * nameOfFile, int numberOfValues);

void readVectorFromFile(std::vector<realtype> & v, const char * fileName, int n);

/* Returns an array of length n with all values set to initializeValue. */
void initializeArray(realtype * array, int n, realtype initializeValue);

/* builds energies for a quasicontinuum (evenly spaced) */
void buildContinuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop);

void buildParabolicBand(realtype * energies, int n, double bandEdge, int flag, PARAMETERS * p);

/* populates a set of states according to a Fermi-Dirac distribution.
 * I'm not sure where the actual Fermi level is, so it defaults to 0.01 Eh
 * below the lowest-energy state in the set of states being populated.
 */
void buildKPops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp, int Nk);

/* populates a set of states according to a Gaussian distribution. */
void buildKPopsGaussian(realtype * kPops, realtype * kEnergies, realtype kBandEdge, double sigma, double mu, int Nk);

/* returns the coupling as a function of energy E given that the middle of the
 * band is at position mid.
 * Eq. 31 in Ramakrishna et al., JCP 2001, 115, 2743-2756
 */
double parabolicV(double Vee, double E, double bandEdge, double bandTop);

realtype gaussPulse(realtype t, double pumpFWHM, double pumpAmpl, double pumpPeak, double pumpFreq, double pumpPhase);

/* normalizes an N_Vector so that the populations of all states are
 * normalized to value 'total'
 */
int Normalize_NV(N_Vector nv, realtype total);

/* compute the six-point derivative of an array.
 * Assumes array elements are evenly spaced.
 * Output array is six elements shorter than input.
 */
int Derivative(double *inputArray, int inputLength, double *outputArray, double timestep);

void arrayDeriv(double * in, int n, int m, int dim, double * out, double dt);

/* Riemann sum of an array (values) at time points (time).
 * Does not assume equal spacing in time.
 */
realtype integrateArray(realtype * values, realtype * time, int num);

/* Returns maximum element in an array. */
realtype findArrayMaximum(realtype * inputArray, int num);

/* Finds the first maximum in an array (the first point where the next
 * point is smaller in value).
 */
realtype findFirstArrayMaximum(realtype * inputArray, int num);

/* This function returns the index of the first maximum in an array.
 * Warning: will return 0 as index of first maximum if the second array element
 * is less than the first.  This may not be what you want.
 */
int findFirstArrayMaximumIndex(realtype * inputArray, int num);

/* returns index of first maximum in an array. */
int findArrayMaximumIndex(realtype * inputArray, int num);

void buildCoupling (realtype ** vArray, struct PARAMETERS * p,
    std::map<const std::string, bool> &outs);

/* builds a Hamiltonian from site energies and couplings. */
void buildHamiltonian(realtype * H, std::vector<realtype> & energy, realtype ** V, struct PARAMETERS * p);

void updateDM(N_Vector dm, realtype * dmt, int timeStep, struct PARAMETERS * p);

void updateWfn(N_Vector wfn, realtype * wfnt, int timeStep, struct PARAMETERS * p);

/* returns sign of number (+1/-1/0) */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif
