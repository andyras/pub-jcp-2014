#ifndef __NUMERICAL_H__
#define __NUMERICAL_H__

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <nvector/nvector_serial.h>
#include <iostream>

/* returns the number of numbers in a file.  This way, it doesn't matter if
 * they are one per line or multiple per line.
 */
int Number_of_values (const char * nameOfFile);

/* reads in the values from file; returns an array the length of the number of 
 * numbers in the file
 */
void Read_array_from_file (realtype * array, const char * nameOfFile, int numberOfValues);

/* Returns an array of length n with all values set to initializeValue. */
void Initialize_array(realtype * array, int n, realtype initializeValue);

/* builds energies for a quasicontinuum (evenly spaced) */
void Build_continuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop);

/* populates a set of states according to a Fermi-Dirac distribution.
 * I'm not sure where the actual Fermi level is, so it defaults to 0.01 Eh
 * below the lowest-energy state in the set of states being populated.
 */
void Build_k_pops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp, int Nk);

/* populates a set of states according to a Gaussian distribution. */
void Build_k_pops_Gaussian(realtype * kPops, realtype * kEnergies, realtype kBandEdge, double sigma, double mu, int Nk);

/* creates a matrix of Franck-Condon factors between two displaced harmonic
 * oscillators.
 */
void Build_Franck_Condon_factors (realtype ** FCmat, double g, int numM, int numN);

/* returns the coupling as a function of energy E given that the middle of the
 * band is at position mid.
 * Eq. 31 in Ramakrishna et al., JCP 2001, 115, 2743-2756
 */
double parabolicV(double Vee, double E, double bandEdge, double bandTop);

/* gives the value of a laser pulse (electric field) at time t */
realtype pump(realtype t, double pumpFWHM, double pumpAmpl, double pumpPeak, double pumpFreq, double pumpPhase);

/* normalizes an N_Vector so that the populations of all states are
 * normalized to value 'total'
 */
int Normalize_NV(N_Vector nv, realtype total);

/* compute the six-point derivative of an array.
 * Assumes array elements are evenly spaced.
 * Output array is six elements shorter than input.
 */
int Derivative(double *inputArray, int inputLength, double *outputArray, double timestep);

/* Finds the first maximum in an array (the first point where the next
 * point is smaller in value).
 */
realtype Find_first_array_maximum (realtype * inputArray, int num);

/* This function returns the index of the first maximum in an array.
 * Warning: will return 0 as index of first maximum if the second array element
 * is less than the first.  This may not be what you want.
 */
int Find_first_array_maximum_index (realtype * inputArray, int num);

/* returns index of first maximum in an array. */
int Find_array_maximum_index (realtype * inputArray, int num);

#endif
