#include "numerical.h"

//#define DEBUG_NUMERICAL

/* returns the number of numbers in a file.  This way, it doesn't matter if
 * they are one per line or multiple per line.
 */
int Number_of_values (const char * nameOfFile) {
 FILE * input;
 double value;
 int numberOfValues = 0;

 input = fopen(nameOfFile, "r");

 if (input != NULL) {
  while (fscanf(input, "%lf", &value) != EOF) { numberOfValues++; }
  if (numberOfValues == 0 ) {
   fprintf(stderr, "WARNING: file %s is empty.\n", nameOfFile);
  }
 }
 else {
  fprintf(stderr, "WARNING [Number_of_values]: file %s does not exist.\n", nameOfFile);
  return -1;
 }

 return numberOfValues;
}

/* reads in the values from file; returns an array the length of the number of 
 * numbers in the file
 */
void Read_array_from_file (realtype * array, const char * nameOfFile, int numberOfValues) {
 FILE * input;
 int i = 0;

 input = fopen(nameOfFile,"r");

 if (input != NULL) {
  while (fscanf(input, "%lf", &array[i]) != EOF && i < numberOfValues) {
   i++;
  }
 }
 else {
  fprintf(stderr, "ERROR [Read_array_from_file]: file %s does not exist.\n", nameOfFile);
 }

 fclose(input);
}

/* Returns an array of length n with all values set to initializeValue. */
void Initialize_array(realtype * array, int n, realtype initializeValue) {
#ifdef DEBUG
 cout << "initializeValue is " << initializeValue << endl;
#endif

 int i;

 for (i = 0; i < n; i++) {
  array[i] = initializeValue;
 }
}

/* builds energies for a quasicontinuum (evenly spaced) */
void Build_continuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop) {
 
 int i;

 Energies[0] = BandEdge;	// the bottom of the conduction band is set
 
 // loop over the remaining states.  This way the top of the band will be at BandTop
 for (i = 1; i < numberOfStates; i++) {
  Energies[i] = Energies[i-1] + (BandTop-BandEdge)/(numberOfStates-1);
 }
}

/* populates a set of states according to a Fermi-Dirac distribution.
 * I'm not sure where the actual Fermi level is, so it defaults to 0.01 Eh
 * below the lowest-energy state in the set of states being populated.
 */
void Build_k_pops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp, int Nk) {

 int i;

 for (i = 0; i < Nk; i++) {
  kPops[i] = sqrt(1.0/(1.0 + exp((kEnergies[i]-kBandEdge+0.01)*3.185e5/(temp))));
#ifdef DEBUG
 cout << "\nk population at state " << i << " is: "
      << sqrt(1.0/(1.0 + exp((kEnergies[i]-kBandEdge+0.01)*3.185e5/(temp))));
#endif
 }
#ifdef DEBUG
 cout << endl;
#endif
}

/* populates a set of states according to a Gaussian distribution. */
void Build_k_pops_Gaussian(realtype * kPops, realtype * kEnergies, realtype kBandEdge, double sigma, double mu, int Nk) {
#ifdef DEBUG_NUMERICAL
  std::cout << "Conduction band edge is " << kBandEdge << std::endl;
  std::cout << "Gaussian mu is          " << mu << std::endl;
  std::cout << "Gaussian sigma is       " << sigma << std::endl;
  std::cout << "Number of k states is   " << Nk << std::endl;
#endif

 for (int ii = 0; ii < Nk; ii++) {
  kPops[ii] = sqrt((1/(sigma*sqrt(2*3.1415926535)))*exp(-pow((kEnergies[ii]-(kBandEdge+mu)),2)/(2*pow(sigma,2))));
#ifdef DEBUG_NUMERICAL
  std::cout << "\nk population at state " << ii << " is: " << kPops[ii];
#endif
 }
#ifdef DEBUG_NUMERICAL
 std::cout << std::endl;
#endif
}

/* creates a matrix of Franck-Condon factors between two displaced harmonic
 * oscillators.
 */
void Build_Franck_Condon_factors (realtype ** FCmat, double g, int numM, int numN) {
 int m, n;

 FCmat[0][0] = exp(-pow(g,2)/2);	// first element
 for (m = 1; m < numM; m++)		// first column
  FCmat[m][0] = g/sqrt((double)m)*FCmat[m-1][0];
 for (n = 1; n < numN; n++)		// first row
  FCmat[0][n] = -1*g/sqrt((double)n)*FCmat[0][n-1];
 for (m = 1; m < numM; m++)		// recursion for the rest of the matrix
  for (n = 1; n < numN; n++)
   FCmat[m][n] = g/sqrt((double)m)*FCmat[m-1][n] + sqrt((double)n/(double)m)*FCmat[m-1][n-1];
}

/* returns the coupling as a function of energy E given that the middle of the
 * band is at position mid.
 * Eq. 31 in Ramakrishna et al., JCP 2001, 115, 2743-2756
 */
double parabolicV(double Vee, double E, double bandEdge, double bandTop) {
// set band edge as zero energy
E = E - bandEdge;
double mid = (bandTop - bandEdge)/2.0;

#ifdef DEBUG
 fprintf(stdout, "coupling at (E - band edge) = %.9e: %.9e\n", E, Vee*sqrt(sqrt(pow(mid,2) - pow((E - mid),2))/mid));
#endif
 // test whether at the very top or bottom of band
 if (abs(abs(E - mid) - mid) < 1e-10) {
#ifdef DEBUG
  fprintf(stdout, "(at bottom or top of band edge, returning 0.0)\n");
#endif
  return 0.0;
 }
 return Vee*sqrt(sqrt(pow(mid,2) - pow((E - mid),2))/mid);
}

/* gives the value of a laser pulse (electric field) at time t */
realtype pump(realtype t, double pumpFWHM, double pumpAmpl, double pumpPeak, double pumpFreq, double pumpPhase) {
 double sigma = pumpFWHM/2.35482005;
 return pumpAmpl*exp((-pow(t-pumpPeak, 2))/(2*pow(sigma, 2)))*cos(pumpFreq*t + pumpPhase);
}

/* normalizes an N_Vector so that the populations of all states are
 * normalized to value 'total'
 */
int Normalize_NV(N_Vector nv, realtype total) {
 int i;
 realtype summ = 0;

 for (i = 0; i < NV_LENGTH_S(nv); i++) {
  summ += (NV_Ith_S(nv, i)*NV_Ith_S(nv, i));
 }
 summ = sqrt(summ);
 for (i = 0; i < NV_LENGTH_S(nv); i++) {
  NV_Ith_S(nv, i) = total*NV_Ith_S(nv,i)/summ;
 }

 return 0;
}

/* compute the six-point derivative of an array.
 * Assumes array elements are evenly spaced.
 * Output array is six elements shorter than input.
 */
int Derivative(double *inputArray, int inputLength, double *outputArray, double timestep) {
 int i;		// counter

 if (inputLength < 6 ) {
  fprintf(stderr, "ERROR [Derivative]: array has length less than 6 elements, cannot proceed");
  return -1;
 }

 for (i = 2; i < inputLength-3; i++) {
  outputArray[i-2] = (2* inputArray[i+3]
		     -15*inputArray[i+2]
		     +60*inputArray[i+1]
		     -20*inputArray[i]
		     -30*inputArray[i-1]
		     +3 *inputArray[i-2])/(60*timestep);
 }

 return 0;
}

/* Finds the first maximum in an array (the first point where the next
 * point is smaller in value).
 */
realtype Find_first_array_maximum (realtype * inputArray, int num) {
 int i;
 realtype currentMax = inputArray[0];

 for (i = 1; i < num; i++) {
  if (inputArray[i] > currentMax)
   currentMax = inputArray[i];
  if (inputArray[i] < currentMax)
   break;
 }

 return currentMax;
}

/* This function returns the index of the first maximum in an array.
 * Warning: will return 0 as index of first maximum if the second array element
 * is less than the first.  This may not be what you want.
 */
int Find_first_array_maximum_index (realtype * inputArray, int num) {
 int i;
 realtype currentMax = inputArray[0];
 int currentMax_index = 0;

 for (i = 1; i < num; i++) {
  if (inputArray[i] > currentMax) {
   currentMax = inputArray[i];
   currentMax_index = i;
  }
  if (inputArray[i] < currentMax)
   break;
 }

 return currentMax_index;
}

/* returns index of first maximum in an array. */
int Find_array_maximum_index (realtype * inputArray, int num) {
 int i;
 realtype currentMax = inputArray[0];
 int currentMax_index = 0;

 for (i = 1; i < num; i++) {
  if (inputArray[i] > currentMax) {
   currentMax = inputArray[i];
   currentMax_index = i;
  }
 }

 return currentMax_index;
}

