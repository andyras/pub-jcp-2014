#include "numerical.hpp"
//#define DEBUG

//#define DEBUG_NUMERICAL
//#define DEBUG_BUILDCOUPLING
//#define DEBUG_UPDATEDM
//#define DEBUG_UPDATEWFN
//#define DEBUG_BUILDHAMILTONIAN
//#define DEBUG_READVECTOR

/* Updates the Hamiltonian with the time-dependent torsional coupling
 * and laser field.
 */
void updateHamiltonian(PARAMETERS * p, realtype t) {
  // TODO unpack NEQ from p

  // get pointer to H
  realtype * H = &(p->H)[0];

  //// first handle torsion
  if (p->torsion) {
    double torsionValue = p->torsionV->value(t);
#ifdef DEBUG_TORSION
    std::cout << "Value of torsion-mediated coupling is " << torsionValue << std::endl;
#endif

    // bridge is off, coupling is between k and c states
    if (!(p->bridge_on)) {
#ifdef DEBUG_RHS
      std::cout << "torsion between k and c states" << std::endl;
#endif
      for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
	for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
	  H[ii*p->NEQ + jj] = torsionValue;
	  H[jj*p->NEQ + ii] = torsionValue;
	}
      }
    }
    // torsion is at first bridge coupling
    else if (p->torsionSite == 0) {
#ifdef DEBUG_RHS
      std::cout << "torsion between k states and bridge" << std::endl;
#endif
      for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
	H[ii*p->NEQ + p->Ib] = torsionValue;
	H[p->Ib*p->NEQ + ii] = torsionValue;
      }
    }
    // torsion is at last bridge coupling
    else if (p->torsionSite == p->Nb) {
#ifdef DEBUG_RHS
      std::cout << "torsion between bridge and c states" << std::endl;
#endif
      for (int ii = p->Ic; ii < (p->Ic + p->Nc); ii++) {
	H[ii*p->NEQ + p->Ib + p->Nb - 1] = torsionValue;
	H[(p->Ib + p->Nb - 1)*p->NEQ + ii] = torsionValue;
      }
    }
    // torsion is between bridge sites
    else {
#ifdef DEBUG_RHS
      std::cout << "torsion between bridge sites " << p->torsionSite - 1
	<< " and " << p->torsionSite << "." << std::endl;
#endif
      H[(p->Ib + p->torsionSite - 1)*p->NEQ + p->Ib + p->torsionSite] = torsionValue;
      H[(p->Ib + p->torsionSite)*p->NEQ + p->Ib + p->torsionSite - 1] = torsionValue;
    }
  }

  //// now handle pump pulse
  double laserCoupling = 0.0;
  if (p->laser_on) {
    laserCoupling = gaussPulse(t, p->pumpFWHM, p->pumpAmpl, p->pumpPeak, p->pumpFreq, p->pumpPhase);
#ifdef DEBUG_RHS
    std::cout << "Value of laser coupling between valence and conduction bands is " << laserCoupling << std::endl;
#endif
    // coupling is between valence and conduction bands
    for (int ii = p->Il; ii < (p->Il + p->Nl); ii++) {
      for (int jj = p->Ik; jj < (p->Ik + p->Nk); jj++) {
	H[(ii)*p->NEQ + jj] = laserCoupling;
	H[(jj)*p->NEQ + ii] = laserCoupling;
      }
    }
  }

  // make sparse version of Hamiltonian
  int job [6] = {0, 0, 0, 2, p->NEQ2, 1};
  int info = 0;

  mkl_ddnscsr(&job[0], &(p->NEQ), &(p->NEQ), &(p->H)[0], &(p->NEQ), &(p->H_sp)[0],
      &(p->H_cols)[0], &(p->H_rowind)[0], &info);

  return;
}

/* returns the number of numbers in a file.  This way, it doesn't matter if
 * they are one per line or multiple per line.
 */
int numberOfValuesInFile(const char * nameOfFile) {
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
    fprintf(stderr, "WARNING [numberOfValuesInFile]: file %s does not exist.\n", nameOfFile);
    return -1;
  }

  fclose(input);

  return numberOfValues;
}

/* reads in the values from file; returns an array the length of the number of 
 * numbers in the file
 */
void readArrayFromFile(realtype * array, const char * nameOfFile, int numberOfValues) {
  FILE * input;
  int i = 0;

  input = fopen(nameOfFile,"r");

  if (input != NULL) {
    while (fscanf(input, "%lf", &array[i]) != EOF && i < numberOfValues) {
      i++;
    }
  }
  else {
    fprintf(stderr, "ERROR [readArrayFromFile]: file %s does not exist.\n", nameOfFile);
  }

  fclose(input);

  return;
}

/* reads in values from a file to a vector. */
void readVectorFromFile(std::vector<realtype> & v, const char * fileName, int n) {
  // resize vector according to number of lines in file
  v.resize(numberOfValuesInFile(fileName));
#ifdef DEBUG_READVECTOR
  std::cout << numberOfValuesInFile(fileName) << " values in " << fileName << std::endl;
#endif
  if (numberOfValuesInFile(fileName) != n) {
    std::cerr << "ERROR reading in vector from " << fileName << ": wrong number of values in file." << std::endl;
    _exit(1);
  }

  // read in the file
  std::ifstream in(fileName);

  for (int ii = 0; ii < n; ii++) {
    in >> v[ii];
  }

  in.close();

  return;
}

/* Returns an array of length n with all values set to initializeValue. */
void initializeArray(realtype * array, int n, realtype initializeValue) {
#ifdef DEBUG
  std::cout << "initializeValue is " << initializeValue << std::endl;
#endif

  for (int ii = 0; ii < n; ii++) {
    array[ii] = initializeValue;
  }

  return;
}

/* builds energies for a quasicontinuum (evenly spaced) */
void buildContinuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop) {
  if (numberOfStates == 0) {
    return;
  }

  Energies[0] = BandEdge;	// the bottom of the conduction band is set

  // loop over the remaining states.  This way the top of the band will be at BandTop
  for (int ii = 1; ii < numberOfStates; ii++) {
    Energies[ii] = Energies[ii-1] + (BandTop-BandEdge)/(numberOfStates-1);
  }

  return;
}

void buildParabolicBand(realtype * energies, int n, double bandEdge, int flag, PARAMETERS * p) {
  int s;	// sign +/-
  double m;	// mass of electron/hole

  // determine conduction vs. valence band, electron/hole masses
  if (flag == CONDUCTION) {
    s = 1;
    m = p->me;
  }
  else if (flag == VALENCE) {
    s = -1;
    m = p->mh;
  }

  // assign energies
  for (int ii = 0; ii < n; ii++) {
    energies[ii] = bandEdge + s*ii*ii/(2*m*pow(p->X2,2));
  }

  return;
}

/* populates a set of states according to a Fermi-Dirac distribution.
 * I'm not sure where the actual Fermi level is, so it defaults to 0.01 Eh
 * below the lowest-energy state in the set of states being populated.
 */
void buildKPops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp, int Nk) {

  for (int ii = 0; ii < Nk; ii++) {
    kPops[ii] = sqrt(1.0/(1.0 + exp((kEnergies[ii]-kBandEdge+0.01)*3.185e5/(temp))));
#ifdef DEBUG
    std::cout << "\nk population at state " << ii << " is: "
      << sqrt(1.0/(1.0 + exp((kEnergies[ii]-kBandEdge+0.01)*3.185e5/(temp))));
#endif
  }
#ifdef DEBUG
  std::cout << std::endl;
#endif
}

/* populates a set of states according to a Gaussian distribution. */
void buildKPopsGaussian(realtype * kPops, realtype * kEnergies, realtype kBandEdge, double sigma, double mu, int Nk) {
#ifdef DEBUG_NUMERICAL
  std::cout << "Conduction band edge is " << kBandEdge << std::endl;
  std::cout << "Gaussian mu is          " << mu << std::endl;
  std::cout << "Gaussian sigma is       " << sigma << std::endl;
  std::cout << "Number of k states is   " << Nk << std::endl;
#endif

  for (int ii = 0; ii < Nk; ii++) {
    // take the square root so that populations have proper width given by sigma
    kPops[ii] = sqrt((1/(sigma*sqrt(2*3.1415926535)))*exp(-pow((kEnergies[ii]-(kBandEdge+mu)),2)/(2*pow(sigma,2))));
#ifdef DEBUG_NUMERICAL
    std::cout << "\nk population at state " << ii << " is: " << kPops[ii];
#endif
  }
#ifdef DEBUG_NUMERICAL
  std::cout << std::endl;
#endif
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

/* gives the value of a Gaussian laser pulse (electric field) at time t */
realtype gaussPulse(realtype t, double pumpFWHM, double pumpAmpl, double pumpPeak, double pumpFreq, double pumpPhase) {
  double sigma = pumpFWHM/2.35482005;	// conversion between FWHM and Gaussian 2nd moment
  return pumpAmpl*exp((-pow(t-pumpPeak, 2))/(2*pow(sigma, 2)))*cos(pumpFreq*t + pumpPhase);
}

/* normalizes an N_Vector so that the populations of all states are
 * normalized to value 'total'
 */
int Normalize_NV(N_Vector nv, realtype total) {
  realtype summ = 0;

  for (int ii = 0; ii < NV_LENGTH_S(nv); ii++) {
    summ += (NV_Ith_S(nv, ii)*NV_Ith_S(nv, ii));
  }
  summ = sqrt(summ);
  for (int ii = 0; ii < NV_LENGTH_S(nv); ii++) {
    NV_Ith_S(nv, ii) = total*NV_Ith_S(nv,ii)/summ;
  }

  return 0;
}

/* compute the six-point derivative of an array.
 * Assumes array elements are evenly spaced.
 * Output array is six elements shorter than input.
 */
int Derivative(double *inputArray, int inputLength, double *outputArray, double timestep) {
  if (inputLength < 6 ) {
    fprintf(stderr, "ERROR [Derivative]: array has length less than 6 elements, cannot proceed");
    return -1;
  }

  for (int ii = 2; ii < inputLength-3; ii++) {
    outputArray[ii-2] = (2* inputArray[ii+3]
	-15*inputArray[ii+2]
	+60*inputArray[ii+1]
	-20*inputArray[ii]
	-30*inputArray[ii-1]
	+3 *inputArray[ii-2])/(60*timestep);
  }

  return 0;
}

/* Compute the six-point time derivative of a 2D array.
 * Assumes array elements are evenly spaced in time.
 * Output array is six elements shorter than the input, skipping first two and
 * last three points.  The output array must be preallocated.
 */
void arrayDeriv(double * in, int nt, int m, int dim, double * out, double dt) {
  // n is length of array (number of time points)
  // m is width of array (number of populations)
  // dim is p->NEQ
  if (nt < 6) {
    std::cerr << "Error [" << __FUNCTION__ << "]: array must be at least six elements." << std::endl;
    _exit(-1);
  }

  for (int ii = 0; ii < m; ii++) {
    for (int jj = 2; jj < (nt-3); jj++) {
      out[ii*nt + (jj-2)] = (2*in[ii*nt + jj+3]
	  - 15*in[ii*nt + jj+2]
	  + 60*in[ii*nt + jj+1]
	  - 20*in[ii*nt + jj]
	  - 30*in[ii*nt + jj-1]
	  +  3*in[ii*nt + jj-2])
	/
	(60*dt);
    }
  }

  return;
}

/* Riemann sum of an array (values) at time points (time).
 * Does not assume equal spacing in time.
 */
realtype integrateArray(realtype * values, realtype * time, int num) {
  realtype riemann = 0;

  for (int ii = 0; ii < num-1; ii++)
    riemann += (values[ii+1] + values[ii])*(time[ii+1]-time[ii])/2;

  return riemann;
}

/* Returns maximum element in an array. */
realtype findArrayMaximum(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];

  for (int ii = 1; ii < num; ii++)
    if (inputArray[ii] > currentMax)
      currentMax = inputArray[ii];

  return currentMax;
}

/* Finds the first maximum in an array (the first point where the next
 * point is smaller in value).
 */
realtype findFirstArrayMaximum(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];

  for (int ii = 1; ii < num; ii++) {
    if (inputArray[ii] > currentMax)
      currentMax = inputArray[ii];
    if (inputArray[ii] < currentMax)
      break;
  }

  return currentMax;
}

/* This function returns the index of the first maximum in an array.
 * Warning: will return 0 as index of first maximum if the second array element
 * is less than the first.  This may not be what you want.
 */
int findFirstArrayMaximumIndex(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];
  int currentMax_index = 0;

  for (int ii = 1; ii < num; ii++) {
    if (inputArray[ii] > currentMax) {
      currentMax = inputArray[ii];
      currentMax_index = ii;
    }
    if (inputArray[ii] < currentMax)
      break;
  }

  return currentMax_index;
}

/* returns index of first maximum in an array. */
int findArrayMaximumIndex(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];
  int currentMax_index = 0;

  for (int ii = 1; ii < num; ii++) {
    if (inputArray[ii] > currentMax) {
      currentMax = inputArray[ii];
      currentMax_index = ii;
    }
  }

  return currentMax_index;
}

/* assign coupling constants to global array V */
void buildCoupling (realtype ** vArray, struct PARAMETERS * p,
    std::map<const std::string, bool> &outs) {

  double Vkc;	// coupling between bulk and QD
  double Vkb1;	// coupling between bulk and first bridge
  double VbNc;	// coupling between last bridge and QD

  // initialize the coupling array
  for (int ii = 0; ii < p->NEQ; ii++) {
    for (int jj = 0; jj < p->NEQ; jj++) {
      vArray[ii][jj] = 0.0;
    }
  }

#ifdef DEBUG_BUILDCOUPLING
  if (p->bridge_on) {
    for (int ii = 0; ii < p->Nb + 1; ii++) {
      std::cout << "p->Vbridge[" << ii << "] is ";
      std::cout << p->Vbridge[ii] << "\n";
    }
  }
#endif

  // bridge
  if (p->bridge_on) {
    // coupling between k and b1
    if ((p->scale_bubr) && (p->Nk > 1)) {
      Vkb1 = sqrt(p->Vbridge[0]*(p->kBandTop-p->kBandEdge)/(p->Nk-1));
    }
    else {
      Vkb1 = p->Vbridge[0];
    }
    if (p->parabolicCoupling) {
      for (int ii = 0; ii < p->Nk; ii++) {
	vArray[p->Ik+ii][p->Ib] = parabolicV(Vkb1, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
	vArray[p->Ib][p->Ik+ii] = parabolicV(Vkb1, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
      }
    }
    else {
      for (int ii = 0; ii < p->Nk; ii++) {
	vArray[p->Ik+ii][p->Ib] = Vkb1;
	vArray[p->Ib][p->Ik+ii] = Vkb1;
      }
    }

    // coupling between bN and c
    if ((p->scale_brqd) && (p->Nc > 1)) {
      VbNc = p->Vbridge[p->Nb]/sqrt(p->Nc-1);
    }
    else {
      VbNc = p->Vbridge[p->Nb];
    }
    for (int ii = 0; ii < p->Nc; ii++) {
      vArray[p->Ic+ii][p->Ib+p->Nb-1] = VbNc;
      vArray[p->Ib+p->Nb-1][p->Ic+ii] = VbNc;
    }

    // coupling between bridge states
    for (int ii = 0; ii < p->Nb - 1; ii++) {
      vArray[p->Ib+ii][p->Ib+ii+1] = p->Vbridge[ii+1];
      vArray[p->Ib+ii+1][p->Ib+ii] = p->Vbridge[ii+1];
    }
  }
  // no bridge
  else {				
    // scaling
    if ((p->scale_buqd) && (p->Nk > 1)) {
      Vkc = sqrt(p->Vnobridge[0]*(p->kBandTop-p->kBandEdge)/(p->Nk-1));
    }
    else {
      Vkc = p->Vnobridge[0];
    }

    // parabolic coupling of bulk band to QD
    if (p->parabolicCoupling) {
      for (int ii = 0; ii < p->Nk; ii++) {
	for (int jj = 0; jj < p->Nc; jj++) {
	  vArray[p->Ik+ii][p->Ic+jj] = parabolicV(Vkc, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
	  vArray[p->Ic+jj][p->Ik+ii] = parabolicV(Vkc, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
	}
      }
    }
    else {
      for (int ii = 0; ii < p->Nk; ii++) {
	for (int jj = 0; jj < p->Nc; jj++) {
	  vArray[p->Ik+ii][p->Ic+jj] = Vkc;
	  vArray[p->Ic+jj][p->Ik+ii] = Vkc;
	}
      }
    }
  }

#ifdef DEBUG
  std::cout << "\nCoupling matrix:\n";
  for (int ii = 0; ii < p->NEQ; ii++) {
    for (int jj = 0; jj < p->NEQ; jj++)
      std::cout << std::scientific << vArray[ii][jj] << " ";
    std::cout << std::endl;
  }
#endif
}

/* builds a Hamiltonian from site energies and couplings. */
void buildHamiltonian(realtype * H, std::vector<realtype> & energy, realtype ** V, struct PARAMETERS * p) {
  // indices
  int idx1, idx2;
  int N = p->NEQ;

#ifdef DEBUG_BUILDHAMILTONIAN
  fprintf(stderr, "Assigning diagonal elements of Hamiltonian.\n");
#endif
  for (int ii = 0; ii < N; ii++) {
    // diagonal
    H[ii*N + ii] = energy[ii];
#ifdef DEBUG_BUILDHAMILTONIAN
    std::cout << "diagonal element " << ii << " of H is " << energy[ii] << "\n";
#endif
  }

  if (p->bridge_on) {
    // assign bulk-bridge coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bulk-bridge coupling elements in Hamiltonian.\n");
#endif
    idx2 = p->Ib;
    for (int ii = 0; ii < p->Nk; ii++) {
      idx1 = p->Ik + ii;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
      H[idx1*N + idx2] = V[idx1][idx2];
      H[idx2*N + idx1] = V[idx2][idx1];
    }
    fprintf(stderr, "Done assigning bulk-bridge coupling elements in Hamiltonian.\n");
    // assign bridge-bridge couplings
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bridge-bridge coupling elements in Hamiltonian.\n");
#endif
    for (int ii = 0; ii < (p->Nb-1); ii++) {
      idx1 = p->Ib + ii;
      idx2 = p->Ib+ ii + 1;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
      H[idx1*N + idx2] = V[idx1][idx2];
      H[idx2*N + idx1] = V[idx2][idx1];
    }
    fprintf(stderr, "Done assigning bridge-bridge coupling elements in Hamiltonian.\n");
    // assign bridge-QD coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bridge-QD coupling elements in Hamiltonian.\n");
#endif
    idx2 = p->Ib + p->Nb - 1;
    for (int ii = 0; ii < p->Nc; ii++) {
      idx1 = p->Ic + ii;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
      H[idx1*N + idx2] = V[idx1][idx2];
      H[idx2*N + idx1] = V[idx2][idx1];
    }
    fprintf(stderr, "Done assigning bridge-QD coupling elements in Hamiltonian.\n");
  }
  // no bridge
  else {
    // assign bulk-QD coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bulk-QD coupling elements in Hamiltonian.\n");
#endif
    for (int ii = 0; ii < p->Nk; ii++) {
      idx1 = p->Ik + ii;
      for (int jj = 0; jj < p->Nc; jj++) {
	idx2 = p->Ic + jj;
#ifdef DEBUG_BUILDHAMILTONIAN
	fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
	fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
	fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
	H[idx1*N + idx2] = V[idx1][idx2];
	H[idx2*N + idx1] = V[idx2][idx1];
      }
    }
  }
}

/* Updates \rho(t) at each time step. */
void updateDM(N_Vector dm, realtype * dmt, int timeStep, struct PARAMETERS * p) {
#ifdef DEBUG_UPDATEDM
  std::cout << "Updating DM at step " << timeStep << "...";
#endif
  int N = 2*p->NEQ2;
  memcpy(&dmt[N*timeStep], N_VGetArrayPointer(dm), N*sizeof(realtype));
  /*
     realtype * nv = N_VGetArrayPointer(dm);
     nv_tp = &dmt[N*timeStep];
     for (int ii = 0; ii < N; ii++) {
     nv_tp[ii] = nv[ii];
     }
     */
#ifdef DEBUG_UPDATEDM
  std::cout << "done.\n";
#endif

  return;
}

/* Updates \psi(t) at each time step. */
void updateWfn(N_Vector wfn, realtype * wfnt, int timeStep, struct PARAMETERS * p) {
#ifdef DEBUG_UPDATEWFN
  std::cout << "Updating wavefunction at time step " << timeStep << "..." << std::endl;
  std::cout << "Wavefunction is " << std::endl;
  N_VPrint_Serial(wfn);
#endif
  int N = 2*p->NEQ;
  memcpy(&wfnt[N*timeStep], N_VGetArrayPointer(wfn), N*sizeof(realtype));
  /*
     for (int ii = 0; ii < N; ii++) {
     wfnt[N*timeStep + ii] = NV_Ith_S(wfnt, ii);
     }
     */
#ifdef DEBUG_UPDATEWFN
  std::cout << "done updating wavefunction." << std::endl;
#endif
  return;
}
