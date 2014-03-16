#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <numeric>
#include <complex>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_diag.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>
#include <map>
#include <fftw3.h>
#include <omp.h>

#include "libdynamix_input_parser.hpp"
#include "libdynamix_outputs.hpp"
#include "output.hpp"
#include "numerical.hpp"
#include "params.hpp"
#include "userdata.hpp"
#include "rhs.hpp"
#include "plots.hpp"
#include "constants.hpp"
#include "conversions.hpp"
#include "analytic.hpp"

// DEBUG compiler flag: turn on to generate basic debug outputs.
// #define DEBUG

// DEBUG2 flag: turn on for more numerical output
// #define DEBUG2

// DEBUGf: debug inner CVode loop
// #define DEBUGf


int main (int argc, char * argv[]) {



  //// DECLARING VARIABLES


  // Struct of parameters
  PARAMETERS p;
  // CVode variables
  void * cvode_mem = NULL;			// pointer to block of CVode memory
  N_Vector y, yout;			// arrays of populations

  // arrays for energetic parameters
  realtype ** V = NULL;				// pointer to k-c coupling constants
  realtype * Vbridge = NULL;			// pointer to array of bridge coupling constants.
  // first element [0] is Vkb1, last [Nb] is VcbN
  realtype * Vnobridge = NULL;			// coupling constant when there is no bridge

  //// Setting defaults for parameters to be read from input
  //// done setting defaults

  int flag;
  realtype * k_pops = NULL;				// pointers to arrays of populations
  realtype * l_pops = NULL;
  realtype * c_pops = NULL;
  realtype * b_pops = NULL;
  realtype * ydata = NULL;				// pointer to ydata (contains all populations)
  realtype * wavefunction = NULL;			// (initial) wavefunction
  realtype * dm = NULL;					// density matrix
  realtype * dmt = NULL;				// density matrix in time
  realtype * wfnt = NULL;				// density matrix in time
  realtype * k_energies = NULL;				// pointers to arrays of energies
  realtype * c_energies = NULL;
  realtype * b_energies = NULL;
  realtype * l_energies = NULL;
  realtype t0 = 0.0;				// initial time
  realtype t = 0;
  realtype tret = 0;					// time returned by the solver
  time_t startRun;				// time at start of log
  time_t endRun;					// time at end of log
  struct tm * currentTime = NULL;			// time structure for localtime
#ifdef DEBUG
  FILE * realImaginary;				// file containing real and imaginary parts of the wavefunction
#endif
  FILE * log;					// log file with run times
  realtype * tkprob = NULL; 				// total probability in k, l, c, b states at each timestep
  realtype * tlprob = NULL;
  realtype * tcprob = NULL;
  realtype * tbprob = NULL;
  double ** allprob = NULL;				// populations in all states at all times
  realtype * times = NULL;
  realtype * qd_est = NULL;
  realtype * qd_est_diag = NULL;
  std::string inputFile = "ins/parameters.in";			// name of input file
  std::string cEnergiesInput = "ins/c_energies.in";
  std::string cPopsInput = "ins/c_pops.in";
  std::string bEnergiesInput = "ins/b_energies.in";
  std::string VNoBridgeInput = "ins/Vnobridge.in";
  std::string VBridgeInput = "ins/Vbridge.in";
  std::map<const std::string, bool> outs;	// map of output file names to bool

  // default output directory
  p.outputDir = "outs/";

  double summ = 0;			// sum variable

  // ---- process command line flags ---- //
  opterr = 0;
  int c;
  std::string insDir;
  /* process command line options */
  while ((c = getopt(argc, argv, "i:o:")) != -1) {
    switch (c) {
      case 'i':
	// check that it ends in a slash
	insDir = optarg;
	if (strcmp(&(insDir.at(insDir.length() - 1)), "/")) {
	  std::cerr << "ERROR: option -i requires argument ("
	            << insDir << ") to have a trailing slash (/)." << std::endl;
	  return 1;
	}
	else {
	  // ---- assign input files ---- //
	  inputFile = insDir + "parameters.in";
	  cEnergiesInput = insDir + "c_energies.in";
	  cPopsInput = insDir + "c_pops.in";
	  bEnergiesInput = insDir + "b_energies.in";
	  VNoBridgeInput = insDir + "Vnobridge.in";
	  VBridgeInput = insDir + "Vbridge.in";
	}
	break;
      case 'o':
	p.outputDir = optarg;
	break;
      case '?':
	if (optopt == 'i') {
	  fprintf(stderr, "Option -%c requires a directory argument.\n", optopt);
	}
	else if (isprint(optopt)) {
	  fprintf(stderr, "Unknown option -%c.\n", optopt);
	}
	else {
	  fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
	}
	return 1;
      default:
	continue;
    }
  }

  //// ASSIGN PARAMETERS FROM INPUT FILE


  // ---- TODO create output directory if it does not exist ---- //
  flag = mkdir(p.outputDir.c_str(), 0755);

  assignParams(inputFile.c_str(), &p);

  // Decide which output files to make
#ifdef DEBUG
  std::cout << "Assigning outputs as specified in " << inputFile << "\n";
#endif
  assignOutputs(inputFile.c_str(), outs, &p);

#ifdef DEBUG
  // print out which outputs will be made
  for (std::map<const std::string, bool>::iterator it = outs.begin(); it != outs.end(); it++) {
    std::cout << "Output file: " << it->first << " will be created.\n";
  }
#endif

  // OPEN LOG FILE; PUT IN START TIME //
  if (isOutput(outs, "log.out")) {
    log = fopen("log.out", "w");			// note that this file is closed at the end of the program
  }
  time(&startRun);
  currentTime = localtime(&startRun);
  if (isOutput(outs, "log.out")) {
    fprintf(log, "Run started at %s\n", asctime(currentTime));
  }

  if (isOutput(outs, "log.out")) {
    // make a note about the laser intensity.
    fprintf(log,"The laser intensity is %.5e W/cm^2.\n\n",pow(p.pumpAmpl,2)*3.5094452e16);
  }


  //// READ DATA FROM INPUTS


  p.Nc = numberOfValuesInFile(cEnergiesInput.c_str());
  p.Nb = numberOfValuesInFile(bEnergiesInput.c_str());
  k_pops = new realtype [p.Nk];
  c_pops = new realtype [p.Nc];
  b_pops = new realtype [p.Nb];
  l_pops = new realtype [p.Nl];
  k_energies = new realtype [p.Nk];
  c_energies = new realtype [p.Nc];
  b_energies = new realtype [p.Nb];
  l_energies = new realtype [p.Nl];
  if (numberOfValuesInFile(cPopsInput.c_str()) != p.Nc) {
    fprintf(stderr, "ERROR [Inputs]: c_pops and c_energies not the same length.\n");
    return -1;
  }
  readArrayFromFile(c_energies, cEnergiesInput.c_str(), p.Nc);
  if (p.bridge_on) {
    if (p.bridge_on && (p.Nb < 1)) {
      std::cerr << "\nERROR: bridge_on but no bridge states.  The file b_energies.in is probably empty.\n";
      return -1;
    }
    p.Vbridge.resize(p.Nb+1);
    readArrayFromFile(b_energies, bEnergiesInput.c_str(), p.Nb);
    readVectorFromFile(p.Vbridge, VBridgeInput.c_str(), p.Nb + 1);
#ifdef DEBUG
    std::cout << "COUPLINGS:";
    for (int ii = 0; ii < p.Nb+1; ii++) {
      std::cout << " " << p.Vbridge[ii];
    }
    std::cout << std::endl;
#endif
  }
  else {
    p.Nb = 0;
    p.Vnobridge.resize(1);
    readVectorFromFile(p.Vnobridge, VNoBridgeInput.c_str(), 1);
  }

#ifdef DEBUG
  std::cout << "\nDone reading things from inputs.\n";
#endif


  //// PREPROCESS DATA FROM INPUTS


  // check torsion parameters, set up torsion spline
  if (p.torsion) {
#ifdef DEBUG
    std::cout << "Torsion is on." << std::endl;
#endif

    // error checking
    if (p.torsionSite > p.Nb) {
      std::cerr << "ERROR: torsion site (" << p.torsionSite
	<< ") is larger than number of bridge sites (" << p.Nb << ")." << std::endl;
      exit(-1);
    }
    else if (p.torsionSite < 0) {
      std::cerr << "ERROR: torsion site is less than zero." << std::endl;
      exit(-1);
    }

    if (!fileExists(p.torsionFile)) {
      std::cerr << "ERROR: torsion file " << p.torsionFile << " does not exist." << std::endl;
    }

    // create spline
    p.torsionV = new Spline(p.torsionFile.c_str());
    if (p.torsionV->getFirstX() != 0.0) {
      std::cerr << "ERROR: time in " << p.torsionFile << " should start at 0.0." << std::endl;
      exit(-1);
    }
    if (p.torsionV->getLastX() < p.tout) {
      std::cerr << "ERROR: time in " << p.torsionFile << " should be >= tout." << std::endl;
      exit(-1);
    }
  }


  // set number of processors for OpenMP
  //omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);

  p.NEQ = p.Nk+p.Nc+p.Nb+p.Nl;				// total number of equations set
  p.NEQ2 = p.NEQ*p.NEQ;				// number of elements in DM
#ifdef DEBUG
  std::cout << "\nTotal number of states: " << p.NEQ << std::endl;
  std::cout << p.Nk << " bulk, " << p.Nc << " QD, " << p.Nb << " bridge, " << p.Nl << " bulk VB.\n";
#endif
  tkprob = new realtype [p.numOutputSteps+1];	// total population on k, b, c at each timestep
  tcprob = new realtype [p.numOutputSteps+1];
  tbprob = new realtype [p.numOutputSteps+1];
  tlprob = new realtype [p.numOutputSteps+1];
  allprob = new double * [p.numOutputSteps+1];
  for (int ii = 0; ii <= p.numOutputSteps; ii++) {
    allprob[ii] = new double [p.NEQ];
  }
  // assign times.
  p.times.resize(p.numOutputSteps+1);
  for (int ii = 0; ii <= p.numOutputSteps; ii++) {
    p.times[ii] = float(ii)/p.numOutputSteps*p.tout;
  }
  qd_est = new realtype [p.numOutputSteps+1];
  qd_est_diag = new realtype [p.numOutputSteps+1];
  p.Ik = 0;					// set index start positions for each type of state
  p.Ic = p.Nk;
  p.Ib = p.Ic+p.Nc;
  p.Il = p.Ib+p.Nb;

  // assign bulk conduction and valence band energies
  // for RTA, bulk and valence bands have parabolic energies
  if (p.rta) {
    buildParabolicBand(k_energies, p.Nk, p.kBandEdge, CONDUCTION, &p);
    buildParabolicBand(l_energies, p.Nl, p.lBandTop, VALENCE, &p);
  }
  else {
    buildContinuum(k_energies, p.Nk, p.kBandEdge, p.kBandTop);
    buildContinuum(l_energies, p.Nl, p.kBandEdge - p.valenceBand - p.bulk_gap, p.kBandEdge - p.bulk_gap);
  }
  // calculate band width
  p.kBandWidth = k_energies[p.Nk - 1] - k_energies[0];


  //// BUILD INITIAL WAVEFUNCTION


  // bridge states (empty to start)
  initializeArray(b_pops, p.Nb, 0.0);

  // coefficients in bulk and other states depend on input conditions in bulk
  if (!p.rta) {
#ifdef DEBUG
    std::cout << "\ninitializing k_pops\n";
#endif
    if (p.bulk_constant) {
      initializeArray(k_pops, p.Nk, 0.0);
#ifdef DEBUG
      std::cout << "\ninitializing k_pops with constant probability in range of states\n";
#endif
      initializeArray(k_pops+p.Nk_first-1, p.Nk_final-p.Nk_first+1, 1.0);
      initializeArray(l_pops, p.Nl, 0.0);		// populate l states (all 0 to start off)
      initializeArray(c_pops, p.Nc, 0.0);		// QD states empty to start
    }
    else if (p.bulk_Gauss) {
      buildKPopsGaussian(k_pops, k_energies, p.kBandEdge,
	  p.bulkGaussSigma, p.bulkGaussMu, p.Nk);   // populate k states with FDD
      initializeArray(l_pops, p.Nl, 0.0);		// populate l states (all 0 to start off)
      initializeArray(c_pops, p.Nc, 0.0);		// QD states empty to start
    }
    else if (p.qd_pops) {
      readArrayFromFile(c_pops, cPopsInput.c_str(), p.Nc);	// QD populations from file
      initializeArray(l_pops, p.Nl, 0.0);		// populate l states (all 0 to start off)
      initializeArray(k_pops, p.Nk, 0.0);             // populate k states (all zero to start off)
    }
    else {
      initializeArray(k_pops, p.Nk, 0.0);             // populate k states (all zero to start off)
      initializeArray(l_pops, p.Nl, 1.0);		// populate l states (all populated to start off)
      initializeArray(c_pops, p.Nc, 0.0);		// QD states empty to start
    }
#ifdef DEBUG
    std::cout << "\nThis is k_pops:\n";
    for (int ii = 0; ii < p.Nk; ii++) {
      std::cout << k_pops[ii] << std::endl;
    }
    std::cout << "\n";
#endif
  }
  // with RTA, use different set of switches
  else {
    // bulk valence band
    if (p.VBPopFlag == POP_EMPTY) {
#ifdef DEBUG
      std::cout << "Initializing empty valence band" << std::endl;
#endif
      initializeArray(l_pops, p.Nl, 0.0);
    }
    else if (p.VBPopFlag == POP_FULL) {
#ifdef DEBUG
      std::cout << "Initializing full valence band" << std::endl;
#endif
      initializeArray(l_pops, p.Nl, 1.0);
    }
    else {
      std::cerr << "ERROR: unrecognized VBPopFlag " << p.VBPopFlag << std::endl;
    }

    // bulk conduction band
    if (p.CBPopFlag == POP_EMPTY) {
#ifdef DEBUG
      std::cout << "Initializing empty conduction band" << std::endl;
#endif
      initializeArray(k_pops, p.Nk, 0.0);
    }
    else if (p.CBPopFlag == POP_FULL) {
#ifdef DEBUG
      std::cout << "Initializing full conduction band" << std::endl;
#endif
      initializeArray(k_pops, p.Nk, 1.0);
    }
    else if (p.CBPopFlag == POP_CONSTANT) {
#ifdef DEBUG
      std::cout << "Initializing constant distribution in conduction band" << std::endl;
#endif
      initializeArray(k_pops, p.Nk, 0.0);
      initializeArray(k_pops, p.Nk, 1e-1); // FIXME
      initializeArray(k_pops+p.Nk_first-1, p.Nk_final-p.Nk_first+1, 1.0);
    }
    else if (p.CBPopFlag == POP_GAUSSIAN) {
#ifdef DEBUG
      std::cout << "Initializing Gaussian in conduction band" << std::endl;
#endif
      buildKPopsGaussian(k_pops, k_energies, p.kBandEdge,
	  p.bulkGaussSigma, p.bulkGaussMu, p.Nk);
    }
    else {
      std::cerr << "ERROR: unrecognized CBPopFlag " << p.CBPopFlag << std::endl;
    }

    //// QD
    if (p.QDPopFlag == POP_EMPTY) {
      initializeArray(c_pops, p.Nc, 0.0);
    }
    else if (p.QDPopFlag == POP_FULL) {
      initializeArray(c_pops, p.Nc, 1.0);
    }
    else {
      std::cerr << "ERROR: unrecognized QDPopFlag " << p.QDPopFlag << std::endl;
    }
  }

  // create empty wavefunction
  wavefunction = new realtype [2*p.NEQ];
  initializeArray(wavefunction, 2*p.NEQ, 0.0);

  // assign real parts of wavefunction coefficients (imaginary are zero)
  for (int ii = 0; ii < p.Nk; ii++) {
    wavefunction[p.Ik + ii] = k_pops[ii];
  }
  for (int ii = 0; ii < p.Nc; ii++) {
    wavefunction[p.Ic + ii] = c_pops[ii];
  }
  for (int ii = 0; ii < p.Nb; ii++) {
    wavefunction[p.Ib + ii] = b_pops[ii];
  }
  for (int ii = 0; ii < p.Nl; ii++) {
    wavefunction[p.Il + ii] = l_pops[ii];
  }

  if (isOutput(outs, "psi_start.out")) {
    outputWavefunction(wavefunction, p.NEQ);
  }

  // Give all coefficients a random phase
  if (p.random_phase) {
    float phi;
    // set the seed
    if (p.random_seed == -1) { srand(time(NULL)); }
    else { srand(p.random_seed); }
    for (int ii = 0; ii < p.NEQ; ii++) {
      phi = 2*3.1415926535*(float)rand()/(float)RAND_MAX;
      wavefunction[ii] = wavefunction[ii]*cos(phi);
      wavefunction[ii + p.NEQ] = wavefunction[ii + p.NEQ]*sin(phi);
    }
  }

#ifdef DEBUG
  // print out details of wavefunction coefficients
  std::cout << std::endl;
  for (int ii = 0; ii < p.Nk; ii++) {
    std::cout << "starting wavefunction: Re[k(" << ii << ")] = " << wavefunction[p.Ik + ii] << std::endl;
  }
  for (int ii = 0; ii < p.Nc; ii++) {
    std::cout << "starting wavefunction: Re[c(" << ii << ")] = " << wavefunction[p.Ic + ii] << std::endl;
  }
  for (int ii = 0; ii < p.Nb; ii++) {
    std::cout << "starting wavefunction: Re[b(" << ii << ")] = " << wavefunction[p.Ib + ii] << std::endl;
  }
  for (int ii = 0; ii < p.Nl; ii++) {
    std::cout << "starting wavefunction: Re[l(" << ii << ")] = " << wavefunction[p.Il + ii] << std::endl;
  }
  for (int ii = 0; ii < p.Nk; ii++) {
    std::cout << "starting wavefunction: Im[k(" << ii << ")] = " << wavefunction[p.Ik + ii + p.NEQ] << std::endl;
  }
  for (int ii = 0; ii < p.Nc; ii++) {
    std::cout << "starting wavefunction: Im[c(" << ii << ")] = " << wavefunction[p.Ic + ii + p.NEQ] << std::endl;
  }
  for (int ii = 0; ii < p.Nb; ii++) {
    std::cout << "starting wavefunction: Im[b(" << ii << ")] = " << wavefunction[p.Ib + ii + p.NEQ] << std::endl;
  }
  for (int ii = 0; ii < p.Nl; ii++) {
    std::cout << "starting wavefunction: Im[l(" << ii << ")] = " << wavefunction[p.Il + ii + p.NEQ] << std::endl;
  }
  std::cout << std::endl;
  summ = 0;
  for (int ii = 0; ii < 2*p.NEQ; ii++) {
    summ += pow(wavefunction[ii],2);
  }
  std::cout << "\nTotal population is " << summ << "\n\n";
#endif


  //// ASSEMBLE ARRAY OF ENERGIES


  // TODO TODO
  p.energies.resize(p.NEQ);
  for (int ii = 0; ii < p.Nk; ii++) {
    p.energies[p.Ik + ii] = k_energies[ii];
  }
  for (int ii = 0; ii < p.Nc; ii++) {
    p.energies[p.Ic + ii] = c_energies[ii];
  }
  for (int ii = 0; ii < p.Nb; ii++) {
    p.energies[p.Ib + ii] = b_energies[ii];
  }
  for (int ii = 0; ii < p.Nl; ii++) {
    p.energies[p.Il + ii] = l_energies[ii];
  }

#ifdef DEBUG
  for (int ii = 0; ii < p.NEQ; ii++) {
    std::cout << "p.energies[" << ii << "] is " << p.energies[ii] << "\n";
  }
#endif


  //// ASSIGN COUPLING CONSTANTS


  V = new realtype * [p.NEQ];
  for (int ii = 0; ii < p.NEQ; ii++) {
    V[ii] = new realtype [p.NEQ];
  }
  buildCoupling(V, &p, outs);

  if (isOutput(outs, "log.out")) {
    // make a note in the log about system timescales
    double tau = 0;		// fundamental system timescale
    if (p.Nk == 1) {
      fprintf(log, "\nThe timescale (tau) is undefined (Nk == 1).\n");
    }
    else {
      if (p.bridge_on) {
	if (p.scale_bubr) {
	  tau = 1.0/(2*p.Vbridge[0]*M_PI);
	}
	else {
	  tau = ((p.kBandTop - p.kBandEdge)/(p.Nk - 1))/(2*pow(p.Vbridge[0],2)*M_PI);
	}
      }
      else {
	if (p.scale_buqd) {
	  tau = 1.0/(2*p.Vnobridge[0]*M_PI);
	}
	else {
	  tau = ((p.kBandTop - p.kBandEdge)/(p.Nk - 1))/(2*pow(p.Vnobridge[0],2)*M_PI);
	}
      }
      fprintf(log, "\nThe timescale (tau) is %.9e a.u.\n", tau);
    }
  }

  //// CREATE DENSITY MATRIX
  if (! p.wavefunction) {
    // Create the initial density matrix
    dm = new realtype [2*p.NEQ2];
    initializeArray(dm, 2*p.NEQ2, 0.0);
#pragma omp parallel for
    for (int ii = 0; ii < p.NEQ; ii++) {
      // diagonal part
      dm[p.NEQ*ii + ii] = pow(wavefunction[ii],2) + pow(wavefunction[ii + p.NEQ],2);
      if (p.coherent) {
	// off-diagonal part
	for (int jj = 0; jj < ii; jj++) {
	  // real part of \rho_{ii,jj}
	  dm[p.NEQ*ii + jj] = wavefunction[ii]*wavefunction[jj] + wavefunction[ii+p.NEQ]*wavefunction[jj+p.NEQ];
	  // imaginary part of \rho_{ii,jj}
	  dm[p.NEQ*ii + jj + p.NEQ2] = wavefunction[ii]*wavefunction[jj+p.NEQ] - wavefunction[jj]*wavefunction[ii+p.NEQ];
	  // real part of \rho_{jj,ii}
	  dm[p.NEQ*jj + ii] = dm[p.NEQ*ii + jj];
	  // imaginary part of \rho_{jj,ii}
	  dm[p.NEQ*jj + ii + p.NEQ2] = -1*dm[p.NEQ*ii + jj + p.NEQ*p.NEQ];
	}
      }
    }

    // Create the array to store the density matrix in time
    dmt = new realtype [2*p.NEQ2*(p.numOutputSteps+1)];
    initializeArray(dmt, 2*p.NEQ2*(p.numOutputSteps+1), 0.0);

#ifdef DEBUG2
    // print out density matrix
    std::cout << "\nDensity matrix without normalization:\n\n";
    for (int ii = 0; ii < p.NEQ; ii++) {
      for (int jj = 0; jj < p.NEQ; jj++) {
	fprintf(stdout, "(%+.1e,%+.1e) ", dm[p.NEQ*ii + jj], dm[p.NEQ*ii + jj + p.NEQ2]);
      }
      fprintf(stdout, "\n");
    }
#endif

    // Normalize the DM so that populations add up to 1.
    // No normalization if RTA is on.
    if (!p.rta) {
      summ = 0.0;
      for (int ii = 0; ii < p.NEQ; ii++) {
	// assume here that diagonal elements are all real
	summ += dm[p.NEQ*ii + ii];
      }
      if ( summ == 0.0 ) {
	std::cerr << "\nFATAL ERROR [populations]: total population is 0!\n";
	return -1;
      }
      if (summ != 1.0) {
	// the variable 'summ' is now a multiplicative normalization factor
	summ = 1.0/summ;
	for (int ii = 0; ii < 2*p.NEQ2; ii++) {
	  dm[ii] *= summ;
	}
      }
#ifdef DEBUG
      std::cout << "\nThe normalization factor for the density matrix is " << summ << "\n\n";
#endif
    }

    // Error checking for total population; recount population first
    summ = 0.0;
    for (int ii = 0; ii < p.NEQ; ii++) {
      summ += dm[p.NEQ*ii + ii];
    }
    if ( fabs(summ-1.0) > 1e-12  && (!p.rta)) {
      std::cerr << "\nWARNING [populations]: After normalization, total population is not 1, it is " << summ << "!\n";
    }
#ifdef DEBUG
    std::cout << "\nAfter normalization, the sum of the populations in the density matrix is " << summ << "\n\n";
#endif
    // Add initial DM to parameters.
    p.startDM.resize(2*p.NEQ2);
    memcpy(&(p.startDM[0]), &(dm[0]), 2*p.NEQ2*sizeof(double));
  }
  // wavefunction
  else {

    // Create the array to store the wavefunction in time
    wfnt = new realtype [2*p.NEQ*(p.numOutputSteps+1)];
    initializeArray(wfnt, 2*p.NEQ*(p.numOutputSteps+1), 0.0);

    // normalize
    summ = 0.0;
    for (int ii = 0; ii < p.NEQ; ii++) {
      summ += pow(wavefunction[ii],2) + pow(wavefunction[ii+p.NEQ],2);
    }
#ifdef DEBUG
    std::cout << "Before normalization, the total population is " << summ << std::endl;
#endif
    summ = 1.0/sqrt(summ);
    for (int ii = 0; ii < 2*p.NEQ; ii++) {
      wavefunction[ii] *= summ;
    }

    // check total population
    summ = 0.0;
    for (int ii = 0; ii < p.NEQ; ii++) {
      summ += pow(wavefunction[ii],2) + pow(wavefunction[ii+p.NEQ],2);
    }
#ifdef DEBUG
    std::cout << "After normalization, the total population is " << summ << std::endl;
#endif
    if (fabs(summ - 1.0) > 1e-12) {
      std::cerr << "WARNING: wavefunction not normalized!  Total density is " << summ << std::endl;
    }

    // Add initial wavefunction to parameters.
    p.startWfn.resize(2*p.NEQ);
    memcpy(&(p.startWfn[0]), &(wavefunction[0]), 2*p.NEQ*sizeof(double));
  }


  //// BUILD HAMILTONIAN


  // //TODO TODO
#ifdef DEBUG
  fprintf(stderr, "Building Hamiltonian.\n");
#endif
  realtype * H = NULL;
  H = new realtype [p.NEQ2];
  for (int ii = 0; ii < p.NEQ2; ii++) {
    H[ii] = 0.0;
  }
  buildHamiltonian(H, p.energies, V, &p);
  // add Hamiltonian to p
  p.H.resize(p.NEQ2);
  for (int ii = 0; ii < p.NEQ2; ii++) {
    p.H[ii] = H[ii];
  }
  // create sparse version of H
  p.H_sp.resize(p.NEQ2);
  p.H_cols.resize(p.NEQ2);
  p.H_rowind.resize(p.NEQ2 + 1);
  int job [6] = {0, 0, 0, 2, p.NEQ2, 1};
  int info = 0;

  mkl_ddnscsr(&job[0], &(p.NEQ), &(p.NEQ), &(p.H)[0], &(p.NEQ), &(p.H_sp)[0],
      &(p.H_cols)[0], &(p.H_rowind)[0], &info);


  //// SET UP CVODE VARIABLES


#ifdef DEBUG
  std::cout << "\nCreating N_Vectors.\n";
  if (p.wavefunction) {
    std::cout << "\nProblem size is " << 2*p.NEQ << " elements.\n";
  }
  else {
    std::cout << "\nProblem size is " << 2*p.NEQ2 << " elements.\n";
  }
#endif
  // Creates N_Vector y with initial populations which will be used by CVode//
  if (p.wavefunction) {
    y = N_VMake_Serial(2*p.NEQ, wavefunction);
  }
  else {
    y = N_VMake_Serial(2*p.NEQ2, dm);
  }
  // put in t = 0 information
  if (! p.wavefunction) {
    updateDM(y, dmt, 0, &p);
  }
  else {
    updateWfn(y, wfnt, 0, &p);
  }
  // the vector yout has the same dimensions as y
  yout = N_VClone(y);

#ifdef DEBUG
  realImaginary = fopen("real_imaginary.out", "w");
#endif

  // Make plot files
  makePlots(outs, &p);

  // only do propagation if not just making plots
  if (! p.justPlots) {
    // Make outputs independent of time propagation
    computeGeneralOutputs(outs, &p);

    // create CVode object
    // this is a stiff problem, I guess?
#ifdef DEBUG
    std::cout << "\nCreating cvode_mem object.\n";
#endif
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    flag = CVodeSetUserData(cvode_mem, (void *) &p);

#ifdef DEBUG
    std::cout << "\nInitializing CVode solver.\n";
#endif
    // initialize CVode solver //

    if (p.wavefunction) {
      //flag = CVodeInit(cvode_mem, &RHS_WFN, t0, y);
      flag = CVodeInit(cvode_mem, &RHS_WFN_SPARSE, t0, y);
    }
    else {
      if (p.kinetic) {
	flag = CVodeInit(cvode_mem, &RHS_DM_RELAX, t0, y);
      }
      else if (p.rta) {
	flag = CVodeInit(cvode_mem, &RHS_DM_RTA, t0, y);
	//flag = CVodeInit(cvode_mem, &RHS_DM_RTA_BLAS, t0, y);
      }
      else if (p.dephasing) {
	flag = CVodeInit(cvode_mem, &RHS_DM_dephasing, t0, y);
      }
      else {
	//flag = CVodeInit(cvode_mem, &RHS_DM, t0, y);
	flag = CVodeInit(cvode_mem, &RHS_DM_BLAS, t0, y);
      }
    }

#ifdef DEBUG
    std::cout << "\nSpecifying integration tolerances.\n";
#endif
    // specify integration tolerances //
    flag = CVodeSStolerances(cvode_mem, p.reltol, p.abstol);

#ifdef DEBUG
    std::cout << "\nAttaching linear solver module.\n";
#endif
    // attach linear solver module //
    if (p.wavefunction) {
      flag = CVDense(cvode_mem, 2*p.NEQ);
    }
    else {
      // Diagonal approximation to the Jacobian saves memory for large systems
      flag = CVDiag(cvode_mem);
    }

    //// CVODE TIME PROPAGATION


#ifdef DEBUG
    std::cout << "\nAdvancing the solution in time.\n";
#endif
    for (int ii = 1; ii <= p.numsteps; ii++) {
      t = (p.tout*((double) ii)/((double) p.numsteps));
      flag = CVode(cvode_mem, t, yout, &tret, 1);
#ifdef DEBUGf
      std::cout << std::endl << "CVode flag at step " << ii << ": " << flag << std::endl;
#endif
      if ((ii % (p.numsteps/p.numOutputSteps) == 0) || (ii == p.numsteps)) {
	// show progress in stdout
	if (p.progressStdout) {
	  fprintf(stdout, "\r%-.2lf percent done", ((double)ii/((double)p.numsteps))*100);
	  fflush(stdout);
	}
	// show progress in a file
	if (p.progressFile) {
	  std::ofstream progressFile("progress.tmp");
	  progressFile << ((double)ii/((double)p.numsteps))*100 << " percent done." << std::endl;
	  progressFile.close();
	}
	if (p.wavefunction) {
	  updateWfn(yout, wfnt, ii*p.numOutputSteps/p.numsteps, &p);
	}
	else {
	  updateDM(yout, dmt, ii*p.numOutputSteps/p.numsteps, &p);
	}
      }
    }

#ifdef DEBUG
    fclose(realImaginary);
#endif


    //// MAKE FINAL OUTPUTS


    // finalize log file //
    time(&endRun);
    currentTime = localtime(&endRun);
    if (isOutput(outs, "log.out")) {
      fprintf(log, "Final status of 'flag' variable: %d\n\n", flag);
      fprintf(log, "Run ended at %s\n", asctime(currentTime));
      fprintf(log, "Run took %.3g seconds.\n", difftime(endRun, startRun));
      fclose(log);					// note that the log file is opened after variable declaration
    }
    if (p.progressStdout) {
      printf("\nRun took %.3g seconds.\n", difftime(endRun, startRun));
    }

    // Compute density outputs.
#ifdef DEBUG
    std::cout << "Computing outputs..." << std::endl;
#endif
    if (p.wavefunction) {
      computeWfnOutput(wfnt, outs, &p);
    }
    else {
      computeDMOutput(dmt, outs, &p);
    }
#ifdef DEBUG
    std::cout << "done computing outputs" << std::endl;
#endif

    // do analytical propagation
    if (p.analytical && (! p.bridge_on)) {
      computeAnalyticOutputs(outs, &p);
    }
  }


  //// CLEAN UP


#ifdef DEBUG
  fprintf(stdout, "Deallocating N_Vectors.\n");
#endif
  // deallocate memory for N_Vectors //
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(yout);

#ifdef DEBUG
  fprintf(stdout, "Freeing CVode memory.\n");
#endif
  // free solver memory //
  CVodeFree(&cvode_mem);

#ifdef DEBUG
  fprintf(stdout, "Freeing memory in main.\n");
#endif
  // delete all these guys
  delete [] tkprob;
  delete [] tlprob;
  delete [] tcprob;
  delete [] tbprob;
  for (int ii = 0; ii <= p.numOutputSteps; ii++) {
    delete [] allprob[ii];
  }
  delete [] allprob;
  delete [] k_pops;
  delete [] c_pops;
  delete [] b_pops;
  delete [] l_pops;
  if (p.bridge_on) {
    delete [] Vbridge;
  }
  else {
    delete [] Vnobridge;
  }
  delete [] k_energies;
  delete [] c_energies;
  delete [] b_energies;
  delete [] l_energies;
  delete [] wavefunction;
  delete [] H;
  for (int ii = 0; ii < p.NEQ; ii++) {
    delete [] V[ii];
  }
  delete [] V;
  if (p.wavefunction) {
    delete [] wfnt;
  }
  else {
    delete [] dm;
    delete [] dmt;
  }
  delete [] times;
  delete [] qd_est;
  delete [] qd_est_diag;

  std::cout << "whoo" << std::endl;

  return 0;
}
