#ifndef __PARAMS__
#define __PARAMS__

#include <vector>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#include "constants.hpp"
#include "spline.hpp"

// Struct of parameters
struct PARAMETERS {
  int nproc = 0;			// number of processors
  bool wavefunction = 1;		// 1 => wavefunction; 0 => density matrix
  bool coherent = 1;			// 1 => \rho = \ket{\psi}\bra{\psi}; 0 => \rho_{ii} only
  bool justPlots = false;		// just make plots, no propagation or other output
  bool timedepH = true;			// if H is TD, use CVODE, else diag H and propogate
  bool analytical = false;		// turn on analytical propagation
  bool kinetic = false;			// kinetic relaxation model
  bool kineticQD = false;		// kinetic relaxation model in QD
  bool dynamicMu = false;		// for kinetic model, calculate Fermi level dynamically
  bool dynamicMuQD = false;		// for kinetic model, calculate Fermi level dynamically on QD
  bool rta = true;			// turn on relaxation time approximation (RTA)
  bool rtaQD = true;			// turn on relaxation time approximation (RTA)
  bool dephasing = false;		// turn on dephasing
  bool progressFile = false;		// create a file to show progress of the run
  bool progressStdout = false;		// show progress in stdout
  realtype abstol = 1e-10;		// absolute tolerance (for SUNDIALS)
  realtype reltol = 1e-10;		// relative tolerance (for SUNDIALS)
  realtype tout = 10000;		// final time reached by solver in atomic units
  int numsteps = 10000;			// number of time steps
  int numOutputSteps = 1000;		// number of timesteps
  int NEQ = 1;				// total number of states
  int NEQ2 = 1;				// total number of coefficients in density matrix
  realtype kBandEdge = 0.0;		// lower edge of bulk conduction band
  realtype kBandTop = 0.01;		// upper edge of bulk conduction band
  realtype lBandTop = -0.01;		// upper edge of bulk valence band
  int Nk_first = 1;			// first k state initially populated
  int Nk_final = 1;			// final k state initially populated
  realtype bulk_gap = 0.001;		// bulk band gap
  double valenceBand = 0.01;		// valence band width
  double bulkGaussSigma = 0.001;	// width of initial Gaussian in bulk
  double bulkGaussMu = 0.01;		// position of initial Gaussian above band edge
  double me = 1.0;			// effective mass of electron
  double mh = 1.0;			// effective mass of hole
  double me_c = 1.0;			// effective mass of electron on QD
  double mh_c = 1.0;			// effective mass of hole on QD
  double X2 = 1512.2873345935727;	// "Bohr radius" of material, inverse spacing in k-space
  realtype temperature = 3e2;		// temperature of the system
  double EF = 0.0;			// Fermi level in bulk
  realtype gamma1 = 1e-3;		// \gamma_1 in RTA (relaxation rate)
  realtype gamma1_c = 1e-3;		// \gamma_1 in RTA (relaxation rate) on QD
  realtype gamma2 = 1e-3;		// \gamma_2 in RTA (dephasing rate)
  double muLK = 1.0;                    // transition dipole moment from l to k (energy a.u.)
  double pumpFWHM = 1000;               // FWHM of pump pulse (time a.u.)
  double pumpPeak = 2000;               // time of peak of pump pulse (a.u.)
  double pumpFreq = 0.01;               // frequency of pump pulse (energy a.u.)
  double pumpAmpl = 1.0;                // intensity of pump pulse (electric field a.u.)
  double pumpPhase = 0.0;               // pump pulse phase (in units of radians)
  int CBPopFlag = 0;			// flag for starting condition in conduction band
  int VBPopFlag = 0;			// flag for starting condition in valence band
  int QDPopFlag = 0;			// flag for starting condition in QD
  bool bulk_FDD = false;		// switches for starting conditions
  bool bulk_Gauss = false;
  bool bulk_constant = false;
  bool qd_pops = false;
  bool laser_on = false;
  bool parabolicCoupling = false;
  bool scale_bubr = false;
  bool scale_brqd = false;
  bool scale_buqd = false;
  bool scale_laser = false;
  bool bridge_on = false;
  bool random_phase = false;
  int random_seed = 0;
  bool torsion = false;
  std::string torsionFile;
  int torsionSite = 0;

  int Nk = 1;				// number of each type of state
  int Nc = 0;				// number of each type of state
  int Nb = 0;				// number of each type of state
  int Nl = 0;				// number of each type of state
  int Ik = 0;				// index starters for each type of state
  int Ic = 0;				// index starters for each type of state
  int Ib = 0;				// index starters for each type of state
  int Il = 0;				// index starters for each type of state

  realtype kBandWidth;

  std::vector<realtype> energies;
  std::vector<realtype> Vbridge;
  std::vector<realtype> Vnobridge;
  std::vector<realtype> H;
  std::vector<realtype> H_sp;
  std::vector<int> H_cols;
  std::vector<int> H_rowind;
  std::vector<realtype> times;
  std::vector<realtype> startWfn;
  std::vector<realtype> startDM;

  std::string outputDir;

  Spline * torsionV;

  double lastTime;			// value of most recently calculated timepoint

  double lastMu;			// value of Fermi level at last time point
  double lastMuQD;			// value of Fermi level in QD at last time point
};

int bandStartIdx(int bandFlag, PARAMETERS * p);

int bandEndIdx(int bandFlag, PARAMETERS * p);

int bandNumStates(int bandFlag, PARAMETERS * p);

#endif
