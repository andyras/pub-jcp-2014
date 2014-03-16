#include "rhs.hpp"

// #define DEBUG_RHS
//#define DEBUG_RTA
//
// DEBUGf flag: general output at each CVode step
//#define DEBUGf
//
// DEBUGf_DM flag: DEBUGf for density matrix EOM
// #define DEBUGf_DM

//#define DEBUG_TORSION
//#define DEBUG_DYNAMIC_MU


/* Right-hand-side equation for wavefunction */
int RHS_WFN(realtype t, N_Vector y, N_Vector ydot, void * data) {
  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // extract parameters from p
  realtype * H = &(p->H)[0];
  //realtype * H = &(p->H_lo)[0];
  int N = p->NEQ;

  // get pointer to y, ydot data
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // set up BLAS variables
  const char TRANS = 'n';
  const char UPLO = 'l';
  double beta = 0.0;
  double alpha_re = 1.0;	// alpha value for real part of wfn derivative
  double alpha_im = -1.0;	// alpha value for imag part of wfn derivative
  int inc = 1;

  // Re(\dot{\psi}) = \hat{H}Im(\psi)
  //DGEMV(&TRANS, &N, &N, &alpha_re, &H[0], &N, &yp[N], &inc, &beta, &ydotp[0], &inc);
  DSYMV(&UPLO, &N, &alpha_re, &H[0], &N, &yp[N], &inc, &beta, &ydotp[0], &inc);

  // Im(\dot{\psi}) = -i\hat{H}Re(\psi)
  //DGEMV(&TRANS, &N, &N, &alpha_im, &H[0], &N, &yp[0], &inc, &beta, &ydotp[N], &inc);
  DSYMV(&UPLO, &N, &alpha_im, &H[0], &N, &yp[0], &inc, &beta, &ydotp[N], &inc);

  return 0;
}

int RHS_WFN_SPARSE(realtype t, N_Vector y, N_Vector ydot, void * data) {
  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // extract parameters from p
  realtype * H = &(p->H_sp)[0];
  int * columns = &(p->H_cols)[0];
  int * rowind = &(p->H_rowind)[0];

  //realtype * H = &(p->H_lo)[0];
  int N = p->NEQ;

  // set up MKL variables
  char transa = 'n';
  double alpha_re = 1.0;	// alpha value for real part of wfn derivative
  double alpha_im = -1.0;	// alpha value for imag part of wfn derivative
  double beta = 0.0;
  char matdescra [6] = {'S', // symmetric matrix
			'L', // lower triangle
			'N', // non-unit on diagonal
			'C', // zero-based indexing (C-style)
			'*', '*'}; // extra characters

  mkl_set_num_threads(1);

  // Re(\dot{\psi}) = \hat{H}Im(\psi)
  mkl_dcsrmv(&transa, &N, &N, &alpha_re, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[N], &beta, &ydotp[0]);

  // Im(\dot{\psi}) = -i\hat{H}Re(\psi)
  mkl_dcsrmv(&transa, &N, &N, &alpha_im, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[0], &beta, &ydotp[N]);

  return 0;
}

/* apply the kinetic relaxation model to one band of the system */
void RELAX_KINETIC(int bandFlag, realtype * yp, realtype * ydotp, PARAMETERS * p) {
  double pop = 0.0;
  int start = bandStartIdx(bandFlag, p);
  int end = bandEndIdx(bandFlag, p);
  int Ni = bandNumStates(bandFlag, p);
  int N = p->NEQ;
  double g1;
  if (bandFlag == CONDUCTION) {
    g1 = p->gamma1;
  }
  else if (bandFlag == QD_CONDUCTION) {
    g1 = p->gamma1_c;
  }
  else {
    std::cerr << "WARNING [" << __FUNCTION__ << "]: unexpected bandFlag " << bandFlag << std::endl;
    std::cerr << "        setting g1 to 1.0" << std::endl;
    g1 = 1.0;
  }
  double g2 = p->gamma2;
  double mu = p->EF;
  double T = p->temperature;
  double * E = &(p->energies[start]);

  // sum current populations in band
  pop = 0.0;
  for (int ii = start; ii < end; ii++) {
    pop += yp[ii*N + ii];
  }

  // do the (N-1) pairs of states along diagonal
  double ePi, ePj;	// equilibrium populations, i and j indices
  double rel;		// relaxation term
  int Ii, Ij;		// indices

  // find equilibrium FDD
  double * fdd = new double [Ni];

  if ((p->dynamicMu) && (pop > 0.0)) {
    if (pop > 1.001) {	// test is against 1.001 since there may be some numerical drift
      std::cout << "WARNING [" << __FUNCTION__
	<< "]: population in band is > 1; dynamic Fermi level may be spurious" << std::endl;
    }

    //// find bounds for Fermi level
    mu = findDynamicMu(pop, T, CONDUCTION, p);
#ifdef DEBUG_DYNAMIC_MU
    std::cout << "mu at time " << t << " is " << mu << std::endl;
#endif
  }

  FDD(mu, T, fdd, E, Ni, pop);

  for (int ii = 0; ii < (Ni-1); ii++) {
    // precalculate indices and such
    Ii = (start + ii)*N + start + ii;
    Ij = Ii + N + 1;		// this index is the next diagonal element, so N+1 places up
    ePi = fdd[ii];
    ePj = fdd[ii+1];

    //// calculate contribution from relaxation

    // assuming \Gamma = k_{ij} + k_{ji}, "unimolecular" model
    rel = g1*(ePi*yp[Ij] - ePj*yp[Ii])/(ePi + ePj);

    // assuming \Gamma = k_{ij} + k_{ji}, "bimolecular" model
    // rel = g1*(yp[Ij]*(1-yp[Ii])*ePi*(1-ePj) - yp[Ii]*(1-yp[Ij])*ePj*(1-ePi))/(ePi+ePj-2*ePi*ePj);

    // assuming downward rates (j-->i) are the same, "bimolecular" model
    // rel = g1*(yp[Ij]*(1 - yp[Ii]) - ePj*(1-ePi)/(ePi*(1-ePj))*yp[Ii]*(1-yp[Ij]));

    // equal and opposite for the interaction of the two states
    ydotp[Ii] += rel;
    ydotp[Ij] -= rel;
  }

  delete [] fdd;

  return;
}

void RELAX_RTA(int bandFlag, realtype * yp, realtype * ydotp, PARAMETERS * p) {
  return;
}

/* Right-hand-side equation for density matrix */
int RHS_DM_RELAX(realtype t, N_Vector y, N_Vector ydot, void * data) {

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  double g2 = p->gamma2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// relaxation in bulk conduction band

  if (p->kinetic) {
    RELAX_KINETIC(CONDUCTION, yp, ydotp, p);
  }
  else if (p->rta) {
    RELAX_RTA(CONDUCTION, yp, ydotp, p);
  }

  //// relaxation in QD conduction band

  if (p->kineticQD) {
    RELAX_KINETIC(QD_CONDUCTION, yp, ydotp, p);
  }
  else if (p->rtaQD) {
    RELAX_RTA(QD_CONDUCTION, yp, ydotp, p);
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
	ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
	ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];

      // dephasing
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  dmf = fopen("dmf.out", "a");
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  fclose(dmf);
#endif

  return 0;
}

/* Right-hand-side equation for density matrix */
int RHS_DM_KINETIC(realtype t, N_Vector y, N_Vector ydot, void * data) {

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  int Nk = p->Nk;
  int Ik = p->Ik;
  double * E = &(p->energies[Ik]);
  double g1 = p->gamma1;
  double g2 = p->gamma2;
  double T = p->temperature;
  double mu = p->EF;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// relaxation along diagonal
  // sum current populations in band
  double CBPop = 0.0;
  for (int ii = Ik; ii < (Ik + Nk); ii++) {
    CBPop += yp[ii*N + ii];
  }

  // do the (N-1) pairs of states along diagonal
  double ePi, ePj;	// equilibrium populations, i and j indices
  double rel;		// relaxation term
  int Ii, Ij;		// indices

  //// Kinetic relaxation model here

  // find equilibrium FDD
  double * fdd = new double [Nk];

  if ((p->dynamicMu) && (CBPop > 0.0)) {
    if (CBPop > 1.001) {	// test is against 1.001 since there may be some numerical drift
      std::cout << "WARNING [" << __FUNCTION__
	<< "]: population in band is > 1; dynamic Fermi level may be spurious" << std::endl;
    }

    //// find bounds for Fermi level
    mu = findDynamicMu(CBPop, T, CONDUCTION, p);
#ifdef DEBUG_DYNAMIC_MU
    std::cout << "mu at time " << t << " is " << mu << std::endl;
#endif
  }

  FDD(mu, T, fdd, E, Nk, CBPop);

  for (int ii = 0; ii < (Nk-1); ii++) {
    // precalculate indices and such
    Ii = (Ik + ii)*N + Ik + ii;
    Ij = Ii + N + 1;		// this index is the next diagonal element, so N+1 places up
    ePi = fdd[ii];
    ePj = fdd[ii+1];

    // calculate contribution from relaxation
    rel = g1*(ePi*yp[Ij] - ePj*yp[Ii])/(ePi + ePj);

    // equal and opposite for the interaction of the two states
    ydotp[Ii] += rel;
    ydotp[Ij] -= rel;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
	ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
	ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];

      // relaxation
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  dmf = fopen("dmf.out", "a");
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  fclose(dmf);
#endif

  // free fdd
  delete [] fdd;

  return 0;
}

/* Right-hand-side equation for density matrix */
int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
	ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
	ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  std::cout << "Closing output file for density matrix coefficients in time.\n";
  fclose(dmf);
#endif

  return 0;
}

/* Right-hand-side equation for density matrix using BLAS */
int RHS_DM_BLAS(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  double * H = &(p->H)[0];
  int N = p->NEQ;
  int N2 = p->NEQ2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot

#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  char TRANSA = 'n';
  char TRANSB = 'n';
  char LEFT = 'l';
  char RGHT = 'r';
  char UPLO = 'l';
  double ONE = 1.0;
  double NEG = -1.0;

  realtype * H_sp = &(p->H_sp)[0];
  int * columns = &(p->H_cols)[0];
  int * rowind = &(p->H_rowind)[0];

  // set up MKL variables
  double beta = 0.0;
  char matdescra [6] = {'T', // symmetric matrix
			'L', // lower triangle
			'N', // non-unit on diagonal
			'C', // zero-based indexing (C-style)
			'*', '*'};

			// set beta to zero for first call

  mkl_set_num_threads(p->nproc);
  // Re(\dot{\rho}) += H*Im(\rho)
  //DGEMM(&TRANSA, &TRANSB, &N, &N, &N, &ONE, &H[0], &N, &yp[N2], &N, &ONE, &ydotp[0], &N);
  //DSYMM(&LEFT, &UPLO, &N, &N, &ONE, &H[0], &N, &yp[N2], &N, &ONE, &ydotp[0], &N);
  //N_VPrint_Serial(ydot);
  // Re(\dot{\psi}) = \hat{H}Im(\psi)
  mkl_dcsrmm(&TRANSA, &N, &N, &N, &ONE, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[N2], &N, &beta, &ydotp[0], &N);
  
  // Re(\dot{\rho}) -= Im(\rho)*H
  //DGEMM(&TRANSA, &TRANSB, &N, &N, &N, &NEG, &yp[N2], &N, &H[0], &N, &ONE, &ydotp[0], &N);
  DSYMM(&RGHT, &UPLO, &N, &N, &NEG, &H[0], &N, &yp[N2], &N, &ONE, &ydotp[0], &N);

  //N_VPrint_Serial(ydot);
  // Im(\dot{\rho}) += i*Re(\rho)*H
  //DGEMM(&TRANSA, &TRANSB, &N, &N, &N, &ONE, &yp[0], &N, &H[0], &N, &ONE, &ydotp[N2], &N);
  DSYMM(&RGHT, &UPLO, &N, &N, &ONE, &H[0], &N, &yp[0], &N, &ONE, &ydotp[N2], &N);

  //N_VPrint_Serial(ydot);
  // Im(\dot{\rho}) -= i*H*Re(\rho)
  //DGEMM(&TRANSA, &TRANSB, &N, &N, &N, &NEG, &H[0], &N, &yp[0], &N, &ONE, &ydotp[N2], &N);
  //DSYMM(&LEFT, &UPLO, &N, &N, &NEG, &H[0], &N, &yp[0], &N, &ONE, &ydotp[N2], &N);
  //N_VPrint_Serial(ydot);

  // Im(\dot{\psi}) = -i\hat{H}Re(\psi)
  mkl_dcsrmm(&TRANSA, &N, &N, &N, &NEG, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[0], &N, &beta, &ydotp[N2], &N);

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  std::cout << "Closing output file for density matrix coefficients in time.\n";
  fclose(dmf);
#endif

  return 0;
}

/* find Fermi level based on sum of population in band */
double findDynamicMu(double pop, double T, int bandFlag, PARAMETERS * p) {
  double summ, lower, upper;

  // starting point
  double mu = 0.0;
  if (bandFlag == CONDUCTION) {
    mu = p->lastMu;
  }
  else if (bandFlag == QD_CONDUCTION) {
    mu = p->lastMuQD;
  }

  summ = FDDSum(mu, T, bandFlag, p);

  // just in case mu is exactly zero
  if (fabs((summ - pop) > 1e-10)) {
    // otherwise search in increments of 1 a.u. energy
    if (summ < pop) {
      lower = -1.0;
      upper = 0.0;
      while (summ < pop) {
	lower++;
	upper++;
#ifdef DEBUG_DYNAMIC_MU
	std::cout << "Lower bound for mu: " << lower << std::endl;
	std::cout << "Upper bound for mu: " << upper << std::endl;
#endif
	summ = FDDSum(upper, T, bandFlag, p);
      }
    }
    else {
      lower = 0.0;
      upper = 1.0;
      while (summ > pop) {
	lower--;
	upper--;
#ifdef DEBUG_DYNAMIC_MU
	std::cout << "Lower bound for mu: " << lower << std::endl;
	std::cout << "Upper bound for mu: " << upper << std::endl;
#endif
	summ = FDDSum(lower, T, bandFlag, p);
      }
    }

    // do a binary search for mu
    mu = FDDBinarySearch(lower, upper, T, pop, bandFlag, p);
  }

  // store the mu for the next time step
  if (bandFlag == CONDUCTION) {
    p->lastMu = mu;
  }
  else if (bandFlag == QD_CONDUCTION) {
    p->lastMuQD = mu;
  }

  return mu;
}

/* Do a binary search to find the value of the Fermi level which makes the
 * sum of populations in a band add up to a certain value.
 */
double FDDBinarySearch(double lower, double upper, double T, double n,
    int bandFlag, PARAMETERS * p) {
  double mid, summ;

  if (fabs(upper - lower) < 1e-10) {
    return lower;
  }
  else {
    mid = (upper + lower)/2.0;
    summ = FDDSum(mid, T, bandFlag, p);
#ifdef DEBUG_DYNAMIC_MU
    std::cout << "Binary search lower bound: " << lower << std::endl;
    std::cout << "Binary search upper bound: " << upper << std::endl;
    std::cout << "Binary search middle: " << mid << std::endl;
    std::cout << "Binary search target value: " << n << std::endl;
    std::cout << "Binary search summ value: " << summ << std::endl;
#endif
    if (summ > n) {
#ifdef DEBUG_DYNAMIC_MU
      std::cout << "value is in lower half of bounds" << std::endl;
      std::cout << std::endl;
#endif
      return FDDBinarySearch(lower, mid, T, n, bandFlag, p);
    }
    else {
#ifdef DEBUG_DYNAMIC_MU
      std::cout << "value is in upper half of bounds" << std::endl;
      std::cout << std::endl;
#endif
      return FDDBinarySearch(mid, upper, T, n, bandFlag, p);
    }
  }
}

/* Add up the populations in a band with a Fermi-Dirac distribution of population
 */
double FDDSum(double mu, double T, int bandFlag, PARAMETERS * p) {
  double summ = 0.0;
  int start = bandStartIdx(bandFlag, p);
  int end = bandEndIdx(bandFlag, p);
  double beta = 3.185e5/T;
  double * E = &(p->energies[start]);

  for (int ii = start; ii < end; ii++) {
    summ += 1.0/(1.0 + exp((E[ii] - mu)*beta));
  }

  return summ;
}

/* fills the array fdd with Fermi-Dirac populations, normalized to a population
 * P.
 */
void FDD(double mu, double T, double * fdd, double * E, int N, double P) {
  double beta = 3.185e5/T;

  for (int ii = 0; ii < N; ii++) {
    fdd[ii] = 1.0/(1.0 + exp((E[ii] - mu)*beta));
  }

  if (P != 0.0) {
    double norm = P/std::accumulate(fdd, fdd+N, 0.0);

    for (int ii = 0; ii < N; ii++) {
      fdd[ii] *= norm;
    }
  }

  return;
}


/* gives the equilibrated FDD for the system */
void FDD_RTA(struct PARAMETERS * p, realtype * y, double * fdd, int flag) {
#ifdef DEBUG_RTA
  //// "fine structure constant" -- conversion from index to wave vector
  std::cout << "p->X2   " << p->X2 << std::endl;
#endif
  // unpack parameters
  int N = p->NEQ;
  int Ni;		// number of states
  int Ii;		// index of states
  if (flag == 1) {	// FDD in conduction band
    Ni = p->Nk;
    Ii = p->Ik;
  }
  else {		// FDD in QD
    Ni = p->Nc;
    Ii = p->Ic;
  }
  double me;		// electron effective mass
  if (flag == 1) {
    me = p->me;
  }
  else {
    me = p->me_c;
  }
  double X2 = p->X2;

  // FDD is unity when number of states is 1
  if (Ni < 2) {
    fdd[0] = 1.0;
    return;
  }

  //// calculate n_e and e_kin
  double ne = 0.0;
  double ekin = 0.0;
  double factor = 1.0/(M_PI*M_PI*pow(X2,3));
#ifdef DEBUG_RTA
  std::cout << "factor   " << factor << std::endl;
#endif
  // assign vector of energies
  std::vector<double> Ev (Ni,0.0);
  realtype * E = &(Ev[0]);		// TODO don't need this if have Intel/GCC flag
  for (int ii = 0; ii < Ni; ii++) {
    E[ii] = pow(ii,2)/(2*me*pow(X2,2));
  }

#ifdef DEBUG_RTA
  std::cout << std::setprecision(28);
#endif
  // Simpson's Rule method
  // skip the first point because the value will be zero
  double SF = 4.0;	// Simpson's factor
  int sign = -1;	// sign
  for (int ii = 1; ii < (Ni-1); ii++) {
    ne += SF*factor*ii*ii*y[(Ii + ii)*N + (Ii + ii)];
    ekin += SF*factor*pow(ii,2)*y[(Ii + ii)*N + (Ii + ii)]*E[ii];
#ifdef DEBUG_RTA
    std::cout << "Ne " << ii*ii << "*" << SF << "/3.0*" << y[ii*N + ii] << "/" << pow(5.29e-11,3)/factor << std::endl;
    std::cout << "ekin " << pow(ii,4) << "*" << SF << "/3.0*" << y[ii*N + ii] << "*" << 4.3597482e-18/(2*me*pow(X2,2)) << "/" << pow(5.29e-11,3)/factor << std::endl;
    std::cout << "ekin " << ekin/pow(5.29e-11,3)*4.3597482e-18/3.0 
      << " += " << SF*factor*pow(ii,4)*y[ii*N + ii]/(2*me*pow(X2,2))/pow(5.29e-11,3)*4.3597482e-18/3.0 << std::endl;
#endif
    SF += sign*2.0;
    sign *= -1;
  }
  // add last point
  ne += factor*pow(Ni-1,2)*y[(Ii + Ni - 1)*N + Ii + Ni - 1];
  ekin += factor*pow(Ni-1,2)*y[(Ii + Ni - 1)*N + Ii + Ni - 1]*E[Ni - 1];
  // divide by three
  ne /= 3.0;
  ekin /= 3.0;

#ifdef DEBUG_RTA
  std::cout << "ne        " << ne << std::endl;
  std::cout << "ne (SI)   " << ne/pow(5.29e-11,3) << std::endl;
  std::cout << "ekin      " << ekin << std::endl;
  std::cout << "ekin (SI) " << ekin/pow(5.29e-11,3)*4.3597482e-18 << std::endl;
#endif

  //// find the inverse temperature (beta)
  int iter = 0;
  const int maxiter = 60;
  double tol = 1e-12;
  double K1 = 4.8966851;		// constants
  double K2 = 0.04496457;
  double K3 = 0.133376;
  double X = 4*ne*pow(M_PI/(2*me),1.5)*6.9608/6.95369; // FIXME conversion at end to match Sai's values...
#ifdef DEBUG_RTA
  std::cout << "XX " << X/pow(2.293710449e+17,1.5) << std::endl;
#endif
  double bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
  double bm = 0.0;
  double vol = pow(1.0/5.29e-11,3);		// volume element, multiply to go from a0^-3 to m^-3

  // loop applies Newton-Raphson method to get zero of function
  double f = 0.0;		// value of function (f)
  double fp = 0.0;		// value of function derivative (f')
#ifdef DEBUG_RTA
  std::cout << "Newton-Raphson to find inverse temperature" << std::endl;
#endif
  while ((fabs(bn - bm)/bm > tol) && (iter < maxiter)) {
    bm = bn;
    f = -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
    fp = -ekin + 2.25*ne*(K1/(K2*X*pow(bm,2.5))*log(1 + K2*X*pow(bm,1.5)) - K1/(bm*(1 + K2*X*pow(bm,1.5))) + 0.5*K3*X*pow(bm,0.5));
#ifdef DEBUG_RTA
    std::cout << "Iteration     " << std::setw(15) << iter << std::endl;
    std::cout << "bm            " << std::setw(15) << bm/4.3597482e-18 << std::endl;
    std::cout << "f(bm) term 1: " << std::setw(15) << vol*-bm*ekin << std::endl;
    std::cout << "f(bm) term 2: " << std::setw(15) << vol*1.5*ne*(1 + K1) << std::endl;
    std::cout << "f(bm) term 3: " << std::setw(15) << vol*1.5*ne*(-1*K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5))) << std::endl;
    std::cout << "f(bm) term 4: " << std::setw(15) << vol*1.5*ne*(0.5*K3*X*pow(bm,1.5)) << std::endl;
    std::cout << "f(bm) (SI)    " << std::setw(15) << f*pow(1.0/5.29e-11,3) << std::endl;
    std::cout << "f(bm) (a.u)   " << std::setw(15) << f << std::endl;
    std::cout << "f'(bm) (SI)   " << std::setw(15) << fp*pow(1.0/5.29e-11,3)*4.3597482e-18 << std::endl;
    std::cout << "f'(bm) (a.u)  " << std::setw(15) << fp << std::endl;
#endif
    bn = bm - f/fp;
    iter++;
  }
#ifdef DEBUG_RTA
  std::cout << std::endl;
#endif

  //// use beta to find chemical potential
  double mue = 0.0;
  double nue = 4*ne*pow(M_PI*bn/(2*me),1.5);	// constant to simplify
  mue = (log(nue) + K1*log(K2*nue + 1) + K3*nue)/bn;
#ifdef DEBUG_RTA
  std::cout << "Chemical potential " << mue*4.3597482e-18 << std::endl;
#endif

  // TODO account for temperature dropping in time
#ifdef DEBUG_RTA
  std::cout << "inverse temp is " << bn << std::endl;
  std::cout << std::endl;
#endif

  //// assign Fermi-Dirac function
  for (int ii = 0; ii < Ni; ii++) {
    // TODO factor in Boltzmann constant?
    fdd[ii] = 1.0/(1.0 + exp((E[ii] - mue)*bn));
#ifdef DEBUG_RTA
    std::cout << "FDD[" << ii << "]: " << std::scientific << fdd[ii] << std::endl;
#endif
  }

  //// normalize FDD to amount of population in conduction band
  double fddSum = 0.0;
  double CBSum = 0.0;
  for (int ii = 0; ii < Ni; ii++) {
    // sum population in FDD
    fddSum += fdd[ii];
    // sum population in states
    CBSum += y[(Ii + ii)*N + Ii + ii];
  }
  double fddNorm = CBSum/fddSum;
#ifdef DEBUG_RTA
  std::cout << "FDD normalization constant is " << fddNorm << std::endl;
#endif
  for (int ii = 0; ii < Ni; ii++) {
    fdd[ii] *= fddNorm;
  }
  if (std::isnan(mue)) {
#ifdef DEBUG_RTA
    std::cout << "mue is NaN!!!!!!!!!!!!!!" << std::endl;
#endif
    for (int ii = 0; ii < Ni; ii++) {
      fdd[ii] = y[(Ii + ii)*N + Ii + ii];
    }
    return;
  }
  else if (std::isnan(bn)) {
#ifdef DEBUG_RTA
    std::cout << "bn is NaN!!!!!!!!!!!!!!" << std::endl;
#endif
    for (int ii = 0; ii < Ni; ii++) {
      fdd[ii] = y[(Ii + ii)*N + Ii + ii];
    }
    return;
  }
  else {
#ifdef DEBUG_RTA
    std::cout << "no NaN++++++++++++" << std::endl;
#endif
    return;
  }
}

/* implements equation B13 from Binder et. al, PRB 1991.
 * bm is the beta (1/kT) value.
 * ekin is the kinetic energy
 * ne is the carrier density
 * K1-3 are constants
 * X is a constant
 */
double b13(double bm, double ekin, double ne, double K1, double K2, double K3, double X) {
  return -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
}


/* Right-hand-side equation for density matrix
 * using relaxation time approximation (RTA) */
int RHS_DM_RTA(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif

#ifdef DEBUG_RHS
  std::cout << "Time " << t << std::endl;
#endif

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  //std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  realtype * H = &(p->H)[0];
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g1 = p->gamma1;
  realtype g2 = p->gamma2;
  realtype g1_c = p->gamma1_c;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  // get bulk CB equilibrium FDD populations
  std::vector<double> fdd(p->Nk);
#ifdef DEBUG_RTA
  std::cout << "POPULATION " << yp[0] << std::endl;
#endif
  FDD_RTA(p, N_VGetArrayPointer(y), &(fdd[0]), 1);
  // force bulk conduction band toward Fermi-Dirac distribution
  for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
    ydotp[ii*N + ii] -= g1*(yp[ii*N + ii] - fdd[ii-p->Ik]);
  }

  if (p->rtaQD && (p->Nc > 1)) {
    //// relaxation in QD band
    fdd.resize(p->Nc);
    FDD_RTA(p, N_VGetArrayPointer(y), &(fdd[0]), 0);
    for (int ii = p->Ic; ii < (p->Ic + p->Nc); ii++) {
      ydotp[ii*N + ii] -= g1_c*(yp[ii*N + ii] - fdd[ii-p->Ic]);
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
	ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
	ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];

      // relaxation
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");
#endif

  return 0;
}

/* Right-hand-side equation for density matrix
 * using relaxation time approximation (RTA) */
int RHS_DM_RTA_BLAS(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif

#ifdef DEBUG_RHS
  std::cout << "Time " << t << std::endl;
#endif

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  //std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  realtype * H = &(p->H)[0];
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g1 = p->gamma1;
  realtype g2 = p->gamma2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }
  char TRANSA = 'n';
  char TRANSB = 'n';
  char LEFT = 'l';
  char RGHT = 'r';
  char UPLO = 'l';
  double ONE = 1.0;
  double NEG = -1.0;

  realtype * H_sp = &(p->H_sp)[0];
  int * columns = &(p->H_cols)[0];
  int * rowind = &(p->H_rowind)[0];

  // set up MKL variables
  double beta = 0.0;
  char matdescra [6] = {'T', // symmetric matrix
			'L', // lower triangle
			'N', // non-unit on diagonal
			'C', // zero-based indexing (C-style)
			'*', '*'};

			// set beta to zero for first call

  mkl_set_num_threads(p->nproc);
  // Re(\dot{\rho}) += H*Im(\rho)
  mkl_dcsrmm(&TRANSA, &N, &N, &N, &ONE, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[N2], &N, &beta, &ydotp[0], &N);
  
  // Re(\dot{\rho}) -= Im(\rho)*H
  DSYMM(&RGHT, &UPLO, &N, &N, &NEG, &H[0], &N, &yp[N2], &N, &ONE, &ydotp[0], &N);

  // Im(\dot{\rho}) += i*Re(\rho)*H
  DSYMM(&RGHT, &UPLO, &N, &N, &ONE, &H[0], &N, &yp[0], &N, &ONE, &ydotp[N2], &N);

  // Im(\dot{\rho}) -= i*H*Re(\rho)
  mkl_dcsrmm(&TRANSA, &N, &N, &N, &NEG, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[0], &N, &beta, &ydotp[N2], &N);

  //// diagonal; no need to calculate the imaginary part
  //   get equilibrium FDD populations
  //std::vector<double> fdd(p->Nk);
  double * fdd = new double [p->Nk];
#ifdef DEBUG_RTA
  std::cout << "POPULATION " << yp[0] << std::endl;
#endif
  FDD_RTA(p, N_VGetArrayPointer(y), fdd, 1);

  //// normalize FDD to amount of population in conduction band
  double fddSum = 0.0;
  double CBSum = 0.0;
  for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
    // sum population in FDD
    fddSum += fdd[ii - p->Ik];
    // sum population in CB
    CBSum += yp[ii*N + ii];
  }
  double fddNorm = CBSum/fddSum;
#ifdef DEBUG_RTA
  std::cout << "FDD normalization constant is " << fddNorm << std::endl;
#endif
  for (int ii = 0; ii < p->Nk; ii++) {
    fdd[ii] *= fddNorm;
  }

  // force conduction band toward Fermi-Dirac distribution
#pragma omp parallel for
  for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
    ydotp[ii*N + ii] -= g1*(yp[ii*N + ii] - fdd[ii]);
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      // relaxation
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");
#endif

  // free fdd
  delete [] fdd;

  return 0;
}

/* Right-hand-side equation for density matrix
 * using dephasing */
int RHS_DM_dephasing(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif


  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g2 = p->gamma2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    updateHamiltonian(p, t);
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
	ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
	ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];

      // relaxation
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");
#endif

  return 0;
}

