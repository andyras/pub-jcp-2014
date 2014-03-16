#include "output.hpp"

//#define DEBUG_OUTPUT
//#define DEBUG_OUTPUTTXPROB

/* returns true if map contains key, otherwise false */
bool isOutput(std::map<const std::string, bool> &myMap, const std::string myStr) {
  try {
    // key exists, return its value
    return myMap.at(myStr);
  }
  catch(const std::out_of_range& oor) {
    // key does not exist
    return false;
  }
}

/* returns output file name as char * */
std::string outputFileName(char * fileName, PARAMETERS * p) {
  // start with output directory name
  std::string fullFileName (p->outputDir);
  // add trailing slash if it is not there
  if (strcmp(&(fullFileName.at(fullFileName.length() - 1)), "/")) {
    fullFileName += "/";
  }
  fullFileName += fileName;

  return fullFileName;
}

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputFFTWVector(const char * fileName, fftw_complex * vec, double * x, int len) {
  // make output file
  FILE * output = NULL;
  output = fopen(fileName, "w");

  // write to the output
  for (int i = 0; i < len; i++) {
    fprintf(output, "%-.7g %-.7g %-.7g\n", x[i], vec[i][0], vec[i][1]);
  }

  // clean up the mess
  fclose(output);

  return;
}

/* Wrapper to outputFFTWVector, which fftshifts the output.
*/
void outputFFTWVectorShift(const char * fileName, fftw_complex * vec, double * x, int len) {
  // make a shifted copy of the vector to be printed
  fftw_complex * vec_shift = NULL;
  vec_shift = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*len);
  fftshift(vec, vec_shift, len);

  /*
  // make a shifted copy of the x array
  double * x_shift = new double [len];
  fftshift_double(x, x_shift, len);
  */

  // make the output
  outputFFTWVector(fileName, vec_shift, x, len);

  // clean up the mess
  fftw_free(vec_shift);
  /*
     delete [] x_shift;
     */

  return;
}

/* prints out initial wave function.  Inputs are the wave function array and
 * the number of equations.
 */
void outputWavefunction(realtype * psi, int n) {
  FILE * psi_start = NULL;
  psi_start = fopen("psi_start.out", "w");
  for (int i = 0; i < n; i++) {
    fprintf(psi_start, "%-.9g %-.9g\n", psi[i], psi[i+n]);
  }
  fclose(psi_start);
}

/* prints a vector W of length N */
void outputVector(realtype * W, int N, char * fileName) {
  FILE * out = NULL;   // output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e\n", W[ii]);
  }

  fclose(out);
}

/* prints the elementwise square of a complex vector */
void outputPsiSquare(complex16 * v, realtype * evals,  int N, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e %-.9e\n", evals[ii], (pow(v[ii].re,2) + pow(v[ii].im,2)));
  }

  fclose(out); 
}

/* prints a complex vector v of length N */
void outputCVector(complex16 * v, int N, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e %-.9e\n", v[ii].re, v[ii].im);
  }

  fclose(out);
}

/* prints a complex vector v of length N with M time steps */
void outputCVectorTime(complex16 * v, int N, int M, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int jj = 0; jj < M; jj++) {
    for (int ii = 0; ii < N; ii++) {
      fprintf(out, "%-.9e %-.9e\n", v[jj*N+ii].re, v[jj*N+ii].im);
    }
    fprintf(out, "\n");
  }

  fclose(out);
}

/* prints a square matrix M of dimension N */
void outputSquareMatrix(realtype * M, int N, char * fileName) {
  FILE * out = NULL; // output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      if (jj == 0) {
	fprintf(out, "%-.9e", M[ii*N + jj]);
      }
      else {
	fprintf(out, " %-.9e", M[ii*N + jj]);
      }
    }
    fprintf(out, "\n");
  }

  fclose(out);
}

/* prints a square matrix stored as a 2D array */
void output2DSquareMatrix(realtype ** M, int N, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e", M[ii][0]);
    for (int jj = 1; jj < N; jj++) {
      fprintf(out, " %-.9e", M[ii][jj]);
    }
    fprintf(out, "\n");
  }

  fclose(out);

  return;
}

/* Output the population in each state over time.  This function takes
 * the indices 'start' and 'end', e.g. Ik and Ik+Nk
 */
void outputXProbsWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << ".\n";
#endif
  std::ofstream output(outputFileName(fileName, p).c_str());

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii];
    for (int jj = start; jj < end; jj++) {
      output << " " << std::setw(8) << std::scientific
	<< (pow(wfnt[ii*p->NEQ*2 + jj],2) + pow(wfnt[ii*p->NEQ*2 + p->NEQ + jj],2));
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Output the integrated population on a set of states. */
void outputIntegralDM(char * fileName, const int start, const int end,
    const realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << std::endl;
#endif
  std::ofstream output(outputFileName(fileName, p).c_str());

  double summ1 = 0.0;	// accumulator variable, population at current time point
  double summ2 = 0.0;	// accumulator variable, population at previous time point
  double total = 0.0;	// actual integral over time
  int N = p->NEQ;
  int N2 = p->NEQ2;

  // get population at zeroth time point
  for (int jj = start; jj < end; jj++) {
    summ2 += dmt[jj*N + jj];
  }

  // output first time point
  output << std::setw(8) << std::scientific << p->times[0] << " " << 0 << std::endl;

  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << " ";
    summ1 = 0.0;
    // sum of population
    for (int jj = start; jj < end; jj++) {
      summ1 += dmt[ii*2*N2 + jj*N + jj];
    }
    // Riemann sum
    total += (p->times[ii+1] - p->times[ii])*(summ1 + summ2)/2.0;
    summ2 = summ1;
    output << std::setw(8) << std::scientific << total << std::endl;
  }

  output.close();

#ifdef DEBUG_OUTPUT
  std::cout << "Done making file " << fileName << std::endl;
#endif

  return;
}

/* Output the integrated population on a set of states. */
void outputIntegratedDM(char * fileName, const int start, const int end,
    const realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << std::endl;
#endif

  double summ1 = 0.0;	// accumulator variable, population at current time point
  double summ2 = 0.0;	// accumulator variable, population at previous time point
  double total = 0.0;	// actual integral over time
  int N = p->NEQ;
  int N2 = p->NEQ2;

  // get population at zeroth time point
  for (int jj = start; jj < end; jj++) {
    summ2 += dmt[jj*N + jj];
  }

  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // sum of population
    summ1 = 0.0;
    for (int jj = start; jj < end; jj++) {
      summ1 += dmt[ii*2*N2 + jj*N + jj];
    }
    // Riemann sum
    total += (p->times[ii+1] - p->times[ii])*(summ1 + summ2)/2.0;
    summ2 = summ1;
  }

  // write to output
  std::ofstream output(outputFileName(fileName, p).c_str());
  output << std::setw(8) << std::scientific << total << std::endl;
  output.close();

#ifdef DEBUG_OUTPUT
  std::cout << "Done making file " << fileName << std::endl;
#endif

  return;
}

/* Output the integrated population on a set of states. */
void outputIntegratedWfn(char * fileName, const int start, const int end,
    const realtype * wfnt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << std::endl;
#endif

  double summ1 = 0.0;	// accumulator variable, population at current time point
  double summ2 = 0.0;	// accumulator variable, population at previous time point
  double total = 0.0;	// actual integral over time
  int N = p->NEQ;

  // get population at zeroth time point
  for (int jj = start; jj < end; jj++) {
    summ2 += pow(wfnt[jj], 2) + pow(wfnt[N + jj], 2);
  }

  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // sum of population
    summ1 = 0.0;
    for (int jj = start; jj < end; jj++) {
      summ1 += pow(wfnt[ii*2*N + jj], 2) + pow(wfnt[ii*2*N + N + jj], 2);
    }
    // Riemann sum
    total += (p->times[ii+1] - p->times[ii])*(summ1 + summ2)/2.0;
    summ2 = summ1;
  }

  // write to output
  std::ofstream output(outputFileName(fileName, p).c_str());
  output << std::setw(8) << std::scientific << total << std::endl;
  output.close();

#ifdef DEBUG_OUTPUT
  std::cout << "Done making file " << fileName << std::endl;
#endif

  return;
}

/* Output the integrated population on a set of states. */
void outputIntegralWfn(char * fileName, const int start, const int end,
    const realtype * wfnt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << std::endl;
#endif
  std::ofstream output(outputFileName(fileName, p).c_str());

  double summ1 = 0.0;	// accumulator variable, population at current time point
  double summ2 = 0.0;	// accumulator variable, population at previous time point
  double total = 0.0;	// actual integral over time
  int N = p->NEQ;

  // get population at zeroth time point
  for (int jj = start; jj < end; jj++) {
    summ2 += pow(wfnt[jj], 2) + pow(wfnt[N + jj], 2);
  }

  // output first time point
  output << std::setw(8) << std::scientific << p->times[0] << " " << 0 << std::endl;

  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << " ";
    summ1 = 0.0;
    // sum of population
    for (int jj = start; jj < end; jj++) {
      summ1 += pow(wfnt[ii*2*N + jj], 2) + pow(wfnt[ii*2*N + N + jj], 2);
    }
    // Riemann sum
    total += (p->times[ii+1] - p->times[ii])*(summ1 + summ2)/2.0;
    summ2 = summ1;
    output << std::setw(8) << std::scientific << total << std::endl;
  }

  output.close();

#ifdef DEBUG_OUTPUT
  std::cout << "Done making file " << fileName << std::endl;
#endif

  return;
}

/* Output the population in each state over time.  This function takes
 * the indices 'start' and 'end', e.g. Ik and Ik+Nk
 */
void outputXProbs(char * fileName, int start, int end, realtype * dmt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << ".\n";
#endif
  std::ofstream output(outputFileName(fileName, p).c_str());

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii];
    for (int jj = start; jj < end; jj++) {
      output << " "
	<< std::setw(8) << std::scientific << dmt[ii*p->NEQ2*2 + jj*p->NEQ + jj];
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Output the total population in a set of states over time */
void outputtXprobWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUTTXPROB
  std::cout << "Making file " << fileName << "..." << std::endl;
  std::cout << "start index is " << start << std::endl;
  std::cout << "end index is " << end << std::endl;
#endif
  int N = p->NEQ;

  std::ofstream output(outputFileName(fileName, p).c_str());

  realtype summ;
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << " ";
    summ = 0.0;
    for (int jj = start; jj < end; jj++) {
      summ += pow(wfnt[ii*2*N + jj],2) + pow(wfnt[ii*2*N + jj + N],2);
    }
    output << std::setw(8) << std::scientific << summ << std::endl;
  }

  output.close();

#ifdef DEBUG_OUTPUTTXPROB
  std::cout << "Done making file " << fileName << "..." << std::endl;
#endif
  return;
}

/* Output the total population in a set of states over time.  This function
 * takes the indices 'start' and 'end', e.g. Ik and Ik+Nk
 */
void outputtXprob(char * fileName, int start, int end, realtype * dmt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << ".\n";
#endif
  std::ofstream output(outputFileName(fileName, p).c_str());
  realtype summ;

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << " ";
    summ = 0.0;
    for (int jj = start; jj < end; jj++) {
      summ += dmt[ii*p->NEQ2*2 + jj*p->NEQ + jj];
    }
    output << std::setw(8) << std::scientific << summ << "\n";
  }

  output.close();

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs the norm of each component of the density matrix */
void outputDMZ(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  fprintf(stderr, "\n\n\noutputting |\\rho| in time\n\n\n");
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  // loop over time steps
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // loop over first index
    for (int jj = 0; jj < p->NEQ; jj++) {
      // first element in row
      output << std::setw(8) << std::scientific
	<< sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj],2)
	    + pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + p->NEQ2],2));
      // loop over second index
      for (int kk = 1; kk < p->NEQ; kk++) {
	output << " " << std::setw(8) << std::scientific
	  << sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk],2)
	      + pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk + p->NEQ2],2));
      }
      output << "\n";
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs the norm of each component of the density matrix */
void outputDMCoherences(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  fprintf(stderr, "\n\n\noutputting |\\rho_{ij}| in time\n\n\n");
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  // loop over time steps
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // loop over first index
    for (int jj = 0; jj < p->NEQ; jj++) {
      // first element in row
      if (jj == 0) {
	output << std::setw(8) << std::scientific << 0;
      }
      else {
	output << std::setw(8) << std::scientific
	  << sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj],2)
	      + pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + p->NEQ2],2));
      }
      // loop over second index
      for (int kk = 1; kk < p->NEQ; kk++) {
	if (jj == kk) {
	  output << " " << std::setw(8) << std::scientific << 0;
	}
	else {
	  output << " " << std::setw(8) << std::scientific
	    << sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk],2)
		+ pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk + p->NEQ2],2));
	}
      }
      output << "\n";
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs the real part of each component of the density matrix */
void outputDMRe(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  fprintf(stderr, "\n\n\noutputting Re(\\rho) in time\n\n\n");
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  // loop over time steps
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // loop over first index
    for (int jj = 0; jj < p->NEQ; jj++) {
      // first element in row
      output << std::setw(8) << std::scientific
	<< dmt[2*p->NEQ2*ii + p->NEQ*jj];
      // loop over second index
      for (int kk = 1; kk < p->NEQ; kk++) {
	output << " " << std::setw(8) << std::scientific
	  << dmt[2*p->NEQ2*ii + p->NEQ*jj + kk];
      }
      output << "\n";
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs the imaginary part of each component of the density matrix */
void outputDMIm(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  fprintf(stderr, "\n\n\noutputting Im(\\rho) in time\n\n\n");
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  // loop over time steps
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // loop over first index
    for (int jj = 0; jj < p->NEQ; jj++) {
      // first element in row
      output << std::setw(8) << std::scientific
	<< dmt[2*p->NEQ2*ii + p->NEQ*jj + p->NEQ2];
      // loop over second index
      for (int kk = 1; kk < p->NEQ; kk++) {
	output << " " << std::setw(8) << std::scientific
	  << dmt[2*p->NEQ2*ii + p->NEQ*jj + kk + p->NEQ2];
      }
      output << "\n";
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs energies of all states */
void outputEnergy(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "\nMaking " << fileName << "\n";
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  for (int ii = 0; ii < p->NEQ; ii++) {
    output << std::setw(8) << std::scientific << p->energies[ii] << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs energy expectation value over time */
void outputEnergyExpWfn(char * fileName, struct PARAMETERS * p,
    double * wfnt) {
  // allocate array to hold matrix-vector product
  std::vector<double> psiHvec (2*p->NEQ, 0.0);
  double * psiH = &(psiHvec[0]);
  // array of energy expectation values
  std::vector<double> ee (p->numOutputSteps, 0.0);

  double * H = &(p->H[0]);
  int N = p->NEQ;

  char TRANS = 'n';
  double ONE = 1.0;
  double ZERO = 0.0;
  int INC = 1;

  std::ofstream output(outputFileName(fileName, p).c_str());

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    // update Hamiltonian
    if (p->laser_on || p->torsion) {
      updateHamiltonian(p, ii*p->tout/p->numOutputSteps);
    }
    // take product \hat{H}(t)\ket{\psi(t)} (real part)
    DGEMV(&TRANS, &N, &N, &ONE, &H[0], &N, &wfnt[ii*N*2], &INC, &ZERO, &psiH[0], &INC);
    // take product \hat{H}(t)\ket{\psi(t)} (imag part)
    DGEMV(&TRANS, &N, &N, &ONE, &H[0], &N, &wfnt[ii*N*2 + N], &INC, &ZERO, &psiH[N], &INC);
    // take dot product of \bra{\psi(t)} with Hpsi(t) (real part)
    ee[ii] = DDOT(&N, &wfnt[ii*N*2], &INC, &psiH[0], &INC);
    // take dot product of \bra{\psi(t)} with Hpsi(t) (imag part)
    ee[ii] += DDOT(&N, &wfnt[ii*N*2 + N], &INC, &psiH[N], &INC);
    output << p->times[ii] << " " << ee[ii] << std::endl;
  }

  output.close();

  return;
}

/* Outputs all the time steps */
void outputTimes(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "\nMaking " << fileName << "\n";
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs expectation value of energy at all times */
void outputEnergyExp(char * fileName, realtype * dmt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "\nMaking " << fileName << "\n";
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  realtype summ ;

  // loop over timesteps
  for (int kk = 0; kk <= p->numOutputSteps; kk++) {
    summ = 0.0;
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += p->H[ii*p->NEQ + ii]*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii];
      for (int jj = 0; jj < ii; jj++) {
	summ += 2*p->H[ii*p->NEQ + jj]*dmt[kk*p->NEQ2*2 + jj*p->NEQ + ii];
      }
    }
    output << std::setw(8) << std::scientific
      << p->times[kk] << " "
      << std::setw(8) << std::scientific
      << summ << std::endl;
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Output Fermi level calculated from total band populations */
void outputDynamicMu(char * fileName, realtype * dmt, int bandFlag, struct PARAMETERS * p) {
  int start = bandStartIdx(bandFlag, p);
  int end = bandEndIdx(bandFlag, p);
  int N = p->NEQ;
  int N2 = p->NEQ2;
  double summ = 0.0;
  double T = p->temperature;

  std::ofstream output(outputFileName(fileName, p).c_str());

  for (int kk = 0; kk < p->numOutputSteps; kk++) {
    summ = 0.0;
    for (int ii = start; ii < end; ii++) {
      summ += dmt[kk*N2*2 + ii*N + ii];
    }
    output << std::setw(8) << std::scientific
      << p->times[kk] << " "
      << std::setw(8) << std::scientific
      << findDynamicMu(summ, T, bandFlag, p) << std::endl;
  }

  return;
}

/* Outputs Fermi level as calculated from populations */
void outputMuFromPops(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "\nMaking " << fileName << "\n";
#endif

  std::ofstream output(outputFileName(fileName, p).c_str());
  std::ofstream average("muFromPopsAvg.out");

  double kT = p->temperature/3.185e5;
  double pop, mu;
  double summ = 0.0;
  int N = p->NEQ;
  int N2 = p->NEQ2;

  // loop over timesteps
  for (int kk = 0; kk <= p->numOutputSteps; kk++) {
    // print time
    output << std::setw(8) << std::scientific << p->times[kk];

    summ = 0.0;
    // loop over k states
    for (int ii = 0; ii < p->Nk; ii++) {
      // inverse of Fermi-Dirac function
      pop = dmt[kk*N2*2 + ii*N + ii];
      // FIXME this is broken when pop > 1.0
      mu = p->energies[p->Ik + ii] - kT*log((1.0 - pop)/pop);
      output << " " << std::setw(8) << std::scientific << mu;
      summ += mu;
    }

    output << std::endl;

    average << std::setw(8) << std::scientific << p->times[kk] << " "
      << std::setw(8) << std::scientific << summ/p->Nk << std::endl;
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs torsional potential at simulation time points */
void outputTorsion(struct PARAMETERS * p, char * fileName) {
  std::ofstream output(outputFileName(fileName, p).c_str());

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific
      << p->times[ii] << " "
      << std::setw(8) << std::scientific
      << p->torsionV->value(p->times[ii]) << std::endl;
  }
  return;
}

/* Outputs quantities from RTA */
void outputRTA(std::map<const std::string, bool> &outs,
    realtype * dmt, struct PARAMETERS * p) {

  double ne = 0.0;
  double ekin = 0.0;
  double factor = 1.0/(M_PI*M_PI*pow(p->X2,3));
  // variables for Simpson's rule integration
  double SF = 4.0;	// Simpson's factor
  int sign = -1;	// sign
  // variables for calculating mu_e and beta
  int iter = 0;
  const int maxiter = 60;
  const double tol = 1e-12;
  const double K1 = 4.8966851;		// constants
  const double K2 = 0.04496457;
  const double K3 = 0.133376;
  const double X = 4*ne*pow(M_PI/(2*p->me),1.5)*6.9608/6.95369; // FIXME conversion at end to match Sai's values...
  double bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
  double bm = 0.0;
  const double vol = pow(1.0/5.29e-11,3);		// volume element, multiply to go from a0^-3 to m^-3
  double f = 0.0;		// value of function (f)
  double fp = 0.0;		// value of function derivative (f')
  double mue = 0.0;
  double nue = 4*ne*pow(M_PI*bm/(2*p->me),1.5);	// constant to simplify

  // vectors for time-dependent quantities
  std::vector<double> mu_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> temp_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> ne_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> ekin_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> fdd_t ((p->numOutputSteps + 1)*p->Nk, 0.0);

  // loop over timesteps
  for (int kk = 0; kk <= p->numOutputSteps; kk++) {
    //// calculate n_e and e_kin

    // assign vector of energies
    std::vector<double> E (p->Nk,0.0);
    for (int ii = 0; ii < p->Nk; ii++) {
      E[ii] = pow(ii,2)/(2*p->me*pow(p->X2,2));
    }

    // Simpson's Rule method
    SF = 4.0;
    sign = -1;
    // skip the first point because the value will be zero
    for (int ii = 1; ii < (p->Nk-1); ii++) {
      ne += SF*factor*ii*ii*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii];
      ekin += SF*factor*pow(ii,2)*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii]*E[ii];
      SF += sign*2.0;
      sign *= -1;
    }
    // add last point
    ne += factor*pow(p->Nk-1,2)*dmt[kk*p->NEQ2*2 + (p->Nk - 1)*p->NEQ + p->Nk - 1];
    ekin += factor*pow(p->Nk-1,2)*dmt[kk*p->NEQ2*2 + (p->Nk - 1)*p->NEQ + p->Nk - 1]*E[p->Nk - 1];
    // divide by three
    ne /= 3.0;
    ekin /= 3.0;

    ne_t[kk] = ne;
    ekin_t[kk] = ekin;

    //// find the inverse temperature (beta)
    bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
    bm = 0.0;

    // loop applies Newton-Raphson method to get zero of function
    f = 0.0;		// value of function (f)
    fp = 0.0;		// value of function derivative (f')
    while ((fabs(bn - bm)/bm > tol) && (iter < maxiter)) {
      bm = bn;
      f = -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
      fp = -ekin + 2.25*ne*(K1/(K2*X*pow(bm,2.5))*log(1 + K2*X*pow(bm,1.5)) - K1/(bm*(1 + K2*X*pow(bm,1.5))) + 0.5*K3*X*pow(bm,0.5));
      bn = bm - f/fp;
      iter++;
    }

    temp_t[kk] = 1.0/bn;
    //// use beta to find chemical potential
    nue = 4*ne*pow(M_PI*bn/(2*p->me),1.5);	// constant to simplify
    mue = (log(nue) + K1*log(K2*nue + 1) + K3*nue)/bn;
#ifdef DEBUG_OUTPUT
    std::cout << "nue is " << nue << std::endl;
    std::cout << "(log(" << nue << ") + " << K1 << "*log(" << K2 << "*" << nue << " + 1) + " << K3 << "*" << nue << ")/" << bn << std::endl;
    std::cout << "mue is " << mue << std::endl;
#endif
    mu_t[kk] = mue;

    // update FDD
    for (int ii = 0; ii < p->Nk; ii++) {
      fdd_t[kk*p->Nk + ii] = 1.0/(1.0 + exp((E[ii] - mue)*bn));
    }
  }
  if (isOutput(outs, "mu.out")) {
    std::ofstream mu_out("mu.out");
    for (int kk = 0; kk <= p->numOutputSteps; kk++) {
      mu_out << p->times[kk] << " " << mu_t[kk] << std::endl;
    }
    mu_out.close();

#ifdef DEBUG_OUTPUT
    std::cout << "\nDone making mu.out" << std::endl;
#endif
  }

  if (isOutput(outs, "temp.out")) {
    std::ofstream temp_out("temp.out");
    for (int kk = 0; kk <= p->numOutputSteps; kk++) {
      temp_out << p->times[kk] << " " << temp_t[kk]*AU2K << std::endl;
    }
    temp_out.close();

#ifdef DEBUG_OUTPUT
    std::cout << "\nDone making temp.out" << std::endl;
#endif
  }

  if (isOutput(outs, "ne.out")) {
    std::ofstream ne_out("ne.out");
    for (int kk = 0; kk <= p->numOutputSteps; kk++) {
      ne_out << p->times[kk] << " " << ne_t[kk] << std::endl;
    }
    ne_out.close();

#ifdef DEBUG_OUTPUT
    std::cout << "\nDone making ne.out" << std::endl;
#endif
  }

  if (isOutput(outs, "ekin.out")) {
    std::ofstream ekin_out("ekin.out");
    for (int kk = 0; kk <= p->numOutputSteps; kk++) {
      ekin_out << p->times[kk] << " " << ekin_t[kk] << std::endl;
    }
    ekin_out.close();

#ifdef DEBUG_OUTPUT
    std::cout << "\nDone making ekin.out" << std::endl;
#endif
  }

  if (isOutput(outs, "fdd.out")) {
    std::ofstream fdd_out("fdd.out");
    for (int kk = 0; kk <= p->numOutputSteps; kk++) {
      fdd_out << p->times[kk];
      for (int ii = 0; ii < p->Nk; ii++) {
	fdd_out << " " << fdd_t[kk*p->Nk + ii];
      }
      fdd_out << std::endl;
    }
    fdd_out.close();

#ifdef DEBUG_OUTPUT
    std::cout << "\nDone making fdd.out" << std::endl;
#endif
  }

  return;
}

/* Outputs the laser pump intensity over time */
void outputPumpIntensity(struct PARAMETERS * p, char * fileName) {
  std::ofstream o(fileName);
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    o << p->times[ii] << " "
      << gaussPulse(p->times[ii], p->pumpFWHM, p->pumpAmpl, p->pumpPeak, p->pumpFreq, p->pumpPhase)
      << std::endl;
  }
  o.close();

  return;
}

/* output couplings as a matrix */
void outputCouplings(struct PARAMETERS * p, char * fileName) {
  std::ofstream o(fileName);
  //// This output is the same as the Hamiltonian, but with diagonal elements zero.

  // first element of first row will be zero
  o << 0;
  // rest of first row
  for (int jj = 1; jj < p->NEQ; jj++) {
    o << " " << p->H[jj];
  }
  o << std::endl;

  // remaining rows
  for (int ii = 1; ii < p->NEQ; ii++) {
    // first element in row
    o << p->H[ii*p->NEQ];
    // rest of row
    for (int jj = 1; jj < p->NEQ; jj++) {
      if (ii == jj) {
	o << " " << 0;
      }
      else {
	o << " " << p->H[ii*p->NEQ + jj];
      }
    }
    o << std::endl;
  }
  o.close();

  return;
}

/* Finds peaks in populations, outputs values, times and differences. */
void findPeaksWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Making file " << fileName << "..." << std::endl;
  std::cout << "start index is " << start << std::endl;
  std::cout << "end index is " << end << std::endl;
#endif
  int N = p->NEQ;
  // vector of populations
  std::vector<double> pops (p->numOutputSteps, 0.0);

  // struct to hold peak info
  struct Peak {
    double peak;
    double time;
    double nextPeakTime;
  };
  std::vector<Peak> peaks;

  // calculate populations
  realtype summ;
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    summ = 0.0;
    for (int jj = start; jj < end; jj++) {
      summ += pow(wfnt[ii*2*N + jj],2) + pow(wfnt[ii*2*N + jj + N],2);
    }
    pops[ii] = summ;
  }

  //// find peaks
  Peak tmpPeak;	// struct to hold peak information before adding to vector
  tmpPeak.nextPeakTime = 0.0;

  // first point can be a peak
  if (pops[0] > pops[1]) {
    tmpPeak.peak = pops[0];
    tmpPeak.time = p->times[0];
    peaks.push_back(tmpPeak);
  }

  for (int ii = 1; ii < p->numOutputSteps; ii++) {
    if ((pops[ii] > pops[ii-1]) && (pops[ii] > pops[ii+1])) {
      tmpPeak.peak = pops[ii];
      tmpPeak.time = p->times[ii];
      peaks.push_back(tmpPeak);
    }
  }

  // last point can be a peak
  if (pops[p->numOutputSteps] > pops[p->numOutputSteps-1]) {
    tmpPeak.peak = pops[p->numOutputSteps];
    tmpPeak.time = p->times[p->numOutputSteps];
    peaks.push_back(tmpPeak);
  }

  //// find time gaps between peaks, last one is zero
#ifdef DEBUG_PEAKS
  std::cout << "size of peaks vector is " << peaks.size() << std::endl;
#endif
  if (peaks.size() > 1) {	// don't do this with zero or one peak
    for (int ii = 0; ii < (peaks.size()-1); ii++) {
      peaks[ii].nextPeakTime = peaks[ii+1].time - peaks[ii].time;
    }
  }

  // create output file
  std::ofstream output(outputFileName(fileName, p).c_str());

  // handle case with no peaks (flat the whole time)
  if (peaks.size() == 0) {
    output << "0 0 0" << std::endl;
#ifdef DEBUG_OUTPUT
    std::cout << "No peaks." << std::endl;
#endif
  }
  else {
    for (int ii = 0; ii < peaks.size(); ii++) {
#ifdef DEBUG_OUTPUT
      std::cout << "Peak (" << peaks[ii].peak << ") at " << peaks[ii].time << std::endl;
#endif
      output << std::setw(8) << std::scientific << peaks[ii].time << " "
	<< std::setw(8) << std::scientific << peaks[ii].peak << " "
	<< std::setw(8) << std::scientific << peaks[ii].nextPeakTime << std::endl;
    }
  }

  output.close();

#ifdef DEBUG_OUTPUT
  std::cout << "Done making file " << fileName << "..." << std::endl;
#endif
  return;
}

void outputDeriv(std::string outFile, int n, realtype * deriv, struct PARAMETERS * p) {
  std::ofstream o (outFile.c_str());
  int nt = p->numOutputSteps;

  // major index of array is state, so this loop outputs the transpose
  for (int ii = 0; ii < (nt-5); ii++) {
    o << std::setw(8) << std::scientific << p->times[ii+2];
    for (int jj = 0; jj < n; jj++) {
      o << " " << std::setw(8) << std::scientific << deriv[jj*nt + ii];
    }
    o << std::endl;
  }

  return;
}

/* yt can be either wavefunction or DM over time */
void outputDerivsWfn(std::map<const std::string, bool> &outs, realtype * yt,
    struct PARAMETERS * p){
  // unpack values from p
  int N = p->NEQ;
  int nt = p->numOutputSteps;

  // create array of populations
  std::vector<realtype> pops (N*nt, 0.0);

  // fill array of populations
  if (p->wavefunction) {
    for (int ii = 0; ii < N; ii++) {
      for (int jj = 0; jj < nt; jj++) {
	pops[ii*nt + jj] = pow(yt[jj*2*N + ii], 2) + pow(yt[jj*2*N + ii + N], 2);
      }
    }
  }
  else { // DM
    for (int ii = 0; ii < N; ii++) {
      for (int jj = 0; jj < nt; jj++) {
	pops[ii*nt + jj] = yt[jj*2*N*N + ii*N + ii];
      }
    }
  }

  // deriv of CB populations
  if (isOutput(outs, "derivKprobs.out")) {
    // create output array
    std::vector<realtype> kderivs (p->Nk*(p->numOutputSteps-5), 0.0);
    // take derivative
    arrayDeriv(&(pops[p->Ik*nt]), nt, p->Nk, N, &(kderivs[0]), p->tout/(nt-1));
    // print derivative
    outputDeriv("derivKprobs.out", p->Nk, &(kderivs[0]), p);
  }

  // deriv of QD populations
  if (isOutput(outs, "derivCprobs.out")) {
    std::vector<realtype> cderivs (p->Nc*(nt-5), 0.0);
    arrayDeriv(&(pops[p->Ic*nt]), nt, p->Nc, N, &(cderivs[0]), p->tout/(nt-1));
    outputDeriv("derivCprobs.out", p->Nc, &(cderivs[0]), p);
  }

  // deriv of bridge populations
  if (isOutput(outs, "derivBprobs.out")) {
    std::vector<realtype> bderivs (p->Nb*(nt-5), 0.0);
    arrayDeriv(&(pops[p->Ib*nt]), nt, p->Nb, N, &(bderivs[0]), p->tout/(nt-1));
    outputDeriv("derivBprobs.out", p->Nb, &(bderivs[0]), p);
  }

  // deriv of VB populations
  if (isOutput(outs, "derivLprobs.out")) {
    std::vector<realtype> lderivs (p->Nl*(nt-5), 0.0);
    arrayDeriv(&(pops[p->Il*nt]), nt, p->Nl, N, &(lderivs[0]), p->tout/(nt-1));
    outputDeriv("derivLprobs.out", p->Nl, &(lderivs[0]), p);
  }

  // deriv of all populations
  if (isOutput(outs, "derivAllprobs.out")) {
    std::vector<realtype> derivs (N*(nt-5));
    arrayDeriv(&(pops[0]), nt, N, N, &(derivs[0]), p->tout/(nt-1));
    outputDeriv("derivAllprobs.out", N, &(derivs[0]), p);
  }

  // deriv of total CB population
  if (isOutput(outs, "derivTkprob.out")) {
    // create population array
    std::vector<realtype> tkpops (nt, 0.0);
    for (int ii = 0; ii < nt; ii++) {
      for (int jj = p->Ik; jj < (p->Ik+p->Nk); jj++) {
	tkpops[ii] += pops[jj*nt + ii];
      }
    }
    // create deriv array
    std::vector<realtype> tkderivs (nt-5);
    // take derivative
    arrayDeriv(&(tkpops[0]), nt, 1, 1, &(tkderivs[0]), p->tout/(nt-1));
    // print derivative
    outputDeriv("derivTkprob.out", 1, &(tkderivs[0]), p);
  }

  // deriv of total QD population
  if (isOutput(outs, "derivTcprob.out")) {
    std::vector<realtype> tcpops (nt, 0.0);
    for (int ii = 0; ii < nt; ii++) {
      for (int jj = p->Ic; jj < (p->Ic+p->Nc); jj++) {
	tcpops[ii] += pops[jj*nt + ii];
      }
    }
    std::vector<realtype> tcderivs (nt-5);
    arrayDeriv(&(tcpops[0]), nt, 1, 1, &(tcderivs[0]), p->tout/(nt-1));
    outputDeriv("derivTcprob.out", 1, &(tcderivs[0]), p);
  }

  // deriv of total bridge population
  if (isOutput(outs, "derivTbprob.out")) {
    std::vector<realtype> tbpops (nt, 0.0);
    for (int ii = 0; ii < nt; ii++) {
      for (int jj = p->Ib; jj < (p->Ib+p->Nb); jj++) {
	tbpops[ii] += pops[jj*nt + ii];
      }
    }
    std::vector<realtype> tbderivs (nt-5);
    arrayDeriv(&(tbpops[0]), nt, 1, 1, &(tbderivs[0]), p->tout/(nt-1));
    outputDeriv("derivTbprob.out", 1, &(tbderivs[0]), p);
  }

  // deriv of total VB population
  if (isOutput(outs, "derivTlprob.out")) {
    std::vector<realtype> tlpops (nt, 0.0);
    for (int ii = 0; ii < nt; ii++) {
      for (int jj = p->Il; jj < (p->Il+p->Nl); jj++) {
	tlpops[ii] += pops[jj*nt + ii];
      }
    }
    std::vector<realtype> tlderivs (nt-5);
    arrayDeriv(&(tlpops[0]), nt, 1, 1, &(tlderivs[0]), p->tout/(nt-1));
    outputDeriv("derivTlprob.out", 1, &(tlderivs[0]), p);
  }

  return;
}

void outputDerivsDM(std::map<const std::string, bool> &outs, realtype * dmt,
    struct PARAMETERS * p){
  return;
}

/* Computes outputs independent of DM or wavefunction propagation*/
void computeGeneralOutputs(std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {
  // torsion-mediated coupling
  if ((p->torsion) && (isOutput(outs, "torsion.out"))) {
    outputTorsion(p, "torsion.out");
  }

  // hamiltonian at time zero
  if (isOutput(outs, "ham.out")) {
    // &(p->H)[0] is address of first element of array in vector
    outputSquareMatrix(&(p->H)[0], p->NEQ, "ham.out");
  }

  // pump intensity
  if (isOutput(outs, "pump_intensity.out")) {
    outputPumpIntensity(p, "pump_intensity.out");
  }

  // energies of all states
  if (isOutput(outs, "energies.out")) {
    outputEnergy("energies.out", p);
  }

  if (isOutput(outs, "couplings.out")) {
    outputCouplings(p, "couplings.out");
  }

  return;
}

/* Computes outputs from \psi(t) */
void computeWfnOutput(realtype * wfnt, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {
  // total population on all sites
  if (isOutput(outs, "totprob.out")) {
    outputtXprobWfn("totprob.out", 0, p->NEQ, wfnt, p);
  }

  // total population on bulk conduction band
  if (isOutput(outs, "tkprob.out")) {
    outputtXprobWfn("tkprob.out", p->Ik, p->Ik + p->Nk, wfnt, p);
  }

  // total population on QD
  if (isOutput(outs, "tcprob.out")) {
    outputtXprobWfn("tcprob.out", p->Ic, p->Ic + p->Nc, wfnt, p);
  }

  // total population on bridge
  if (isOutput(outs, "tbprob.out")) {
    outputtXprobWfn("tbprob.out", p->Ib, p->Ib + p->Nb, wfnt, p);
  }

  // total population on bulk valence band
  if (isOutput(outs, "tlprob.out")) {
    outputtXprobWfn("tlprob.out", p->Il, p->Il + p->Nl, wfnt, p);
  }

  // populations in all states
  if (isOutput(outs, "allprobs.out")) {
    outputXProbsWfn("allprobs.out", 0, p->NEQ, wfnt, p);
  }

  // populations in bulk CB
  if (isOutput(outs, "kprobs.out")) {
    outputXProbsWfn("kprobs.out", p->Ik, p->Ik + p->Nk, wfnt, p);
  }

  // populations in QD
  if (isOutput(outs, "cprobs.out")) {
    outputXProbsWfn("cprobs.out", p->Ic, p->Ic + p->Nc, wfnt, p);
  }

  // populations in bridge states
  if (isOutput(outs, "bprobs.out")) {
    outputXProbsWfn("bprobs.out", p->Ib, p->Ib + p->Nb, wfnt, p);
  }

  // populations in bulk VB
  if (isOutput(outs, "lprobs.out")) {
    outputXProbsWfn("lprobs.out", p->Il, p->Il + p->Nl, wfnt, p);
  }

  // integrated population on QD
  if (isOutput(outs, "cumItkprob.out")) {
    outputIntegralWfn("cumItkprob.out", p->Ik, p->Ik + p->Nk, wfnt, p);
  }
  if (isOutput(outs, "cumItcprob.out")) {
    outputIntegralWfn("cumItcprob.out", p->Ic, p->Ic + p->Nc, wfnt, p);
  }
  if (isOutput(outs, "cumItbprob.out")) {
    outputIntegralWfn("cumItbprob.out", p->Ib, p->Ib + p->Nb, wfnt, p);
  }
  if (isOutput(outs, "cumItlprob.out")) {
    outputIntegralWfn("cumItlprob.out", p->Il, p->Il + p->Nl, wfnt, p);
  }

  // integrated population on QD over all time
  if (isOutput(outs, "Itkprob.out")) {
    outputIntegratedWfn("Itkprob.out", p->Ik, p->Ik + p->Nk, wfnt, p);
  }
  if (isOutput(outs, "Itcprob.out")) {
    outputIntegratedWfn("Itcprob.out", p->Ic, p->Ic + p->Nc, wfnt, p);
  }
  if (isOutput(outs, "Itbprob.out")) {
    outputIntegratedWfn("Itbprob.out", p->Ib, p->Ib + p->Nb, wfnt, p);
  }
  if (isOutput(outs, "Itlprob.out")) {
    outputIntegratedWfn("Itlprob.out", p->Il, p->Il + p->Nl, wfnt, p);
  }

  // energy expectation value
  if (isOutput(outs, "energyexp.out")) {
    outputEnergyExpWfn("energyexp.out", p, wfnt);
  }

  // peaks in populations
  if (isOutput(outs, "peaksTkprob.out")) {
    findPeaksWfn("peaksTkprob.out", p->Ik, p->Ik + p->Nk, wfnt, p);
  }
  if (isOutput(outs, "peaksTcprob.out")) {
    findPeaksWfn("peaksTcprob.out", p->Ic, p->Ic + p->Nc, wfnt, p);
  }
  if (isOutput(outs, "peaksTbprob.out")) {
    findPeaksWfn("peaksTbprob.out", p->Ib, p->Ib + p->Nb, wfnt, p);
  }
  if (isOutput(outs, "peaksTlprob.out")) {
    findPeaksWfn("peaksTlprob.out", p->Il, p->Il + p->Nl, wfnt, p);
  }

  // derivatives of populations
  outputDerivsWfn(outs, wfnt, p);

  return;
}

/* Computes outputs from \rho(t) */
void computeDMOutput(realtype * dmt, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {

  // total population
  if (isOutput(outs, "totprob.out")) {
    outputtXprob("totprob.out", 0, p->NEQ, dmt, p);
  }

  // populations in k states
  if (isOutput(outs, "kprobs.out")) {
    outputXProbs("kprobs.out", p->Ik, p->Ik + p->Nk, dmt, p);
  }
  if (isOutput(outs, "tkprob.out")) {
    outputtXprob("tkprob.out", p->Ik, p->Ik + p->Nk, dmt, p);
  }

  // populations in c states
  if (isOutput(outs, "cprobs.out")) {
    outputXProbs("cprobs.out", p->Ic, p->Ic + p->Nc, dmt, p);
  }
  if (isOutput(outs, "tcprob.out")) {
    outputtXprob("tcprob.out", p->Ic, p->Ic + p->Nc, dmt, p);
  }

  // populations in b states
  if (isOutput(outs, "bprobs.out")) {
    outputXProbs("bprobs.out", p->Ib, p->Ib + p->Nb, dmt, p);
  }
  if (isOutput(outs, "tbprob.out")) {
    outputtXprob("tbprob.out", p->Ib, p->Ib + p->Nb, dmt, p);
  }

  // populations in l states
  if (isOutput(outs, "lprobs.out")) {
    outputXProbs("lprobs.out", p->Il, p->Il + p->Nl, dmt, p);
  }
  if (isOutput(outs, "tlprob.out")) {
    outputtXprob("tlprob.out", p->Il, p->Il + p->Nl, dmt, p);
  }

  // integrated populations
  if (isOutput(outs, "cumItkprob.out")) {
    outputIntegralDM("cumItkprob.out", p->Ik, p->Ik + p->Nk, dmt, p);
  }
  if (isOutput(outs, "cumItcprob.out")) {
    outputIntegralDM("cumItcprob.out", p->Ic, p->Ic + p->Nc, dmt, p);
  }
  if (isOutput(outs, "cumItbprob.out")) {
    outputIntegralDM("cumItbprob.out", p->Ib, p->Ib + p->Nb, dmt, p);
  }
  if (isOutput(outs, "cumItlprob.out")) {
    outputIntegralDM("cumItlprob.out", p->Il, p->Il + p->Nl, dmt, p);
  }

  // integrated populations over all time
  if (isOutput(outs, "Itkprob.out")) {
    outputIntegratedDM("Itkprob.out", p->Ik, p->Ik + p->Nk, dmt, p);
  }
  if (isOutput(outs, "Itcprob.out")) {
    outputIntegratedDM("Itcprob.out", p->Ic, p->Ic + p->Nc, dmt, p);
  }
  if (isOutput(outs, "Itbprob.out")) {
    outputIntegratedDM("Itbprob.out", p->Ib, p->Ib + p->Nb, dmt, p);
  }
  if (isOutput(outs, "Itlprob.out")) {
    outputIntegratedDM("Itlprob.out", p->Il, p->Il + p->Nl, dmt, p);
  }

  // norm of DM elements
  if (isOutput(outs, "dmt_z.out")) {
    outputDMZ("dmt_z.out", dmt, p);
  }

  // norm of DM elements
  if (isOutput(outs, "dmt_re.out")) {
    outputDMRe("dmt_re.out", dmt, p);
  }

  // norm of DM elements
  if (isOutput(outs, "dmt_im.out")) {
    outputDMIm("dmt_im.out", dmt, p);
  }

  // coherences (magnitude)
  if (isOutput(outs, "dmCoherences.out")) {
    outputDMCoherences("dmCoherences.out", dmt, p);
  }

  // all time steps
  if (isOutput(outs, "times.out")) {
    outputTimes("times.out", p);
  }

  // expectation value of energy
  if (isOutput(outs, "energyexp.out")) {
    outputEnergyExp("energyexp.out", dmt, p);
  }

  // Fermi level as calculated from populations
  if (isOutput(outs, "muFromPops.out")) {
    outputMuFromPops("muFromPops.out", dmt, p);
  }

  // Fermi level in bulk as calculated from populations
  if (isOutput(outs, "dynamicMuBulk.out")) {
    outputDynamicMu("dynamicMuBulk.out", dmt, CONDUCTION, p);
  }

  // Fermi level in bulk as calculated from populations
  if (isOutput(outs, "dynamicMuQD.out")) {
    outputDynamicMu("dynamicMuQD.out", dmt, QD_CONDUCTION, p);
  }

  // derivatives of populations
  outputDerivsDM(outs, dmt, p);

  // RTA outputs are tied together
  outputRTA(outs, dmt, p);

  return;
}
