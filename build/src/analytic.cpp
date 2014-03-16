#include "analytic.hpp"

//#define DEBUG_ANALYTIC

/* Compute analytic dynamics */
void computeAnalyticOutputs(std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {

  // energy spacing in bulk
  std::complex <double> dE ((p->kBandTop-p->kBandEdge)/(p->Nk-1), 0);
  // bulk-QD coupling
  std::complex <double> Vee (p->Vnobridge[0], 0);
  // rate constant (can be defined also as K/2)
  std::complex <double> K = std::complex <double> (3.1415926535,0)*pow(Vee,2)/dE;
  // time
  std::complex <double> t (0, 0);
  // energy differences
  std::complex <double> wnm (0, 0);
  std::complex <double> wnnp (0, 0);
  std::complex <double> wnpm (0, 0);
  // coefficients
  std::complex <double> cm (0, 0);
  std::complex <double> cn (0, 0);
  std::complex <double> cn_term1 (0, 0);
  std::complex <double> cn_term2 (0, 0);
  std::complex <double> cn_diag (0, 0);
  std::complex <double> cn_offdiag (0, 0);
  double cn_tot;
  // complex numbers are dumb
  std::complex <double> C0 (0.0, 0.0);
  std::complex <double> C1 (1.0, 0.0);
  std::complex <double> NEGC1 (-1.0, 0.0);
  std::complex <double> CI (0.0, 1.0);
  std::complex <double> NEGCI (0.0, -1.0);

  // unpack params a bit
  int Nk = p->Nk;
  int Nc = p->Nc;
  int Ik = p->Ik;
  int Ic = p->Ic;
  int N = p->NEQ;
  double * energies = &(p->energies[0]);
  double * startWfn = &(p->startWfn[0]);

  // Create matrix of energy differences
  std::vector<std::complex <double>> Elr (Nk*Nc, std::complex <double> (0.0, 0.0));
  for (int ii = 0; ii < Nk; ii++) {
    for (int jj = 0; jj < Nc; jj++) {
      // array follows convention that first index is for QC state
      // e.g. Elr[i*Nc + j] = E_{ij}
      Elr[ii*Nc + jj] = std::complex <double> (energies[Ik + ii] - energies[Ic + jj], 0);
    }
  }
#ifdef DEBUG_ANALYTIC
  std::cout << std::endl;
  std::cout << "Energy gaps:" << std::endl;
  for (int ii = 0; ii < Nc*Nk; ii++) {
    std::cout << Elr[ii] << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
#endif

  // Create matrix of prefactors for each QC (n) state
  std::complex <double> pref;
  std::vector<std::complex <double>> prefQC (Nk*Nc, std::complex <double> (0.0, 0.0));
  for (int ii = 0; ii < Nk; ii++) {
    // V*c_l/(E_{lr} + i\kappa)
    pref = Vee*(std::complex <double> (startWfn[Ik + ii], startWfn[Ik + N + ii]));
    std::cout << startWfn[Ik + ii] << "," << pref << " ";
    for (int jj = 0; jj < Nc; jj++) {
      prefQC[ii*Nc + jj] = pref/(Elr[ii*Nc + jj] + CI*K);
    }
  }
#ifdef DEBUG_ANALYTIC
  std::cout << std::endl;
  for (int ii = 0; ii < Nc*Nk; ii++) {
    std::cout << prefQC[ii] << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
#endif

  // calculate wavefunction coefficients on electron-accepting side over time
  std::vector<std::complex <double>> crt (Nc*p->numOutputSteps, std::complex <double> (0.0, 0.0));

  int timeIndex = 0;
  for (std::complex <double> t = C0; std::real(t) <= p->tout;
       t += std::complex <double> (p->tout/p->numOutputSteps, 0.0), timeIndex++) {
    for (int ii = 0; ii < Nc; ii++) {
      // TODO add bit for multiple state terms
      for (int jj = 0; jj < Nk; jj++) {
	crt[timeIndex*Nc + ii] += prefQC[jj]*(exp(NEGCI*Elr[jj*Nc + ii]*t) - exp(NEGC1*K*t));
      }
    }
  }
  
  // calculate populations on electron-accepting side over time
  std::vector<double> Prt (Nc*p->numOutputSteps, 0.0);
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    for (int jj = 0; jj < Nc; jj++) {
      Prt[ii*Nc + jj] = pow(real(crt[ii*Nc + jj]), 2) + pow(imag(crt[ii*Nc + jj]), 2);
    }
  }

  if (isOutput(outs, "analytic_tcprob.out")) {
    std::ofstream output("analytic_tcprob.out");
    for (int ii = 0; ii <= p->numOutputSteps; ii++) {
      output << p->times[ii];
      for (int jj = 0; jj < Nc; jj++) {
	output << " " << Prt[ii*Nc + jj];
	output << " " << real(crt[ii*Nc + jj]) << " " << imag(crt[ii*Nc + jj]);
      }
      output << std::endl;
    }
    output.close();
  }

  return;
}
