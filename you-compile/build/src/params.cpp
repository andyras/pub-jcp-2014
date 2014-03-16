#include "params.hpp"

/* Get band index based on flag */
int bandStartIdx(int bandFlag, PARAMETERS * p) {
  if (bandFlag == CONDUCTION) {
    return p->Ik;
  }
  else if (bandFlag == VALENCE) {
    return p->Il;
  }
  else if (bandFlag == BRIDGE) {
    return p->Ib;
  }
  else if (bandFlag == QD_CONDUCTION) {
    return p->Ic;
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: unspecified band with flag " << bandFlag << std::endl;
    return 0;
  }
}

/* Get band index based on flag */
int bandEndIdx(int bandFlag, PARAMETERS * p) {
  if (bandFlag == CONDUCTION) {
    return p->Ik + p->Nk;
  }
  else if (bandFlag == VALENCE) {
    return p->Il + p->Nl;
  }
  else if (bandFlag == BRIDGE) {
    return p->Ib + p->Nb;
  }
  else if (bandFlag == QD_CONDUCTION) {
    return p->Ic + p->Nc;
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: unspecified band with flag " << bandFlag << std::endl;
    return 0;
  }
}

/* get number of states in band */
int bandNumStates(int bandFlag, PARAMETERS * p) {
  if (bandFlag == CONDUCTION) {
    return p->Nk;
  }
  else if (bandFlag == VALENCE) {
    return p->Nl;
  }
  else if (bandFlag == BRIDGE) {
    return p->Nb;
  }
  else if (bandFlag == QD_CONDUCTION) {
    return p->Nc;
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: unspecified band with flag " << bandFlag << std::endl;
    return 0;
  }
}
