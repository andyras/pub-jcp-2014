#ifndef __RHS__
#define __RHS__


#include <iostream>
#include <vector>
#include <numeric>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>

#include "params.hpp"
#include "numerical.hpp"
#include "omp.h"

int RHS_WFN(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_WFN_SPARSE(realtype t, N_Vector y, N_Vector ydot, void * data);

void RELAX_KINETIC(int bandFlag, realtype * yp, realtype * ydotp, PARAMETERS * p);

void RELAX_RTA(int bandFlag, realtype * yp, realtype * ydotp, PARAMETERS * p);

int RHS_DM_RELAX(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_KINETIC(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_BLAS(realtype t, N_Vector y, N_Vector ydot, void * data);

double findDynamicMu(double pop, double T, int bandFlag, PARAMETERS * p);

double FDDBinarySearch(double lower, double upper, double T, double n,
    int bandFlag, PARAMETERS * p);

double FDDSum(double mu, double T, int bandFlag, PARAMETERS * p);

void FDD(double mu, double T, double * fdd, double * E, int N, double P);

void FDD_RTA(struct PARAMETERS * p, realtype * y, std::vector<double> & fdd, int flag);

double b13(double bm, double ekin, double ne, double K1, double K2, double K3, double X);

int RHS_DM_RTA(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_RTA_BLAS(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_dephasing(realtype t, N_Vector y, N_Vector ydot, void * data);

#endif
