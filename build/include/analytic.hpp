#ifndef __ANALYTIC__
#define __ANALYTIC__

#include <complex>
#include <vector>

#include "output.hpp"
#include "params.hpp"

void computeAnalyticOutputs(std::map<const std::string, bool> &outs,
    struct PARAMETERS * p);

#endif
