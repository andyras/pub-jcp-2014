#ifndef __LIBDYNAMIX_INPUT_PARSER_H__
#define __LIBDYNAMIX_INPUT_PARSER_H__

#include <fstream>
#include <iostream>
#include <string>
#include <map>

#include "constants.hpp"
#include "output.hpp"
#include "params.hpp"
#include "spline.hpp"

// initiator for a string->bool map
std::map<const std::string, bool> initOutputMap();

const std::string makeConstString(std::string inputString);

// turns on output files in 'outs' map
void assignOutputs(const char * inputFile, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p);

void assignParams(const std::string inputFile, struct PARAMETERS * p);

bool fileExists(std::string fileName);

#endif
