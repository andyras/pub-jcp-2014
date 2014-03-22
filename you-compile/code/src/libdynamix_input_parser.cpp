#include <fstream>
#include <iostream>
#include <string>
#include <map>

#include "libdynamix_input_parser.h"

//#define DEBUG_INPUT_PARSER

// Methods for the outputFile class

// Construct an outputFile object
outputFile::outputFile(char * fileName) {
 name = fileName;
}

// Get the name of an output file.
char * outputFile::getName() {
 return name;
}

// Toggle the creation of an output file
void outputFile::create() {
 createMe = true;
}

//void parseInput(const char * inputFile)

void assignOutputs(const char * inputFile, std::map<std::string, bool> &outputs) {
 std::string line;
 std::ifstream input;
 
 input.open(inputFile, std::ios::in);

#ifdef DEBUG_INPUT_PARSER
 std::cout << "Parsing file " << inputFile << "\n";
#endif

 if (input.good() == false) {
  std::cerr << "[assignOutputs] ERROR: bad input file " << inputFile << "\n";
 }

 // get first line
 getline(input, line);

 // skip lines until header specifying start of outputs
 while (getline(input, line)) {
#ifdef DEBUG_INPUT_PARSER
  std::cout << "line is " << line << "\n";
#endif
  if (line != "[[Output Files]]") {
   continue;
  }
  else {
   break;
  }
 }

 // read inputs until [[End]] header or EOF
 while (getline(input, line)) {
  if (line != "[[End]]") {

   // skip comments
   if (line.substr(0,1) == "#") {
#ifdef DEBUG_INPUT_PARSER
    std::cout << "Skipping comment line: " << line << "\n";
#endif
    continue;
   }

   // skip whitespace (space/tab) and blank lines
   else if ((line.find_first_not_of(' ') == std::string::npos)
            || (line.find_first_not_of('\t') == std::string::npos)) {
#ifdef DEBUG_INPUT_PARSER
    std::cout << "Skipping blank/whitespace line\n";
#endif
    continue;
   }

   // turn on outputs which are in this section of the input
   else {
#ifdef DEBUG_INPUT_PARSER
    std::cout << "Creating output file:  " << line << "\n";
#endif
    outputs[line] = true;
    if ((line.substr(line.length()-4, line.length()) != ".out")
        && (line.substr(line.length()-4, line.length()) != ".plt")) {
     std::cerr << "WARNING [" << __FUNCTION__ << "]: output file extension is not '.out' or '.plt'\n";
    }
   }
  }
 }
}
