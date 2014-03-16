#ifndef __SPLINE__
#define __SPLINE__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>

// this struct holds the data for each point
struct Point {
  double x = 0.0;
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
  double d = 0.0;
};

class Spline {
  public:
    // default constructor
    Spline(const char * dataFile = NULL);
    // method to get the value at a certain point
    double value(double x);
    double getFirstX();
    double getLastX();
    // method to print the contents of the spline (for debug)
    void print();
  private:
    std::vector<Point> s;
};

bool comparePoints(const Point &pa, const Point &pb);

#endif
