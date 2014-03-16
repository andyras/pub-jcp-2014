#include <stdio.h>

#include "libdynamix_outputs.hpp"

/* Writes a scalar to a file */
void outputDScalar(const char * fileName, double scalar) {
  FILE * output;
  output = fopen(fileName, "w");
  fprintf(output, "%-.7g\n", scalar);
  fclose(output);
}

/* Writes a vector of length N to a file, one time point and scalar per row */
void outputTimeDVector(const char * fileName, double * time, double * vec, int n) {
  FILE * output;
  output = fopen(fileName, "w");
  for (int i = 0; i < n; i++) {
    fprintf(output, "%-.7g %-.7g\n", time[i], vec[i]);
  }
  fclose(output);
}

/* Writes transpose of a vector of length N to a file, all on one row*/
void outputDVectorT(const char * fileName, double * vec, int n) {
  FILE * output;
  output = fopen(fileName, "w");
  // put first value
  fprintf(output, "%-.7g", vec[0]);
  // put other values
  for (int i = 1; i < n; i++) {
    fprintf(output, " %-.7g", vec[i]);
  }
  fprintf(output, "\n");
  fclose(output);
}

/* Writes a vector of length N to a file, one scalar per row */
void outputDVector(const char * fileName, double * vec, int n) {
  FILE * output;
  output = fopen(fileName, "w");
  for (int i = 0; i < n; i++) {
    fprintf(output, "%-.7g\n", vec[i]);
  }
  fclose(output);
}

/* Writes a matrix with n rows and m columns to a file. */
void outputDMatrix(const char * fileName, double * mat, int n, int m) {
  FILE * output;
  output = fopen(fileName, "w");
  for (int i = 0; i < n; i++) {
    // write first value of row
    fprintf(output, "%-.7g", mat[i*m]);
    for (int j = 1; j < m; j++) {
      // write other row values
      fprintf(output, " %-.7g", mat[i*m + j]);
    }
    fprintf(output, "\n");
  }
  fclose(output);
}

/* Writes transpose of a matrix with n rows and m columns to a file. */
void outputDMatrixT(const char * fileName, double * mat, int n, int m) {
  FILE * output;
  output = fopen(fileName, "w");
  for (int j = 0; j < m; j++) {
    // write first value of row
    fprintf(output, "%-.7g", mat[j]);
    for (int i = 1; i < n; i++) {
      // write other row values
      fprintf(output, " %-.7g", mat[i*m + j]);
    }
    fprintf(output, "\n");
  }
  fclose(output);
}
