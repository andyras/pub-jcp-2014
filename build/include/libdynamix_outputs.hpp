#ifndef __LIBDYNAMIX_OUTPUTS_H__
#define __LIBDYNAMIX_OUTPUTS_H__

/* Writes a scalar to a file */
void outputDScalar(const char * fileName, double scalar);

/* Writes a vector of length N to a file, one time point and scalar per row */
void outputTimeDVector(const char * fileName, double * time, double * vec, int n);

/* Writes transpose of a vector of length N to a file, all on one row*/
void outputDVectorT(const char * fileName, double * vec, int n);

/* Writes a vector of length N to a file, one scalar per row */
void outputDVector(const char * fileName, double * vec, int n);

/* Writes a matrix with n rows and m columns to a file. */
void outputDMatrix(const char * fileName, double * mat, int n, int m);

/* Writes transpose of a matrix with n rows and m columns to a file. */
void outputDMatrixT(const char * fileName, double * mat, int n, int m);

#endif
