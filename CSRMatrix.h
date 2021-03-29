#pragma once
#include "Matrix.h"

template <class T>
class CSRMatrix: public Matrix<T>
{
public:

   // constructor where we want to preallocate ourselves
   CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
   // constructor where we already have allocated memory outside
   CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index);
   // destructor
   ~CSRMatrix();

   // create a positive definite sparse matrix of size n * n
   virtual void pos_def(int n);
   // Print out the values in our matrix
   virtual void printMatrix();
   // Explicitly print out the values in values array as if they are a matrix
   virtual void printValues();

   // Perform some operations with our matrix
   void matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output);
   // Perform some operations with our matrix
   void matVecMult(double *input, double *output);
   // LU factorisation method for linear solver
   virtual void lu_fac(Matrix<T>& b, Matrix<double>& res);
   // the Jacobi method for linear solver
   virtual void Jacobi(Matrix<T>& b, Matrix<double>& res);
   // the  Gauss-Seidel method for linear solver
   virtual void GS(Matrix<T>& b, Matrix<double>& res);
   // the successive over relaxation method for linear solver
   virtual void SOR(Matrix<T>& b, Matrix<double>& res);

   // Explicitly using the C++11 nullptr here
   // use the unique pointer to store value, in order to avoid leak memory
   std::unique_ptr<int[]> row_position;
   std::unique_ptr<int[]> col_index;

   // How many non-zero entries we have in the matrix
   int nnzs=-1;

// Private variables - there is no need for other classes
// to know about these variables
private:
   
};
