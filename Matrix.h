#pragma once
#include <memory>
template <class T>
class Matrix
{
public:

   // constructor where we want to preallocate ourselves
   Matrix(int rows, int cols, bool preallocate);
   // constructor where we already have allocated memory outside
   Matrix(int rows, int cols, T *values_ptr);
   // destructor
   virtual ~Matrix();

   // Print out the values in our matrix
   virtual void printValues();
   virtual void printMatrix();

   // create a positive definite matrix of size n * n
   virtual void pos_def(int n);
   // fill the matrix with random value
   void random_full();
   

   // Perform some operations with our matrix
   void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);
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
   std::unique_ptr<T[]> values;
   int rows = -1;
   int cols = -1;

// We want our subclass to know about this
protected:
   bool preallocated = false;

// Private variables - there is no need for other classes
// to know about these variables
private:

   int size_of_values = -1;
};

// compute the infinity norm of a matrix
template<class T>
T MAX(T *a, int n);
