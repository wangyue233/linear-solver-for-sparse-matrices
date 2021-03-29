#include <iostream>
// #include <memory>
#include <cmath>
#include "Matrix.h"

// Constructor - using an initialisation list here
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate): rows(rows), cols(cols), size_of_values(rows * cols), preallocated(preallocate)
{
   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // use the unique pointer, no need for delete anymore
      T *input_array = new T[this->size_of_values];
      this->values.reset(input_array);
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
Matrix<T>::Matrix(int rows, int cols, T *values_ptr): rows(rows), cols(cols), size_of_values(rows * cols)
{
   this->values.reset(values_ptr);
}

// destructor
template <class T>
Matrix<T>::~Matrix()
{
   // with unique pointer, no need for the delete command anymore
 
}

// Just print out the values in our values array
template <class T>
void Matrix<T>::printValues()
{
   // Check if our matrix has had space allocated to it
   if (this->values) {
      std::cout << "Printing values" << std::endl;
      for (int i = 0; i< this->size_of_values; i++)
      {
         std::cout << this->values[i] << " ";
      }
      std::cout << std::endl;
   }
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void Matrix<T>::printMatrix()
{
   // Check if our matrix has had space allocated to it
   if (this->values) {
      std::cout << "Printing matrix" << std::endl;
      for (int j = 0; j< this->rows; j++)
      {
         std::cout << std::endl;
         for (int i = 0; i< this->cols; i++)
         {
            // We have explicitly used a row-major ordering here
            std::cout << this->values[i + j * this->cols] << " ";
         }
      }
      std::cout << std::endl;
   }
}

// create a positive definite matrix of size n * n
template <class T>
void Matrix<T>::pos_def(int n){
   this->cols = n;
   this->rows = n;
   this->size_of_values = n * n;
   srand(time(NULL));
// create a random vector
// positive definite matrix =  vec * (vec).T + I
   T *ini_vec =  new T[n];
   for (int i = 0;i < n;i++){
      ini_vec[i] = (rand()%10)/7.0;
   }
   T *input_value = new T[this->size_of_values];
   this->values.reset(input_value);
   for (int i = 0; i<this->rows;i++){
      for(int j = 0;j<this->cols;j++){
         this->values[i*n+j] = ini_vec[i]*ini_vec[j];
         if(i==j){
            this->values[i*n+j]+=this->rows*1.5;
         }
      }
   }
   delete[] ini_vec;

}

// fill the matrix with random value
template <class T>
void Matrix<T>::random_full(){
   // Check if our matrix has had space allocated to it
   if (this->preallocated){
   srand(time(NULL));
   // fill the matrix with random values
   for (int i = 0;i < this->size_of_values;i++){
      this->values[i] =(T) (rand()%10)/3.0;
   }
   }
}

// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_right, Matrix<T>& output)
{

   // Check our dimensions match
   if (this->cols != mat_right.rows)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (output.values != nullptr)
   {
      // Check our dimensions match
      if (this->rows != output.rows || this->cols != output.cols)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      output.values = new T[this->rows * mat_right.cols];
      output.preallocated = true;
   }

   // Set values to zero before hand
   for (int i = 0; i < output.size_of_values; i++)
   {
      output.values[i] = 0;
   }

   // Now we can do our matrix-matrix multiplication
   // CHANGE THIS FOR LOOP ORDERING AROUND
   // AND CHECK THE TIME SPENT
   // Does the ordering matter for performance. Why??
   for(int i = 0; i < this->rows; i++)
   {
      for(int k = 0; k < this->cols; k++)
      {
         for(int j = 0; j < mat_right.cols; j++)
         {
               output.values[i * output.cols + j] += this->values[i * this->cols + k] * mat_right.values[k * mat_right.cols + j];
         }
      }
   }
}

// LU factorisation of linear solver: AX=b
// matrix res stores the value of x, should always be double
//  if matrix A is positive definite symmetric matrix
//  it can be decomposed as A = LD(L.T)
//  where l is a lower triangular matrix and d is a diagonal matrix
template <class T>
void Matrix<T>::lu_fac(Matrix<T>& b, Matrix<double>& res)
{
   int n = this->cols;
   // Check our dimensions match
   if (this->cols!=this->rows)
   {
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix A must be a square matrix"<<std::endl;
      return;
   }
   // Check our dimensions match
   if (b.rows!=this->rows){
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix b must have the same dimension with matrix A"<<std::endl;
      return;
   }

   //  initize matrix
   double *d = new double[n];
   double *l = new double[n*n];
   double *y = new double[n];
   double *x = new double[n];
   //decomposition
   for (int i = 0;i < n; i++)
   {
      double sum1 = 0;
      for (int j = 0; j < i ; j++)
      {
            sum1 = 0;
            for (int k = 0; k < j ; k++)
               sum1 += l[i*n+k] * l[j*n+k] * d[k];
            l[i*n+j] = (this->values[i*n+j] - sum1) / d[j];
      }
      double sum2 = 0;
      for (int k = 0; k < i ; k++)
            sum2 += l[i*n+k] * l[i*n+k] * d[k];
      d[i] = this->values[i*n+i] -sum2;
   }


   //Ly=b
   y[0] = b.values[0];
   for (int i = 1; i < n; i++)
   {
      double sum3 = 0;
      for (int k = 0; k < i ; k++)
            sum3 += l[i*n+k] * y[k];
      y[i] = b.values[i] - sum3;
   }

   //D(L.T)x=y
   x[n-1] = y[n - 1] / d[n - 1];
   for (int i = n - 2; i >= 0; i--)
   {
      double sum4 = 0;
      for (int k = i + 1; k < n; k++)
            sum4 += l[k*n+i] * x[k];
      x[i] = y[i] / d[i] - sum4;
   }
   res.values.reset(x);
   
   delete[] d;
   delete[] l;
   delete[] y;

}

// Jacobi method of linear solver: AX=b
// x: matrix res stores the value of x, should always be double
// attention: Jacobi method does not converge for every symmetric positive-definite matrix.
// it requires that the spectral radius of the iteration matrix is less than 1
// in order to simplify the computation, we stop the program once the iteration diverges
template<class T>
void Matrix<T>::Jacobi(Matrix<T>& b, Matrix<double>& res)
{
   int n = this->rows;
   // Check our dimensions match
   if (this->cols!=this->rows)
   {
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix A must be a square matrix"<<std::endl;
      return;
   }
   // Check our dimensions match
   if (b.rows!=this->rows){
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix b must have the same dimension with matrix A"<<std::endl;
      return;
   }
   double *xk = new double[n](); //kth approximation of x
   double *x0 = new double[n];  //k-1th approximation of x
   double aaa; //the infinity norm of the kth approximation of x
   double bbb; //the infinity norm of the k-1th approximation of x
   double ccc;  //ccc = |xk|-|xk-1|
   // iterated until it converges
   do {
      for (int i = 0; i<n; i++)
         x0[i] = xk[i];

      for (int i = 0; i<n; i++)
      {
         double sum = 0;
         for (int j = 0; j<n; j++)
         {
               if (j == i)
                  continue;
               else
                  sum += this->values[i*n+j] * x0[j];
         }
         xk[i] = (b.values[i] - sum) / this->values[i*n+i];
      }
      aaa = MAX(xk, n);
      bbb = MAX(x0, n);
      ccc = fabs(bbb - aaa);                     // testify the convergence
   } while (ccc>1e-3 && ccc<10*aaa);                          //if the iteration diverges, stop the program
   res.values.reset(xk);
   if(ccc>1e-3){
      double *manual = res.values.release();
      std::cout<<"Sorry, Jacobi method cannot solve this matrix."<<std::endl;
      std::cout<<"Please use other methods."<<std::endl;
      delete[] manual;
   }
   delete[] x0;
}


// compute the infinity norm of a matrix
template<class T>
T MAX(T *a, int n)
{
    T max = 0;
    for (int i = 0; i<n; i++)
    {
        if (fabs(a[i])>max)
            max = fabs(a[i]);
    }
    return max;
}
// Gauss Seidel method of linear solver: AX=b
// x: matrix res stores the value of x, should always be double
// attention: Gauss Seidel method does not converge for every symmetric positive-definite matrix.
// it requires that the spectral radius of the iteration matrix is less than 1
// in order to simplify the computation, we stop the program once the iteration diverges
template<class T>
void Matrix<T>::GS(Matrix<T>& b, Matrix<double>& res)
{
   int n = this->rows;
    // Check our dimensions match
   if (this->cols!=this->rows)
   {
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix A must be a square matrix"<<std::endl;
      return;
   }
   // Check our dimensions match
   if (b.rows!=this->rows){
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix b must have the same dimension with matrix A"<<std::endl;
      return;
   }
   double *xk = new double[n](); //kth approximation of x
   double *x0 = new double[n];  //k-1th approximation of x
   double aaa; //the infinity norm of the kth approximation of x
   double bbb; //the infinity norm of the k-1th approximation of x
   double ccc;  //ccc = |xk|-|xk-1|
   // iterated until it converges
   do {
      for (int i = 0; i<n; i++)
         x0[i] = xk[i];

      for (int i = 0; i<n; i++)
      {
         double sum1 = 0;
         double sum2 = 0;
         for (int j = 0; j<i; j++)
         {
            sum1 += this->values[i*n+j] * xk[j];
         }
         for (int j = i+1; j<n; j++)
         {
            sum2 += this->values[i*n+j] * x0[j];
         }
         xk[i] = (b.values[i] - sum1 - sum2) / this->values[i*n+i];
      }
      aaa = MAX(xk, n);
      bbb = MAX(x0, n);
      ccc = fabs(bbb - aaa);                     // testify the convergence
   } while (ccc>1e-3 && ccc<10*aaa);                          //if the iteration diverges, stop the program
   res.values.reset(xk);
   if(ccc>1e-3){
      double *manual = res.values.release();
      std::cout<<"Sorry, Gauss Seidel method cannot solve this matrix."<<std::endl;
      std::cout<<"Please use other methods."<<std::endl;
      delete[] manual;
   }
   delete[] x0;
}
// Successive over-relaxation(SOR) method of linear solver: AX=b
// x: matrix res stores the value of x, should always be double
// attention: SOR method does not converge for every symmetric positive-definite matrix.
// it requires that the spectral radius of the iteration matrix is less than 1
// in order to simplify the computation, we stop the program once the iteration diverges
template<class T>
void Matrix<T>::SOR(Matrix<T>& b, Matrix<double>& res)
{
   int n = this->rows;
   // Check our dimensions match
   if (this->cols!=this->rows)
   {
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix A must be a square matrix"<<std::endl;
      return;
   }
   // Check our dimensions match
   if (b.rows!=this->rows){
      std::cout<<"Input dimensions for matrices don't match"<<std::endl;
      std::cout<<"matrix b must have the same dimension with matrix A"<<std::endl;
      return;
   }
   double *xk = new double[n](); //kth approximation of x
   double *x0 = new double[n];  //k-1th approximation of x
   double *xk_g = new double[n](); //Gauss-Seidel iterate
   double aaa; //the infinity norm of the kth approximation of x
   double bbb; //the infinity norm of the k-1th approximation of x
   double ccc;  //ccc = |xk|-|xk-1|
   double w = 1.4;  // relaxation factor 
   // iterated until it converges
   do {
      for (int i = 0; i<n; i++)
         x0[i] = xk[i];

      for (int i = 0; i<n; i++)
      {
         double sum1 = 0;
         double sum2 = 0;
         for (int j = 0; j<i; j++)
         {  
            sum1 += this->values[i*n+j] * xk[j];
         }
         for (int j = i+1; j<n; j++)
         {  
            sum2 += this->values[i*n+j] * x0[j];
         }
         xk_g[i] = (b.values[i] - sum1 - sum2) / this->values[i*n+i];
         // xk_g[i] = (b[i] - sum1 - sum2) / a[i][i];
         xk[i] = (1 - w)*x0[i] + w*xk_g[i];
      }
      aaa = MAX(xk, n);                           
      bbb = MAX(x0, n);
      ccc = fabs(bbb - aaa);                     // testify the convergence
   } while (ccc>1e-3 && ccc<10*aaa);                          //if the iteration diverges, stop the program
   res.values.reset(xk);
   if(ccc>1e-3){
      double *manual = res.values.release();
      std::cout<<"Sorry, SOR method cannot solve this matrix."<<std::endl;
      std::cout<<"Please use other methods."<<std::endl;
      delete[] manual;
   }
   delete[] x0;
   delete[] xk_g;
}
