#include <iostream>
#include <vector>
#include "CSRMatrix.h"
#include <cmath>

// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // use the unique pointer, no need for delete
      T *input_array = new T[this->nnzs];
      this->values.reset(input_array);
      int *row_array = new int[this->rows+1];
      this->row_position.reset(row_array);
      int *col_array = new int[this->nnzs];
      this->col_index.reset(col_array);

   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs)
{
   this->row_position.reset(row_position);
   this->col_index.reset(col_index);
}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // with unique pointer, no need for the delete command anymore
}

// Just print out the values in our values array
template <class T>
void CSRMatrix<T>::printMatrix()
{
   std::cout << "Printing matrix" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {
      std::cout << this->values[j] << " ";
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {
      std::cout << this->row_position[j] << " ";
   }
   std::cout << std::endl;
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {
      std::cout << this->col_index[j] << " ";
   }
   std::cout << std::endl;
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printValues(){
   
   T *print_max = new T[this->rows*this->cols]();
   // transform the CSR to matrix format
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         print_max[i*this->rows+this->col_index[val_index]]= this->values[val_index] ;

      }
   }
   std::cout << "Printing matrix" << std::endl;
   for (int j = 0; j< this->rows; j++)
   {
      std::cout << std::endl;
      for (int i = 0; i< this->cols; i++)
      {
         // We have explicitly used a row-major ordering here
         std::cout << print_max[i + j * this->cols] << " ";
      }
   }
   std::cout << std::endl;
   delete[] print_max;
}

// create a positive definite sparse matrix of size n * n
// With the value of <degree> to determine how sparse users want it to be
template<class T>
void CSRMatrix<T>::pos_def(int n)
{
   int degree = 2;
   this->cols = n;
   this->rows = n;
   srand(time(NULL));
// create a random sparse vector
// positive definite matrix =  vec * (vec).T + I
   T *ini_vec =  new T[n]();
   std::vector<int> position; // store the curret number of non zero elements
   int nnz = 0; // store the total number of non zero elements
   // create a random sparse vector
   for (int i = 0;i < n;i++){
      if (rand() % degree == 0){
         ini_vec[i] = (rand()%10)/7.0 + 1;
         nnz++;
         position.push_back(i);
      }
   }
   // positive definite matrix =  vec * (vec).T + I
   this->nnzs = nnz*nnz + (n-nnz);
   T *input_value = new T[this->nnzs];
   this->values.reset(input_value);
   int *col_array = new int[this->nnzs];
   this->col_index.reset(col_array);
   int *row_array = new int[this->rows+1];
   this->row_position.reset(row_array);

   this->row_position[0] = 0;
   int col_pos = 0;
   for (int i = 0;i<this->rows;i++){
      this->row_position[i+1] = this->row_position[i];
      //our solver cannot accept zero value of diagonal elemet
      // for each diagonal element, initialise with non zero value
      if (ini_vec[i]==0){
         this->values[col_pos] = this->cols * 6.5;
         // store the row and position in arrays
         this->col_index[col_pos] = i;
         this->row_position[i+1] += 1;
         col_pos++;
         continue;
      }
      for (int j = 0;j<position.size();j++){
         this->values[col_pos] = ini_vec[i]*ini_vec[position[j]];
          // store the row and position in arrays
         this->col_index[col_pos] = position[j];
         this->row_position[i+1] += 1;
         if (i==position[j])
            // for each diagonal element, initialise with non zero value
            this->values[col_pos] += this->cols * 6.5;
         col_pos++;
      }
      
   }
   delete[] ini_vec;
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(double *input, double *output)
{
   if (input == nullptr || output == nullptr)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   int val_counter = 0;
   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * input[this->col_index[val_index]];

      }
   }
}


// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output)
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
      if (this->rows != output.rows || mat_right.cols != output.cols)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      std::cerr << "OUTPUT HASN'T BEEN ALLOCATED" << std::endl;

   }

   // HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??
}
// Jacobi method of linear solver: Ax=b
// x: matrix res stores the value of x, should always be double
// attention: Jacobi method does not converge for every symmetric positive-definite matrix.
// it requires that the spectral radius of the iteration matrix is less than 1
// in order to simplify the computation, we stop the program once the iteration diverges
template<class T>
void CSRMatrix<T>::Jacobi(Matrix<T>& b, Matrix<double>& res)
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
   double *x0 = new double[n]();  //k-1th approximation of x
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
         int row_index = i;
         for (int val_index = this->row_position[row_index]; val_index < this->row_position[row_index+1]; val_index++)
         {
               if (row_index != this->col_index[val_index])
                  sum += this->values[val_index] * x0[this->col_index[val_index]];

         }
         for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
         {
               if (row_index == this->col_index[val_index])
                  xk[i] = (b.values[i] - sum) / this->values[val_index];
         }
      }
      
      aaa = MAX(xk, n);
      bbb = MAX(x0, n);
      ccc = fabs(bbb - aaa);                     // testify the convergence
   } while (ccc>1e-6 && ccc<10*aaa);                          //if the iteration diverges, stop the program
   res.values.reset(xk);
   if(ccc>1e-6){
      double *manual = res.values.release();
      std::cout<<"Sorry, Jacobi method cannot solve this matrix."<<std::endl;
      std::cout<<"Please use other methods."<<std::endl;
      delete[] manual;
   }
   delete[] x0;
}

// LU factorisation of linear solver: AX=b
// matrix res stores the value of x, should always be double
//  if matrix A is positive definite symmetric matrix
//  it can be decomposed as A = LD(L.T)
//  where l is a lower triangular matrix and d is a diagonal matrix
template <class T>
void CSRMatrix<T>::lu_fac(Matrix<T>& b, Matrix<double>& res)
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
   double *d = new double[n]();
   double *l = new double[n*n]();
   double *y = new double[n]();
   double *x = new double[n]();

    // store diagonal element
   std::vector<T> dia;
   for (int i =0; i<this->rows;i++){
         for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++){
            if (i == this->col_index[val_index]){
               dia.push_back(this->values[val_index]);
            }
         }
   }

   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      double sum1 = 0;
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         int j = this->col_index[val_index]; //column index
         if (j<i){
            sum1 = 0;
            for (int k = 0 ; k < j ; k++)
               sum1 += l[i*n+k] * l[j*n+k] * d[k];
            l[i*n+j] = (this->values[val_index] - sum1) / d[j];
         }
      }
      double sum2 = 0;
      for (int k = 0; k < i ; k++)
            sum2 += l[i*n+k] * l[i*n+k] * d[k];
      // Minus diagonal elements
      d[i] = dia[i] - sum2;
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

// Gauss Seidel method of linear solver: Ax=b
// x: matrix res stores the value of x, should always be double
// attention: Gauss Seidel method does not converge for every symmetric positive-definite matrix.
// it requires that the spectral radius of the iteration matrix is less than 1
// in order to simplify the computation, we stop the program once the iteration diverges
template<class T>
void CSRMatrix<T>::GS(Matrix<T>& b, Matrix<double>& res)
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
   double *x0 = new double[n]();  //k-1th approximation of x
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
         double diag;
         // int row_index = i;
         for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
         {
               if (i > this->col_index[val_index])
                  sum1 += this->values[val_index] * xk[this->col_index[val_index]];
               if (i == this->col_index[val_index])
                  diag = this->values[val_index];
               if (i < this->col_index[val_index])
                  sum2 += this->values[val_index] * x0[this->col_index[val_index]];


         }

         xk[i] = (b.values[i] - sum1 - sum2) / diag;
      }
      
      aaa = MAX(xk, n);
      bbb = MAX(x0, n);
      ccc = fabs(bbb - aaa);                     // testify the convergence
   } while (ccc>1e-6 && ccc<10*aaa);                          //if the iteration diverges, stop the program
   res.values.reset(xk);
   if(ccc>1e-6){
      double *manual = res.values.release();
      std::cout<<"Sorry, Gauss method cannot solve this matrix."<<std::endl;
      std::cout<<"Please use other methods."<<std::endl;
      delete[] manual;
   }
   delete[] x0;
}
// Successive over-relaxation(SOR) method of linear solver: Ax=b
// x: matrix res stores the value of x, should always be double
// attention: SOR method does not converge for every symmetric positive-definite matrix.
// it requires that the spectral radius of the iteration matrix is less than 1
// in order to simplify the computation, we stop the program once the iteration diverges
template<class T>
void CSRMatrix<T>::SOR(Matrix<T>& b, Matrix<double>& res)
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
   double *x0 = new double[n]();  //k-1th approximation of x
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
         double diag;
         // int row_index = i;
         for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
         {
               if (i > this->col_index[val_index])
                  sum1 += this->values[val_index] * xk[this->col_index[val_index]];
               if (i == this->col_index[val_index])
                  diag = this->values[val_index];
               if (i < this->col_index[val_index])
                  sum2 += this->values[val_index] * x0[this->col_index[val_index]];


         }
         xk_g[i] = (b.values[i] - sum1 - sum2) / diag;
         xk[i] = (1 - w)*x0[i] + w*xk_g[i];
      }
      
      aaa = MAX(xk, n);
      bbb = MAX(x0, n);
      ccc = fabs(bbb - aaa);                     // testify the convergence
   } while (ccc>1e-6 && ccc<10*aaa);                          //if the iteration diverges, stop the program
   res.values.reset(xk);
   if(ccc>1e-6){
      double *manual = res.values.release();
      std::cout<<"Sorry, Gauss method cannot solve this matrix."<<std::endl;
      std::cout<<"Please use other methods."<<std::endl;
      delete[] manual;
   }
   delete[] x0;
   delete[] xk_g;
}

