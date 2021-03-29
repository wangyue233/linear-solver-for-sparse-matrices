# acse-5-assignment-awosile 
Ver. 1.0; 

Date: 1st Feb. 2020

Author: Wang, Y., Xiong, X., and Yu,T.

Contact: ty616@ic.ac.uk / +44 7724005943

# 1. Intro & Usage

This programme is used to solve dense & CSR matrices using Jacobi, 
Gauss-Seidel, LU(LDL) decomposition and Successive-over-Relaxation (SOR) methods. 
### CAUTION: No 0 element on diagonal MUST be satisfied.


# 2. Environment 
C++11 standard or above C++ environment.

e.g. Windows: Visual Studio 2017/2019; MacOS: XCode

# 3. Files
## 3.1 Header(.h) files
### 3.1.1 Matrix.h & Matrix.cpp
Class files to implement actions over dense matrices.

### 3.1.2 CSRMatrix.h & CSRMatrix.cpp
Class files to implement actions over CSR matrices.

## 3.2 Test files: main.cpp
A .cpp file that is used to run linear solver  Ax=b over the four methods writen in the program.

For matrix A, generated a random positive definite symmetric matrix from the solver.

For RHS vector b, use the default matrix or user input by her/himself.

Use g++ compiler,run command: g++ -std=c++11 main.cpp .


# 4. Implementation
## 4.1 Classes & Definitions
### 4.1.1 Class Matrix: manipulation over dense matrix
Parameters: 

rows (int, No. of rows)

cols (int, No. of columns)

preallocate (boolean, depicting whether the matrix has been set)

\*values (T pointer, elements in matrix)

---------------------------------------------------------------
Functions: 

2 Constructors (one taking in rows, cols, preallocate; the other taking rows, cols, \*values),

Destructor;

```cpp
virtual void printValues() // show on screen elements in matrix;

virtual void printMatrix() // show on screen elements in a matrix style;

virtual void  pos_def(int n) // to create an n*n matrix;

void  random_full() // fill the matrix with random value; 

void  matMatMult(Matrix<T>& mat_right, Matrix<T>& output) // matrix multiplication;

virtual void lu_fac(Matrix<T>& b, Matrix<double>& res) //  LU factorisation method for linear solver;

virtual void Jacobi(Matrix<T>& b, Matrix<double>& res) //  the Jacobi method for linear solver;

virtual void GS(Matrix<T>& b, Matrix<double>& res) // Gauss-Seidel method for linear solver;

virtual void SOR(Matrix<T>& b, Matrix<double>& res) //SOR method for linear solver;
  
T MAX(T \*a, int n) // compute the infinity norm of a matrix, outputing the maximum value;
```

### 4.1.2 Class CSRMatrix: manipulation over CSR matrix, inherited from Matrix class

Parameters: 

matrix class parameters,

nnzs (int, No. of non-zero values),

\*row_position (int, accumulated No. of non-zero values at that row),

\*col_position (int, respective column index of the non-zero value)

------------------------------------------------------------
Functions: 

2 Constructors (one taking in rows, cols, nnzs, preallocate; the other taking rows, cols, \*values,\*row_index, \*col_index);

Destructor;
```cpp
virtual void printValues() // show on screen elements in matrix;

virtual void printMatrix() // show on screen elements in a CSR matrix style;

virtual void pos_def(int n) // to create a pisitive-definite n\*n matrix;

void random_full() // fill the matrix with random value; 
  
virtual void lu_fac(Matrix<T>& b, Matrix<double>& res) // LU factorisation method for linear solver;

virtual void Jacobi(Matrix<T>& b, Matrix<double>& res) //  the Jacobi method for linear solver;

virtual void GS(Matrix<T>& b, Matrix<double>& res) // Gauss-Seidel method for linear solver;

virtual void SOR(Matrix<T>& b, Matrix<double>& res) //SOR method for linear solver.
```
## 4.2 Main function
```cpp
void solver(Matrix<T>& A, Matrix<T>& b, Matrix<double>& res) 
//a solver function that lets User to implement all the 4 methods to solve the matrix.

int main() //the MAIN function that lets User to decide:

//1. whether to use a dense or CSR matrix;
//2. the dimension of matrix A*;
//3. create a positive-definite matrix A through random generatation and print it out*;
//4. whether to input RHS vector b by User or random generatation,
//then solve the system by the 4 methods.

//* If the matrix's rank is above 20, a separate file is generated to store the matrix and it won't be printed.
```
