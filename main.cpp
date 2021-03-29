#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <ctime>
#include <fstream>


using namespace std;

template<class T>
void solver(Matrix<T>& A, Matrix<T>& b, Matrix<double>& res)
{
  A.GS(b, res);
  cout << "Result obtained from Gauss-Seidel" << endl;
  res.printMatrix();
  A.Jacobi(b, res);
  cout << "Result obtained from Jacobi" << endl;
  res.printMatrix();
  A.lu_fac(b, res);
  cout << "Result obtained from LU factorisation" << endl;
  res.printMatrix();
  A.SOR(b, res);
  cout << "Result obtained from SOR" << endl;
  res.printMatrix();
}

int main()
{

    int nnzs = 11;
    int choice;
    int rows;

    cout<<"Welcome to linear solver Ax = b: "<<endl;
    cout<<"Do you want to use dense matrix or sparse matrix structure for matrix A? "<<endl;
    cout<<"dense matrix: 1, sparse matrix: 2"<<endl;
    cout<<"Please input your choice: "<<endl;
    cin>>choice;
    cout<<"What is the rank of your matrix?"<<endl;
    cin>>rows;

    Matrix<double> *mat_A;

    // auto *dense_mat = new Matrix<double>(rows, cols, false);
    // auto *sparse_mat = new CSRMatrix<double>(rows,cols,nnzs,false);

   
    switch(choice) {
      case 1:
        cout<<"creating dense matrix"<<endl;
        mat_A = new Matrix<double>(rows, rows, false);
        mat_A->pos_def(rows);
        if (rows >= 20){
            fstream matrixA;
            matrixA.open("matrix_A.txt", fstream::out);
            if (matrixA.fail())
            {
                cout<<"Error opening file"<<endl;
                exit(0);
            }
            for(int i = 0;i<mat_A->rows;i++){
                for(int j = 0;j<mat_A->cols;j++){
                    matrixA<<mat_A->values[i*mat_A->rows+j]<<" ";
                }
                matrixA<<endl;
            }
            matrixA.close();
            cout << "This matrix is too big to show, thus we output it to a file called matrix_A.txt, please check the values there." <<endl;
        }else{
            mat_A->printMatrix();
        }
        cout<<"successfully created dense matrix"<<endl;
        break;
      case 2:
        cout<<"creating sparse matrix"<<endl;
        mat_A = new CSRMatrix<double>(rows,rows,nnzs,false);
        mat_A->pos_def(rows);
        if (rows >= 20){
            fstream matrixA;
            matrixA.open("matrix_A.txt", fstream::out);
            if (matrixA.fail())
            {
                cout<<"Error opening file"<<endl;
                exit(0);
            }
            int number = static_cast<CSRMatrix<double>*>(mat_A)->nnzs;
            matrixA<<"value:"<<endl;
            for(int i=0;i<number;i++){
                matrixA<<mat_A->values[i]<<" ";
            }
            matrixA<<endl;
            matrixA<<"col_index:"<<endl;
            for(int i=0;i<number;i++){
               matrixA<<static_cast<CSRMatrix<double>*>(mat_A)->col_index[i]<<" ";
            }
            matrixA<<endl;
            matrixA<<"row_position:"<<endl;
            for(int i=0;i<=mat_A->rows;i++){
            matrixA<<static_cast<CSRMatrix<double>*>(mat_A)->row_position[i]<<" ";
            }
            matrixA<<endl;

            matrixA.close();
           cout << "This matrix is too big to show, thus we output it to a file called matrix_A.txt, please check the values there." <<endl;
        }else{
            mat_A->printMatrix();
        }
        cout<<"successfully created sparse matrix"<<endl;
        break;
      default:
        cout<<"Sorry, I cannot understand."<<endl;
        exit(0);
    }

    auto *dense_vec = new Matrix<double>(rows, 1, true);
    auto *dense_res = new Matrix<double>(rows, 1, true);

    cout<<"If you want to generate your own RHS vector, please press 1. If you want us to generate it for you, please press 2."<<endl;
    int choice_b;
    double ele;
    cin>>choice_b;
    switch(choice_b) {
      case 1:
        for (int i=0;i<rows;i++){
        cout<<"Please enter your " <<i<<"th element."<<endl;
        cin>>ele;
        dense_vec->values[i] = ele;
      }
      cout<<"This is your RHS vector." <<endl;
      dense_vec->printMatrix();
      break;
      case 2:
        dense_vec->random_full();
      cout<<"This is the auto-generated RHS vector." <<endl;
          dense_vec->printMatrix();
        break;
      default:
        cout<<"Sorry, I cannot understand."<<endl;
        delete dense_vec;
      delete dense_res;
        exit(0);
    }
    solver(*mat_A, *dense_vec, *dense_res);

     // Don't forget to call delete (ie the destructor) once we're done, otherwise
     // we will leak memory from the new called inside our matrix
    delete mat_A;
    delete dense_vec;
    delete dense_res;

}
