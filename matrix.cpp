#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>
using namespace std;

class matrix
{
    public:
        matrix( int rows, int cols);
        matrix(int * arrpt , int rows, int cols);
        ~matrix(){};
        void generateRandomMatrix( int rows, int cols , int rangeMax);
        void printMatrix();
        void setValues(int * arrpt);
        matrix operator+(const matrix & mat);
        matrix operator*(const matrix & mat);
        matrix transpose();
        matrix inverse  ();
        matrix gaussElimination();
        int lcm(int x, int y);
    //void print ();
        int rows;
        int cols;
    private:
        int matrank;
        int * arrayPt;
};

matrix matrix::operator+(const matrix & mat){
    try{
        if (matrix::rows != mat.rows && matrix::cols != mat.cols){
            throw(0);
        }
        else{
            int len = matrix::rows * matrix::cols;
            int * arrPt = new int[len];
            for (int i =0; i <len; i++){
                *(arrPt+i)= *(matrix::arrayPt+i) + *(mat.arrayPt+i);
            }
            matrix resultMatrix = matrix(matrix::rows, matrix::cols);
            resultMatrix.setValues(arrPt);
            return resultMatrix;
        }
    }
    catch(int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
    }
}

matrix matrix::operator*(const matrix & mat){
    try{
        if (matrix::cols != mat.rows){
            throw(0);
        }
        else{
            int len = matrix::rows * mat.cols;
            int * arrPt = new int[len];
            //page 38 AB = [ sum from j=1 to n (a_ij*b_jk)]
            int m = matrix::rows;
            int n = matrix::cols;
            int p = mat.cols;

            for (int i = 0; i < m; i++){
                for(int k = 0; k < p; k++){
                    int c_ik = 0;
                    for(int j = 0; j < matrix::cols; j++){
                        c_ik = c_ik + *(matrix::arrayPt + i*n + j) * *(mat.arrayPt + j*mat.cols +k);
                    }
                    * (arrPt + i*mat.cols + k) = c_ik;
                }
            }
            matrix resultMatrix = matrix(matrix::rows, mat.cols);
            resultMatrix.setValues(arrPt);
            return resultMatrix;
        }
    }
    catch(int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
    }
}

matrix::matrix(int * arrpt , int rows, int cols){
   matrix::arrayPt = arrpt;
}

matrix::matrix(int rows, int cols){
    matrix::rows = rows;
    matrix::cols = cols;
    matrix::arrayPt = NULL;
}

matrix matrix::transpose(){
    int len = matrix::rows*matrix::cols;
    int * arrPt = new int[len];
    for (int i = 0; i< matrix::rows; i++){
        for (int j=0; j<matrix::cols; j++){
            *(arrPt + i*matrix::rows + j) = *(matrix::arrayPt + j*matrix::cols + i);
        }
    }
    matrix resultMatrix = matrix(matrix::rows, matrix::cols);
    resultMatrix.setValues(arrPt);
    return resultMatrix;
}

int matrix::lcm(int x,int y){ //least common multiple
    int t;
    while (y != 0){
      t=y;
      y=x%y;
      x=t;
    }
    return x;
}

matrix matrix::gaussElimination(){
    for(int i=0;i<matrix::cols-1;i++){
        for(int j=i+1;j<matrix::rows;j++){
            int l=matrix::lcm(*(matrix::arrayPt+i*matrix::cols+i), *(matrix::arrayPt+j*matrix::cols+i));
            if(l!=0&&(*(matrix::arrayPt+i*matrix::cols+i)!=0&&*(matrix::arrayPt+j*matrix::cols+i)!=0)){
                //l=(a[i][i]*a[j][i])/l;
                l=(*(matrix::arrayPt+i*matrix::cols+i) * *(matrix::arrayPt+j*matrix::cols+i))/l;
                //int d1=l/a[i][i];
                int d1=l/ *(matrix::arrayPt+i*matrix::cols+i);
                //int d2=l/a[j][i];
                int d2=l/ *(matrix::arrayPt+j*matrix::cols+i);
                //a[j][i]=0;
                *(matrix::arrayPt+j*matrix::cols+i)=0;
                for(int k=i+1;k<matrix::cols;k++){
                    //a[j][k]=(d2*a[j][k])-(d1*a[i][k]);
                    *(matrix::arrayPt+j*matrix::cols+k) = d2* *(matrix::arrayPt+j*matrix::cols+k) - d1* *(matrix::arrayPt+i*matrix::cols+k);
                }
            }
        }
    }
    matrix resultMatrix = matrix(matrix::rows, matrix::cols);
    resultMatrix.setValues(matrix::arrayPt);
    return resultMatrix;
}

matrix matrix::inverse(){
}

void matrix::generateRandomMatrix( int rows, int cols, int rangeMax){
    int len = rows*cols;
    matrix::arrayPt = new int[len];
    for(int i=0; i<len; i++){
        int r = rand()%rangeMax;
        *(matrix::arrayPt+i) = r;
    }
}

void matrix::setValues(int * arrpt){
    matrix::arrayPt = arrpt;
}

void matrix::printMatrix(){
    int len = matrix::cols*matrix::rows;
    for (int i =0; i<len ; i++){
        if (i%cols==0){ cout <<"\n"; }
        cout << *(matrix::arrayPt+i) <<" ";
    }
}

int main(){
    cout <<"A";
    int rows =3;
    int cols =3;
    matrix A = matrix(rows, cols);
    A.generateRandomMatrix(rows, cols,3);
    A.printMatrix();

    cout <<"\nB";
    rows =3;
    cols =1;
    matrix B = matrix(rows, cols);
    B.generateRandomMatrix(rows, cols,3);
    B.printMatrix();

    cout <<"\nD";
    rows =3;
    cols =4;
    matrix D = matrix(rows, cols);
    D.generateRandomMatrix(rows, cols, 3);
    D.printMatrix();

    cout <<"\nH= gaussElimination(D)";
    matrix H = D.gaussElimination();
    H.printMatrix();

    /*
    cout <<"\nE";
    matrix E = A*D;
    E.printMatrix();

    cout <<"\nF=A*B";
    matrix F = A*B;
    F.printMatrix();

    cout <<"\nB";
    matrix B = matrix(rows, cols);
    B.generateRandomMatrix(rows, cols,3);
    B.printMatrix();
    cout <<"\nC=A+B";
    matrix C = matrix( rows, cols);
    C= A+B;
    C.printMatrix();
    cout <<"\nD=A*B";
    matrix D = matrix (rows, cols);
    D = A*B;
    D.printMatrix();
    cout <<"\nE=B*A";
    matrix E = matrix (rows, cols);
    E = B*A;
    E.printMatrix();

    cout <<"\nF=transpose(A)";
    matrix F = matrix (rows, cols);
    F = A.transpose();
    F.printMatrix();
    cout <<"\nG=inverse(A)";
    matrix G = matrix (rows, cols)
    return 0;*/
}
