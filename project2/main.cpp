#include <stdio.h>
#include <cmath>
#include <iostream>
#include "armadillo"

using namespace std;
using namespace arma;

void jacobiSolver(mat, int);
double maxOffDiag(mat, int, int, int);
void rotate();

int main()
{
    //defining constants
    double rho_min = 0;
    double rho_max = 100; //must test this value for different values of rho_max

    int N = 10;
    double h = (rho_max - rho_min)/N;
    double e_i = -1./(h*h); //the off diagonal entries to the tri-diagonal matrix

    vec V = zeros<vec>(N);
    vec rho = zeros<vec>(N);
    vec d = zeros<vec>(N);


    for(int i=0; i<N; i++){
        rho(i) = rho_min + i*h;
        V(i) = rho(i)*rho(i); //the potential energy
        d(i) = 2./(h*h) + V(i); //the diagonal entries to the matrix
    }

    //making the tri-diagonal matrix
    mat A = zeros<mat>(N,N);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if (i==j){
                A(i,j) = d(i);
            }
            if(i==j-1){
                A(i,j) = e_i;
            } else if (i==j+1) {
                A(i,j) = e_i;
            }
        }
    }

   return 0;
}

void jacobiSolver(mat A, int n){
    double epsilon = 1e-8;


}

//function to find the maximum value off the diagonal of the matrix
//returns the value of maxOffVal and making pointers to the index of the max value of the diagonal
double maxOffDiag(mat A, int N, int *k, int *l ){
    double maxOffVal = 0;
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++)
                if(j==i){
                    //Skip iteration
                }
                else{
                    if(A(i,j)>maxOffVal){
                        maxOffVal = A(i,j);
                        *l = i;
                        *k = j;
                    }
                }

        }
return maxOffVal;
}

void rotate(){

}




