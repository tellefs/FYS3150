#include <stdio.h>
#include <cmath>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void jacobiSolver(mat, int);
double maxOffDiag(mat, int, int*, int*);
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
    jacobiSolver(A, N);


   return 0;
}


//function to solve the matrix with the jacobi-algorithm
void jacobiSolver(mat A, int N){
    double epsilon = 1e-8;
    int k = 0; int l=0;

    while( maxOffDiag(A, N, &k, &l) > epsilon){
        double t,c,s;
        double tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if(tau<0){
            t = -1./(-tau + sqrt(1 + tau*tau));
        }
        else{
            t = 1./(tau + sqrt(1 + tau*tau));
        }
        c = 1./(1 + t*t);
        s = t*c;

        //now defining the new elements in the matrix A
        for(int i=0; i<0; i++){
            if(i!=k && i!=l){
                A(i,i) = A(i,i);
                A(i,k) = A(i,k)*c - A(i,l)*s;
                A(i,l) = A(i,l)*c + A(i,k)*s;

            }
        }
        A(k,k) = A(k,k)*c*c - 2*A(k,l)*c*s + A(l,l)*s*s;
        A(l,l) = A(l,l)*c*c + 2*A(k,l)*c*s + A(k,k)*s*s;
    }

A.print();
}

//function to find the maximum value off the diagonal of the matrix
//returns the value of maxOffVal and making pointers to the index of the max value of the diagonal
double maxOffDiag(mat A, int N, int* k, int* l){

    double maxOffVal = 0;
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++)
                if(j != i){
                    if(A(i,j)>maxOffVal){
                        maxOffVal = A(i,j);
                        *l = i;
                        *k = j;
                    }
                }
           }

    return maxOffVal;

}

//function to rotate the matrix
void rotate(){

}




