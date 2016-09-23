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
    double rho_max = 5; //must test this value for different values of rho_max

    int N = 50;
    double h = (rho_max - rho_min)/N;
    double e_i = -1./(h*h); //the off diagonal entries to the tri-diagonal matrix

    vec V = zeros<vec>(N);
    vec rho = zeros<vec>(N);
    vec d = zeros<vec>(N);


    for(int i=0; i<N; i++){
        rho(i) = rho_min + (i+1)*h;
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

    //A.print();
    double epsilon = 1e-8;
    int k = 0; int l=0;

    double max_off = maxOffDiag(A, N, &k, &l);


    while(max_off > epsilon){

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
        double a_ll, a_kk, a_ik, a_il;
        a_kk = A(k,k);
        a_ll = A(l,l);

        a_kk = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
        a_ll = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
        A(k,l) = 0.0;
        A(l,k) = 0.0;

        for(int i=0; i<N; i++){
            if(i!=k && i!=l){
                a_ik = A(i,k);
                a_il = A(i,l);
                A(i,k) = c*a_ik - s*a_il;
                A(k,i) = A(i,k);
                A(i,l) = c*a_il + s*a_ik;
                A(l,i) = A(i,l);

            }
        }

        max_off = maxOffDiag(A, N, &k, &l);

    }

    vec eigvals = zeros<vec>(N);
    //printing the diagonal elements aka. the eigenvalues
    for(int i=0; i<N; i++){
            eigvals(i) = A(i,i);
    }

    eigvals = sort(eigvals);

    cout << eigvals(0) << endl;

}

//function to find the maximum value off the diagonal of the matrix
//returns the value of maxOffVal and making pointers to the index of the max value of the diagonal
double maxOffDiag(mat A, int N, int* k, int* l){

    double maxOffVal = 0;

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++)
            if( (j != i) && ( abs(A(i,j)) > maxOffVal ) ){
                maxOffVal = abs(A(i,j));
                *l = i; //row
                *k = j; //col
            }
    }

    return maxOffVal;

}




