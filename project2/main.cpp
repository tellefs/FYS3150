#include <stdio.h>
#include <cmath>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void jacobiSolver(mat &, mat &, int);
double maxOffDiag(mat, int, int &, int &);
void rotate();
void test_eigenvectors(mat, int);
void test_maxOffDiag();

int main()
{
    //defining constants
    double rho_min = 0;
    double rho_max = 5;

    int N = 50;
    double h = (rho_max - rho_min)/N;
    double e_i = -1./(h*h); //the off diagonal entries to the tri-diagonal matrix
    double omega_r = 0.01;


    vec V = zeros<vec>(N);
    vec rho = zeros<vec>(N);
    vec d = zeros<vec>(N);


    for(int i=0; i<N; i++){
        rho(i) = rho_min + (i+1)*h;
        V(i) = rho(i)*rho(i); //the potential energy
        //V(i) = omega_r*omega_r*rho(i)*rho(i) + 1./rho(i); //the new potential for problem 2.d
        d(i) = 2./(h*h) + V(i); //the diagonal entries to the matrix
    }

    //making the tri-diagonal matrix A and the eigevector matrix R
    mat A = zeros<mat>(N,N);
    mat R = zeros<mat>(N,N);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if (i==j){
                A(i,j) = d(i);
                R(i,j) = 1.0;
            }
            if(i==j-1){
                A(i,j) = e_i;
            } else if (i==j+1) {
                A(i,j) = e_i;
            }
        }
    }



    jacobiSolver(A, R, N);


    test_eigenvectors(R, N);
    test_maxOffDiag();

    return 0;
}








//function to solve the matrix with the jacobi-algorithm
void jacobiSolver(mat &A, mat &R, int N){

    double epsilon = 1e-8;
    int k, l;
    double max_off = maxOffDiag(A, N, k, l);


    while(max_off > epsilon){

        double t,c,s;

        double tau = (A(l,l) - A(k,k))/(2*A(k,l));

        if(tau<0){
            t = -1./(-tau + sqrt(1 + tau*tau));
        }
        else{
            t = 1./(tau + sqrt(1 + tau*tau));
        }

        c = 1./sqrt(1 + t*t);
        s = t*c;

        //now defining the new elements in the matrix A
        double a_ll, a_kk, a_ik, a_il, r_ik, r_il;
        a_kk = A(k,k);
        a_ll = A(l,l);

        A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
        A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
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
            r_ik = R(i,k);
            r_il = R(i,l);
            R(i,k) = c*r_ik - s*r_il;
            R(i,k) = c*r_il + s*r_ik;
        }

        max_off = maxOffDiag(A, N, k, l);

    }



    vec eigvals = zeros<vec>(N);
    //sorting the eigenvalues
    for(int i=0; i<N; i++){
        eigvals(i) = A(i,i);
    }
    eigvals = sort(eigvals);

    cout <<"Lowest eigenvalue is: "<< eigvals(0) << endl; //printing the lowest eigenvalue
}


void test_maxOffDiag(){
    mat X = zeros<mat>(3,3);
    double s = 0;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            X(i,j) = s+1;
            s++;

        }
    }
    int k,l;
    double maxVal = maxOffDiag(X, 3, k, l);
    if((maxVal != 8) || (k!=2) || (l!=1)){
        cout<<"fuck off, your maxOffDiag function is wrong and it sucks, get it together Gary"<<endl;
        cout << "maxVal = " << maxVal << endl;
        cout << "k = " << k << endl;
        cout << "l = " << l << endl;
    }

}







//function to find the maximum value off the diagonal of the matrix
//returns the value of maxOffVal and making pointers to the index of the max value of the diagonal
double maxOffDiag(mat A, int N, int& k, int& l){

    double maxOffVal = 0;

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++)
            if( (j != i) && ( abs(A(i,j)) > maxOffVal ) ){
                maxOffVal = abs(A(i,j));
                k = i; //col
                l = j; //row
            }
    }
    return maxOffVal;
}





void test_eigenvectors(mat X, int n){
    //function for testing of the eigenvectors are orthogonal
    double tol = 1e-5;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i!=j && dot(X.row(i), X.row(j)) > tol ){
                cout<<"eigevectors "<<i<<", "<<j<<" are not orthogonal, dot product is: " << dot(X.row(i), X.row(j)) <<endl;
            }
        }
    }
}




