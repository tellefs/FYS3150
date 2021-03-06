#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include "special.cpp"
#include "error.cpp"
#include "LU_decomp.h"
#include "armadillo"

using namespace std;
using namespace arma;

int main()
{
    int const n = 10;
    double h = 1./(n+1);
    double *b_arr = new double[n]; //should contain 1's
    double *a_arr = new double[n];   //should contain -2's
    double *c_arr = new double[n]; //should contain 1's
    double *u_arr = new double[n];
    double *f_arr = new double[n];
    double *f_tilde = new double[n];
    double *a_tilde = new double[n];


    //setting values of known arrays
    for(int i=0; i<n; i++){
        b_arr[i] = 1;
        a_arr[i] = -2;
        c_arr[i] = 1;
        f_arr[i] = -h*h*100*exp(-10*(i+1)*h);
        a_tilde[0] = a_arr[0];
        f_tilde[0] = f_arr[0];
    }

    clock_t start, finish; //to find the time the algorithm takes
    start = clock();

    //forward substitution for a_tilde and f_tilde
    for (int i=0; i < n-1; i++) {
        a_tilde[i+1] = a_arr[i+1] - (b_arr[i+1]*c_arr[i+1])/a_tilde[i];
        f_tilde[i+1] = f_arr[i+1] - f_tilde[i]*(c_arr[i]/a_tilde[i]);
    }

    //backward substitution for u_arr
    u_arr[n-1] = f_tilde[n-1] / a_tilde[n-1];
    for(int i=n-2; i>=0; i--){ //loop to calculate u
        u_arr[i] = (f_tilde[i] - b_arr[i]*u_arr[i+1])/a_tilde[i];
    }
    finish = clock();
    cout<<"time for gaussian elimination is "<<((double) (finish - start)/CLOCKS_PER_SEC)<<" sec."<<endl;

    //writing to file
    ofstream outFile;
    outFile.open("../../project1/data.dat", ios::out);
    outFile << 0 << endl;
    for (int i =0; i < n; i++) {
        outFile << u_arr[i] << endl;
    }
    outFile << 0 <<endl;
    outFile.close();

    //the rest of the programs to the project
    special(n);
    error(n, u_arr, h);
    LU_decomp(n, f_arr);


    return 0;

}

