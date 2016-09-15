#ifndef SPECIAL
#define SPECIAL
#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

int special(int n)
{
    double h = 1./(n+1);
    double *a_tilde = new double[n];
    double *f_tilde = new double[n];
    double *f_arr = new double[n];
    double *u_arr = new double[n];
    double *i_array = new double[n];


    for(int i=0; i<n; i++){
        f_arr[i] = h*h*100*exp(-10*(i+1)*h);
        a_tilde[i] = (i+3.)/(i+2.);
        i_array[i] = (i+1.)/(i+2);
    }
    //initial conditions
    a_tilde[0] = 2;
    f_tilde[0] = f_arr[0];

    clock_t start, finish;
    start = clock();
    //forward sub.

    for(int i=0; i<n; i++){
        f_tilde[i+1] = f_arr[i+1] + f_tilde[i]*i_array[i];
    }

    //backward sub.
    u_arr[n-1] = f_tilde[n-1] / a_tilde[n-1]; //initial conditions
    for(int i=n-2; i>=0; i--){ //loop to calculate u
        u_arr[i] = i_array[i]*(f_tilde[i]+u_arr[i+1]);
    }

    finish = clock();
    cout<<"time for special gaussian is "<< ((double) (finish-start)/CLOCKS_PER_SEC)<<" sec."<<endl;

    //writing to file
    ofstream outFile;
    outFile.open("../../project1/special_data.dat", ios::out);
    outFile << 0 << endl;
    for (int i =0; i < n; i++) {
        outFile << u_arr[i] << endl;
    }
    outFile << 0 <<endl;
    outFile.close();

    return 0;
}
#endif
