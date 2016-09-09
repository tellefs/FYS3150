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

    for(int i=0; i<n; i++){
        f_arr[i] = h*h*100*exp(-10*(i+1)*h);
    }
    //initial conditions
    a_tilde[0] = 2;
    f_tilde[0] = f_arr[0];

    //forward sub.

    for(int i=0; i<n; i++){
        a_tilde[i+1] = (i+2.)/(i+1);
        f_tilde[i+1] = f_arr[i+1] + f_tilde[i]*(i+1.)/(i+2.);
        cout<<"f_tilde= "<<f_tilde[i]<<" a_tilde= "<<a_tilde[i]<<endl;
    }

    //backward sub.
    u_arr[n-1] = f_tilde[n-1] / a_tilde[n-1]; //initial conditions
    for(int i=n-2; i>=1; i--){ //loop to calculate u
        u_arr[i] = ((i)/(i+1.))*(f_tilde[i]+u_arr[i+1]);
        cout<<"u_arr = "<<u_arr[i]<<endl;
    }

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
