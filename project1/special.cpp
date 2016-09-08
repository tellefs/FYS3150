#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

int special()
{
    int const n = 10;
    double h = 1./(n+1);
    double *a_tilde = new double[n];
    double *f_tilde = new double[n];
    double *f_arr = new double[n];
    double *u_arr = new double[n];

    for(int i=0; i<n; i++){
        f_arr[i] = -h*h*100*exp(-10*i*h);
    }

    a_tilde[0] = 2;
    f_tilde[0] = f_arr[0];

    //forward sub.
    for(int i=0; i<n; i++){
        a_tilde[i+1] = (i+1)/i;
        f_tilde[i+1] = f_arr[i+1] + (i*f_tilde[i])/(i+1);
    }

    //backward substitution for u_arr
    u_arr[n-1] = f_tilde[n-1] / a_tilde[n-1];
    for(int j=n-2; j>=0; j--){ //loop to calculate u
        u_arr[j] = (j/(j+1))*(f_tilde[j]-u_arr[j+1]);
    }

    return 0;
}
