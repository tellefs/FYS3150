#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

int error(int n, double u_arr)
{
    double *error_array = new double[n];
    double *v_array = new double[n];

    for(int i=0; i<n-1; i++){
        v_array[i+1] = 1 - (1-exp(-10))*(i+1) - exp(-10*(i+1));
        error_array[i+1] = log10(abs((v_array[i+1]-u_arr[i+1])/v_array[i+1]));




    }

    return 0;
}
