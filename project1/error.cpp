#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

int error(int n, double *u_arr, double h)
{
    double *error_array = new double[n];
    double *exact_array = new double[n];

    for(int i=0; i<n-1; i++){
        double x = (i+1)*h;
        exact_array[i+1] = 1 - (1-exp(-10))*x - exp(-10*x); //the exact solution
        error_array[i+1] = log10(abs((u_arr[i+1]- exact_array[i+1])/exact_array[i+1]));
    }
    //writing to file
    ofstream outFile;
    outFile.open("../../project1/error_data.dat", ios::out);
    outFile << 0 << endl;
    for (int i =0; i < n; i++) {
        outFile << error_array[i] << endl;
    }
    outFile << 0 <<endl;
    outFile.close();


    return 0;
}
