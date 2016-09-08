#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

//function to make matrix
double **createMatrix(int n, int m, double initialValue) {
    double **matrix = new double*[n];
    for(int i=0; i<n; i++) {
        matrix[i] = new double[m];
        for(int j=0; j<m; j++) {
            matrix[i][j] = initialValue;
        }
    }

    return matrix;
}

int main()
{
    int const n = 10;
    double h = 1./(n+1);

    cout << "h=" << h << endl;
    double *b_arr = new double[n]; //should contain 1's
    double *a_arr = new double[n];   //should contain -2's
    double *c_arr = new double[n]; //should contain 1's
    double *u_arr = new double[n];
    double *f_arr = new double[n];
    double *f_tilde = new double[n];
    double *a_tilde = new double[n];

    for(int i=0; i<n; i++){
        b_arr[i] = 1;
        a_arr[i] = -2;
        c_arr[i] = 1;
        f_arr[i] = -h*h*100*exp(-10*i*h);
        a_tilde[0] = a_arr[0];
        f_tilde[0] = f_arr[0];
    }

    for (int i=0; i < n-1; i++) {
        a_tilde[i+1] = a_arr[i+1] - (b_arr[i]*c_arr[i])/a_tilde[i];
        f_tilde[i+1] = f_arr[i+1] - f_tilde[i+1]*(c_arr[i]/a_arr[i]);

        cout<<"f_tilde= "<<f_tilde[i]<<"  a_tilde= "<<a_tilde[i]<<endl;
    }

    u_arr[n-1] = f_tilde[n-1] / a_tilde[n-1];
    for(int j=n-2; j>0; j--){ //loop to calculate u
        u_arr[j] = (f_tilde[j] - b_arr[j]*u_arr[j+1])/a_tilde[j];
        cout << "u_Arr= " << u_arr[j]<<endl;
    }

    ofstream outFile;
    outFile.open("../../project1/data.dat", ios::out);

    for (int i =0; i < n; i++) {
        outFile << u_arr[i] << endl;
    }
    outFile.close();

    return 0;
}

