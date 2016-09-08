#include <iostream>

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
    double h = 1/(n+1);
    double *b_arr = new double[n-1]; //should contain 1's
    double *a_arr = new double[n];   //should contain -2's
    double *c_arr = new double[n-1]; //should contain 1's
    double *v_arr = new double[n];
    double *f_arr = new double[n];
    double *a_tilde = new double[n];

    for(int i=0; i<n; i++){
        b_arr[i] = 1;
        a_arr[i] = -2;
        c_arr[i] = 1;
        a_tilde[i] = a_arr[i] - (c_arr[i-1]*b_arr[i-1])/a_tilde[i-1];
        cout<<a_tilde[i]<<endl;
    }


    return 0;
}

