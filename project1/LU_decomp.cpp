#include "LU_decomp.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include "armadillo"


using namespace std;
using namespace arma;

int LU_decomp(int n, double *f_arr){
    mat A = mat(n,n, fill::eye)*-2;
    vec b = zeros<vec>(n);

    for (int i = 0; i < n; i++) {
        b(i) = f_arr[i];
    }


    //Making the A array which contains -2's at the diagonal and
    // -1's next to the diagonal
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j-1){
                A(i,j) = 1;
            } else if (i==j+1) {
                A(i,j) = 1;
            }
        }
    }

    clock_t start, finish;
    start = clock();

    mat L, U, P;
    lu(L,U,P,A);
    vec x = solve(trimatu(U), solve(trimatl(L), P*b) );

    finish = clock();
    cout<<"time for LU-decomp is "<< ((double) (finish-start)/CLOCKS_PER_SEC)<<" sec."<<endl;



    return 0;
}

