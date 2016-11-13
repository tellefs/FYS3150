#include <iostream>
#include <armadillo>
#include <time.h>
#include <fstream>
#include <cmath>


using namespace std;
using namespace arma;
double matrixEnergy(const mat &M, int L);
double matrixMagnetization(const mat &M, int L);
void updateBoundary(mat &M, int L);
void updateExpectationValues(mat M, int L);

mat randomMatrix(int L);
//mat constantMatrix(int L, double constant);


int main()
{

    int L = 20; //number of rows/cols in the matrix
     srand(time(NULL));

    /* Making matrix with random values 1 or -1.
     I have also made the matrix so that the matrix elements we are looking at are at index 1-L+1 (L elements in row and col.
     Index 0 and L+2 is a copy of the matrix element on the other side of the matrix (makes sense to me).
     This is so i dont have to check if the values are at the end or not later in the program, now I instead loop through the
     values from 1-L+1.
     */
    mat spinMatrix = randomMatrix(L);
    //mat spinMatrix = constantMatrix(L, 1);



    double beta_array [2] = {1.0, 1./2.4};
    const string tempFiles[] = {"dataT1.dat", "dataT24.dat"};
    int mccMax = 1e6;
    double w [17];



    //loop over temp
    for(int tempIndex=0; tempIndex<2; tempIndex++){
        ofstream outFile;
        string file = tempFiles[tempIndex];
        outFile.open(file.c_str(), ios::out);

        //go through MC cycles
        double beta = beta_array[tempIndex];
        w[0] = exp(8*beta);
        w[4] = exp(4*beta);
        w[8] = 1.;
        w[12] = exp(beta*-4);
        w[16] = exp(beta*-8);
        double meanE = 0;
        double meanE2 = 0;
        double meanM = 0;
        double meanM2 = 0;
        double meanMAbs = 0;
        double success = 0;

        double E = matrixEnergy(spinMatrix, L);
        double M = matrixMagnetization(spinMatrix, L);

        for(int mcc=0; mcc<mccMax; mcc++){
            //L*L metropolis algos
            success = 0;
            for(int sweep = 0; sweep < L*L; sweep++){
                int i = rand() % L + 1;
                int j = rand() % L + 1;

                double deltaEnergy = 2*spinMatrix(i,j)*
                        (spinMatrix(i+1, j) + spinMatrix(i-1, j) + spinMatrix(i, j+1) + spinMatrix(i, j-1));
                if(deltaEnergy <= 0){
                    spinMatrix(i,j) *= -1;
                    updateBoundary(spinMatrix, L);
                    E += (double) deltaEnergy;
                    M += (double) 2*spinMatrix(i,j);
                    success += 1;
                }
                else if(deltaEnergy > 0){
                    double wValue = w[(int)deltaEnergy+8];
                    double r = rand()/ (double)RAND_MAX;
                    if(r<= wValue){
                        spinMatrix(i,j) *= -1;
                        updateBoundary(spinMatrix, L);
                        E += (double) deltaEnergy;
                        M += (double) 2*spinMatrix(i,j);
                        success += 1;
                    }
                }
            }

            meanE += E;
            meanM += M;
            meanE2 += E*E;
            meanM2 += M*M;
            meanMAbs += fabs(M);

            outFile << meanE <<" "<< meanE2<<" "<< meanM <<" "<< meanM2 << " "<< success<<" "<< E << endl;
        }
        outFile.close();
        cout<<meanE/mccMax<<" "<<meanM/mccMax<<endl;

    }

    return 0;
}



//functions

double matrixEnergy(const mat &M, int L){
    double energy = 0;
    double J = 1;
    for(int i=1; i<L+1; i++){
        for(int j=1; j<L+1; j++){
            double matrixValue = M(i,j);
            energy += (matrixValue*M(i-1, j) + matrixValue*M(i, j-1));
        }
    }
    return -J*energy;
}


double matrixMagnetization(const mat &M, int L){
    double magnetization = 0;
    for(int i=1; i<L+1; i++){
        for(int j=1; j<L+1; j++){
            magnetization += M(i,j);
        }
    }
    return magnetization;
}


void updateBoundary(mat& M, int L){
    for(int i=1 ; i<+L+1; i++){
        M(i,L+1) = M(i, 0+1);
        M(L+1,i) = M(0+1, i);
        M(0,i) = M(L, i);
        M(i,0) = M(i, L);
    }
}

mat randomMatrix(int L){
    mat M = mat(L+2,L+2);

    for(int i=1; i<L+1; i++){
        for(int j=1; j<L+1; j++){
            int randSeed = rand() %2;
            if(randSeed == 1){
                M(i,j) = 1;
            }
            else{
                M(i,j) = -1;
            }
        }
    }
    M(0,0) = NAN;
    M(0,L+1) = NAN;
    M(L+1,0) = NAN;
    M(L+1,L+1) = NAN;

    updateBoundary(M, L);

    return M;
}

mat constantMatrix(int L, double constant){
    mat M = mat(L+2,L+2);
    for(int i=1; i<L+1; i++){
        for(int j=1; j<L+1; j++){
            M(i, j) = constant;
        }
    }
    M(0,0) = NAN;
    M(0,L+1) = NAN;
    M(L+1,0) = NAN;
    M(L+1,L+1) = NAN;

    updateBoundary(M, L);


    return M;
}
