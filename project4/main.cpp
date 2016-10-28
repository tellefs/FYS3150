#include <iostream>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;


int main()
{
    int L = 4;
    int n = L*L;

    /* Making matrix with random values 1 or -1.
     I have also made the matrix so that the matrix elements we are looking at are at index 1-L+1 (L elements in row and col.
     Index 0 and L+2 is a copy of the matrix element on the other side of the matrix (makes sense to me).
     This is so i dont have to check if the values are at the end or not later in the program, now I instead loop through the
     values from 1-L+1.
     */
    mat randomMatrix = mat(L+2,L+2);
    srand(time(NULL));

    for(int i=1; i<L+1; i++){
        for(int j=1; j<L+1; j++){
            int randSeed = rand() %2;
            if(randSeed == 1){
                randomMatrix(i,j) = 1;
            }
            else{
                randomMatrix(i,j) = -1;
            }

    }
    }


    for(int i=1 ; i<+L+1; i++){
        randomMatrix(i,L+1) = randomMatrix(i, 0+1);
        randomMatrix(L+1,i) = randomMatrix(0+1, i);
        randomMatrix(0,i) = randomMatrix(L, i);
        randomMatrix(i,0) = randomMatrix(i, L);
    }


    //pick random position in matrix
    int iRandom = rand() % L +1;
    int jRandom = rand() % L +1;

    //calculate energy at random position

    double J = 1;
    double beta = 1;

    double matrixValue = randomMatrix(iRandom, jRandom);
    double energyB = -J*(matrixValue*randomMatrix(iRandom-1, jRandom)
                        + matrixValue*randomMatrix(iRandom+1)
                        + matrixValue*randomMatrix(iRandom, jRandom-1)
                        + matrixValue*randomMatrix(iRandom, jRandom+1));

    //calculate trial energy

    double newMatrixValue = 0;
    if(randomMatrix(iRandom, jRandom) == 1){
        newMatrixValue = -1;

    }
    else{
        newMatrixValue = 1;
    }

    double energyTrial = -J*(newMatrixValue*randomMatrix(iRandom-1, jRandom)
                        + newMatrixValue*randomMatrix(iRandom+1)
                        + newMatrixValue*randomMatrix(iRandom, jRandom-1)
                        + newMatrixValue*randomMatrix(iRandom, jRandom+1));

    double deltaEnergy = energyTrial - energyB;

    if(deltaEnergy <= 0){
        randomMatrix(iRandom, jRandom) = newMatrixValue;
    }
    else if(deltaEnergy > 0){
        double w = exp(-beta*deltaEnergy);
        double r = rand()/ (double)RAND_MAX;
        cout<<r<<endl;
        if(r<= w){
            randomMatrix(iRandom, jRandom) = newMatrixValue;
        }
    }

    return 0;
}

