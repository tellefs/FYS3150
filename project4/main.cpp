#include <iostream>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;


int main()
{
    int L = 4;
    int n = L*L;

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

    cout<<randomMatrix<<endl;

    for(int i=1 ; i<+L+1; i++){
        randomMatrix(i,L+1) = randomMatrix(i, 0+1);
        randomMatrix(L+1,i) = randomMatrix(0+1, i);
        randomMatrix(0,i) = randomMatrix(L, i);
        randomMatrix(i,0) = randomMatrix(i, L);
    }


    cout<<randomMatrix<<endl;

    //pick random position in matrix
    int iRandom = rand() % L +1;
    int jRandom = rand() % L +1;

    //calculate energy at random position

    double J = 1;

    double energy = -J*(randomMatrix(iRandom,jRandom)*randomMatrix(iRandom-1, jRandom)
                        + randomMatrix(iRandom,jRandom)*randomMatrix(iRandom+1)
                        + randomMatrix(iRandom,jRandom)*randomMatrix(iRandom, jRandom-1)
                        + randomMatrix(iRandom,jRandom)*randomMatrix(iRandom, jRandom+1));
    cout<<energy<<endl;



    return 0;
}

