#include <iostream>
#include <armadillo>
#include <time.h>
#include <fstream>
#include <cmath>
#include <mpi.h>

using namespace std;
using namespace arma;

const double RAND_MAX_INVERSE = 1.0/RAND_MAX;

//functions
double matrixEnergy(const mat &M, int L);
double matrixMagnetization(const mat &M, int L);
void updateBoundary(mat &matrix, int N, int i, int j);
void updateExpectationValues(mat M, int L);
mat randomMatrix(int L);
mat constantMatrix(int L, double constant);
void oneSweep(mat &matrix, int N, double &energy, double &mag, vec &w);
void monteCarlo(mat &matrix, int N, int maxMcc, double &meanE, double &meanM, double &meanE2, double &meanM2, vec &w);

int main(int nargs, char* args[])
{
    int numprocs, my_rank;
    //   MPI initializations
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    srand(my_rank);


    const string files[] = {"dataL40.dat", "dataL60.dat", "dataL100.dat", "dataL140.dat"};
    int L_array [] = {40, 60, 100, 140};

    for(int L_index=0; L_index<4; L_index++){
    //if(L_index==1) break; //REMOVE THIS WHEN RUNNING FULL PROGRAM!!!!!!!!!!!!!!
    clock_t start, finish;
    start = clock();

    int L = L_array[L_index];
    mat spinMatrix = randomMatrix(L);

    int mccMax = 1e6;
    vec w = vec(17);
    double norm = 1./mccMax;
    ofstream outFile;
    string file = files[L_index];
    if(my_rank == 0) {
        outFile.open(file.c_str(), ios::out);
    }

    //loop over temp
    for(double temp=2.2; temp<=2.35; temp+=0.01){
        double beta = 1./temp;
        w(0) = exp(8*beta);
        w(4) = exp(4*beta);
        w(8) = 1.;
        w(12) = exp(beta*-4);
        w(16) = exp(beta*-8);

        double meanE = 0;
        double meanM = 0;
        double meanE2 = 0;
        double meanM2 = 0;

        if(temp == 2.2) monteCarlo(spinMatrix, L, 1e3, meanE, meanM, meanE2, meanM2, w);
        meanE = 0;
        meanM = 0;
        meanE2 = 0;
        meanM2 = 0;
        monteCarlo(spinMatrix, L, mccMax, meanE, meanM, meanE2, meanM2, w);

        double Cv = (meanE2*norm - meanE*meanE*norm*norm)/(temp*temp);
        double chi = (meanM2*norm - meanM*meanM*norm*norm)/temp;

        double meanE_tot = 0;
        MPI_Reduce(&meanE, &meanE_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        double meanM_tot = 0;
        MPI_Reduce(&meanM, &meanM_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        double Cv_tot = 0;
        MPI_Reduce(&Cv, &Cv_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        double chi_tot = 0;
        MPI_Reduce(&chi, &chi_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(my_rank == 0) {
            outFile << meanE_tot*norm/numprocs << " " << meanM_tot*norm/numprocs << " " << Cv_tot/numprocs << " " << chi_tot/numprocs << " " << temp << endl;
        }

    }
    outFile.close();
    finish = clock();
    cout<<"L= "<< L_array[L_index]<<" done, rank "<< my_rank+1<<"/4. Time = "
       <<((double)(finish - start)/CLOCKS_PER_SEC)<<" sec."<<endl;
    }
    MPI_Finalize ();
    return 0;
}



//functions

void monteCarlo(mat &matrix, int N, int mccMax, double &meanE, double &meanM, double &meanE2, double &meanM2, vec& w){

    double norm = (double) 1./mccMax;
    double E = matrixEnergy(matrix, N);
    double M = matrixMagnetization(matrix, N);
    for(int mcc=0; mcc<mccMax; mcc++){
        oneSweep(matrix, N, E, M, w);
        meanE += E;
        meanM += fabs(M);
        meanE2 += E*E;
        meanM2 += M*M;

        }
}

void oneSweep(mat &matrix, int N, double &energy, double &mag, vec &w){
for(int sweep = 0; sweep < N*N; sweep++){
    int i = rand() % N + 1;
    int j = rand() % N + 1;

    double deltaEnergy = 2*matrix(i,j)*
            (matrix(i+1, j) + matrix(i-1, j) + matrix(i, j+1) + matrix(i, j-1));
    if(deltaEnergy <= 0){
        matrix(i,j) *= -1;
        if (i==1 || i==N || j==1 || j==N) updateBoundary(matrix, N, i, j);
        energy += (double) deltaEnergy;
        mag += (double) 2*matrix(i,j);
    }
    else if(deltaEnergy > 0){
        double wValue = w((int) deltaEnergy+8);
        double r = rand()*RAND_MAX_INVERSE;
        if(r<= wValue){
            matrix(i,j) *= -1;
            if (i==1 || i==N || j==1 || j==N) updateBoundary(matrix, N, i, j);
            energy += (double) deltaEnergy;
            mag += (double) 2*matrix(i,j);
        }
    }
}
}

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

void updateBoundary(mat &matrix, int N, int i, int j){
    if (i==1) {
        matrix(N+1,j) = matrix(i,j);
    }  else if (i==N) {
        matrix(0,j) = matrix(i,j);
    } if (j==1) {
        matrix(i,N+1) = matrix(i,j);
    } else if (j==N) {
        matrix(i,0) = matrix(i,j);
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
        if (i==1 || i==L || j==1 || j==L) updateBoundary(M, L, i, j);
        }
    }
    M(0,0) = NAN;
    M(0,L+1) = NAN;
    M(L+1,0) = NAN;
    M(L+1,L+1) = NAN;

    return M;
}

mat constantMatrix(int L, double constant){
    mat M = mat(L+2,L+2);
    for(int i=1; i<L+1; i++){
        for(int j=1; j<L+1; j++){
            M(i, j) = constant;
            if (i==1 || i==L || j==1 || j==L) updateBoundary(M, L, i, j);
        }
    }
    M(0,0) = NAN;
    M(0,L+1) = NAN;
    M(L+1,0) = NAN;
    M(L+1,L+1) = NAN;

    return M;
}
