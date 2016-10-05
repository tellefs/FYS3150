#include "solver.h"
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

Solver::Solver()
{

}

void Solver::forwardEulerMethod(mat pos, mat vel, double force, int N, double h){
    //function to integratre numerically using eulers forward method
    //the matrix pos should contain 3 N-long arrays representing x,y and z directions
    //same with matrix vel, only with velocity
    // force is the value of the force, N is array length and h is the steplength

    for(i=0; i<N; i++){

        double r = sqrt(pos(0,i)*pos(0,i) + pos(1,i)*pos(1,i) + pos(2,i)*pos(2,i));

        for(k=0; k<3; k++){
            vel(k, i+1) = vel(k,i) + force*pos(k,i)*h;
            pos(k,i+1) = pos(k, i) + vel(k, i+1)*h
        }

    }
}

void Solver::verletMethod(mat pos, vec vel, double force_factor, double h, int N){
    //function to integratre numerically using the verlet method
    //the matrix pos should contain 3 N-long arrays representing x,y and z directions
    //same with matrix vel, only with velocity
    // force is the value of the force, N is array length and h is the steplength

    for(i=1; i<N, i++){
        double r = sqrt(pos(0,i)*pos(0,i) + pos(1,i)*pos(1,i) + pos(2,i)*pos(2,i));
        double pi = 3.14; //how to find the value of pi? cmath?

        for(k=0; k<3; k++){
            force = force_factor*pos(k, i);
            vel_halfStep_negative = vel(k, i-1) + (vel(k, i) - vel(k,i-1))/2.;
            vel_halfStep_positive = vel_halfStep_negative + h*force;
            pos(k, i+1) = pos(k, i) + h + vel_halfStep_positive;
            force = force_factor*pos(k, i+1);
            vel(k, i+1) = vel_halfStep_positive + (h/2.)*force;
        }
    }

}
