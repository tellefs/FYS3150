#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>
using namespace arma;
using namespace std;

class Solver
{
public:
    Solver();
    void forwardEulerMethod(mat pos, mat vel, double force, double N );
    void verletMethod(mat pos, mat vel, double force_factor, double h, int N);
};

#endif // SOLVER_H
