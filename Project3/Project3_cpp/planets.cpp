#include "planets.h"
#include <armadillo>

using namespace std;
using namespace arma;

Planets::Planets()
{ :Planets(vec, vec, double)

}
Planets::Planets(vec position, vec velocity, double mass ){
    m_position = position;
    m_velocity = velocity;
    m_mass = mass;
}
