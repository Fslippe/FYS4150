// PenningTrap class
extern const double k_e;

#ifndef __penningtrap_hpp__
#define __tenningtrap_hpp__

#include <string>
#include <armadillo>
#include "particle.hpp"

class PenningTrap
{
  public:
    // Member variables
    double B_0; // magnetic field strength
    double V_0; // applied potetial
    double d; // characteristic dimension
    //double B_0 = B_0_in;
    std::vector<Particle> p; // conntains Particle objects in the Penning trap

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_objects);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);
};

#endif
