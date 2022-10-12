// Definitions for the functions in the PenningTrap class

#include "penningtrap.hpp"
#include "particle.hpp"

// Constructor
PenningTrap::PenningTrap(double B_0_in, double V_0_in, double d_in, std::vector<Particle> particle_objects_in)
{
  B_0 = B_0_in; 
  V_0 = V_0_in; 
  d = d_in; 
  particle_objects = particle_objects_in;
}


//Method that adds a particle to the Penning trap
void add_p(Particle p_in)
{
    particle_objects.push_back(p_in);
}

// Method that returns the ecternal electric field
arma::vec ext_el_field()
{

}

// Method that returns the ecternal magnetic field
arma::vec ext_mag_field()
{

}

// Method that returns the force due to the interaction among particles
arma::vec p_interaction_force()
{

}
