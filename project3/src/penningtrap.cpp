// Definitions for the functions in the PenningTrap class

#include "penningtrap.hpp"
//#include "particle.hpp"

// Constructor
PenningTrap::PenningTrap(double B_0_in, double V_0_in, double d_in, std::vector<Particle> particle_objects_in)
{
  B_0 = B_0_in; 
  V_0 = V_0_in; 
  d = d_in; 
  particle_objects = particle_objects_in;
}


//Method that adds a particle to the Penning trap
void PenningTrap::add_particle(Particle p_in)
{
  particle_objects.push_back(p_in);
}

//Method that adds several particles in a vector to the Penning trap
void PenningTrap::add_particles(std::vector<Particle> particles_in)
{
  for (p_in : particles_in)
  {
    particle_objects.push_back(p_in);
  }
}

// Method that returns the external electric field
arma::vec PenningTrap::ext_el_field()
{

}

// Method that returns the external magnetic field
arma::vec PenningTrap::ext_mag_field()
{

}

// Method that returns the force due to the interaction among particles
arma::vec PenningTrap::p_interaction_force()
{

}
