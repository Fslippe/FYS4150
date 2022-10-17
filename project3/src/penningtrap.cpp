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
arma::vec PenningTrap::ext_el_field(double x, double y, double z, double V_0, double d)
{
  //Analytic gradient:
  double E_x = -x*V_0/(d*d);
  double E_y = -y*V_0/(d*d);
  double E_z = 2*z*V_0/(d*d);

  return arma::vec(-E_x, -E_y, -E_z);
}

// Method that returns the external magnetic field
arma::vec PenningTrap::ext_mag_field()
{
  return aram::vec(0, 0, B_0)
}

// Method that returns the force due to the interaction among particles
arma::vec PenningTrap::p_interaction_force(arma::vec r, arma::vec q)
{
  int n = r.size()
  for (size_t i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        E[i] = E[i] + k_e*q[j]* (r[i] - r[j]) / std::pow(arma::abs(r[i] - r[j]), 3)
      }
    }
  }

  arma::vec F =
}
