// Definitions for the functions in the PenningTrap class

#include "penningtrap.hpp"
//#include "particle.hpp"

// Constructor
// Constructor
PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_objects_in);
{
  double B_0 = B_0_in;
  double V_0 = V_0_in;
  double d = d_in;
  Particle particle_objects = particle_objects_in;
}
// Add a particle to the trap
void add_particle(Particle p_in);
{
  particle_objects.push_back(p_in);
}
// External electric field at point r=(x,y,z)
arma::vec external_E_field(arma::vec r);
{
  //Analytic gradient:
  double E_x = -r[0]*V_0/(d*d);
  double E_y = -r[1]*V_0/(d*d);
  double E_z = 2*r[2]*V_0/(d*d);

  return arma::vec(-E_x, -E_y, -E_z);
}

// External magnetic field at point r=(x,y,z)
arma::vec external_B_field(arma::vec r);
{
  return arma::vec(0, 0, B_0)
}

// Force on particle_i from particle_j
arma::vec force_particle(int i, int j);
{
  Particle p_i = particle_objects[i];
  Particle p_j = particle_objects[j];

  r_diff = p_i.position - p_j.position;
  E = k_e * p_j.charge * r_diff / std::pow(std::sqrt(r_diff[0]*r_diff[0] + r_diff[0]*r_diff[0] + r_diff[0]*r_diff[0]), 3);

  return p_i.charge * E;
}

// The total force on particle_i from the external fields
arma::vec total_force_external(int i);
{
  Particle p = particle_objects[i];
  arma::vec B = external_B_field(p.position);
  arma::vec E = external_E_field(p.position);
  arma::vec F_ext = p.charge * E + p.charge*arma::cross(p.velocity, B);

  return F_ext;
}

// The total force on particle_i from the other particles
arma::vec total_force_particles(int i);
{
  int n = (particle_objects.charge).size;
  arma::vec F_p = arma::vec(3);

  for (int j = 0; j < n; j++)
  {
    F_p = F_p + force_particle(i, j);
  }

  return F_p;
}

// The total force on particle_i from both external fields and other particles
arma::vec total_force(int i);
{
  arma::vec F_tot = total_force_particles(i) + total_force_external(i);
  return F_tot
}
// Evolve the system one time step (dt) using Runge-Kutta 4th order
void evolve_RK4(double dt);
{

}
// Evolve the system one time step (dt) using Forward Euler
void evolve_forward_Euler(double dt);
{

}
