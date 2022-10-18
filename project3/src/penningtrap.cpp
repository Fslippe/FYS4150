// Definitions for the functions in the PenningTrap class

#include "penningtrap.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_objects_in);
{
  double B_0 = B_0_in;
  double V_0 = V_0_in;
  double d = d_in;
  std::vector<Particle> particle = particle_objects_in;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in);
{
  particle.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r);
{
  //Analytic gradient:
  double E_x = -r[0]*V_0/(d*d);
  double E_y = -r[1]*V_0/(d*d);
  double E_z = 2*r[2]*V_0/(d*d);

  return arma::vec(-E_x, -E_y, -E_z);
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r);
{
  return arma::vec(0, 0, B_0)
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j);
{
  r_diff = particle[i].position - particle[j].position;

  arma::vec E = k_e * particle[j].charge * r_diff /
  std::pow(std::sqrt(r_diff[0]*r_diff[0] + r_diff[0]*r_diff[0] + r_diff[0]*r_diff[0]), 3);

  return particle[i].charge * E;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i);
{
  arma::vec B = external_B_field(particle[i].position);
  arma::vec E = external_E_field(particle[i].position);
  arma::vec F_ext =  particle[i].charge * E
  + particle[i].charge*arma::cross( particle[i].velocity, B);

  return F_ext;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i);
{
  int n = (particle.charge).size;
  arma::vec F_p = arma::vec(3);

  for (int j = 0; j < n; j++)
  {
    F_p = F_p + force_particle(i, j);
  }

  return F_p;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i);
{
  arma::vec F_tot = total_force_particles(i) + total_force_external(i);
  return F_tot;
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt);
{
  int n = (particle.charge).size;
  arma::mat k = arma::mat(4, 3);
  double dt2 = dt/2; //perform 4n less FLOPs by cutting divide in k2 and k3
  double dt6 = dt/6; //perform 2n less FLOPs by cutting divide in evolve
  arma::vec tmp_pos;
  arma::vec tmp_vel;


  for (int i = 0; i < n; i++)
  {
    //Temporary positions and velocities used to evolve
    tmp_pos = particle[i].position;
    tmp_vel = particle[i].velocity;

    //K_1
    arma::vec kr_1 = particle[i].velocity;
    arma::vec kv_1 = total_force(i) / particle[i].mass;

    //K_2
    particle[i].position += dt2*kr_1;
    particle[i].velocity += dt2*kv_1;
    arma::vec kr_2 = particle[i].velocity;
    arma::vec kv_2 = total_force(i) / particle[i].mass;

    //K_3
    particle[i].position += dt2*kr_2
    particle[i].velocity += dt2*kv_2;
    arma::vec kr_3 = particle[i].velocity;
    arma::vec kv_3 = total_force(i) / particle[i].mass;

    //K_4
    particle[i].position += dt*kr_3
    particle[i].velocity += dt*kv_3;
    arma::vec kr_4 = particle[i].velocity;
    arma::vec kv_4 = total_force(i) / particle[i].mass;

    //evolve
    particle[i].position = tmp_pos + dt6 * (kr_1 + 2*kr_2 + 2*kr_3 + kr_4);
    particle[i].velocity = tmp_vel + dt6 * (kv_1 + 2*kv_2 + 2*kv_3 + kv_4);
  }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt);
{
  for (int i = 0; i < n; i++)
  {
    particle[i].velocity += particle[i].velocity[i]*dt;
    particle[i].position += particle[i].velocity[i]*dt;
  }
}
