// Definitions for the functions in the PenningTrap class

#include "penningtrap.hpp"
#include "particle.hpp"

const double k_e = 1.38064903*std::pow(10, -23);

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_objects_in)
{
  double B_0 = B0_in;
  double V_0 = V0_in;
  double d = d_in;
  std::vector<Particle> p = particle_objects_in;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
  p.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
  //Analytic gradient:
  double E_x = -r[0]*V_0/(d*d);
  double E_y = -r[1]*V_0/(d*d);
  double E_z = 2*r[2]*V_0/(d*d);
  arma::vec E_ext = arma::vec((-E_x, -E_y, -E_z));
  return E_ext;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
  return arma::vec(0, 0, B_0);
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
  arma::vec r_diff = p[i].position() - p[j].position();

  arma::vec E = k_e * p[j].charge() * r_diff /
  std::pow(std::sqrt(r_diff[0]*r_diff[0] + r_diff[0]*r_diff[0] + r_diff[0]*r_diff[0]), 3);
  arma::vec F_int = p[i].charge() * E;

  return F_int;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i)
{
  arma::vec B = external_B_field(p[i].position());
  arma::vec E = external_E_field(p[i].position());
  arma::vec F_ext =  p[i].charge() * E
  + p[i].charge()*arma::cross( p[i].velocity(), B);

  return F_ext;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
  int n = p.size();
  arma::vec F_p = arma::vec(3);

  for (int j = 0; j < n; j++)
  {
    F_p = F_p + force_particle(i, j);
  }

  return F_p;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{
  arma::vec F_tot = total_force_particles(i) + total_force_external(i);
  return F_tot;
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
  int n = p.size();
  arma::mat k = arma::mat(4, 3);
  double dt2 = dt/2; //perform 4n less FLOPs by cutting divide in k2 and k3
  double dt6 = dt/6; //perform 2n less FLOPs by cutting divide in evolve
  arma::vec tmp_pos;
  arma::vec tmp_vel;


  for (int i = 0; i < n; i++)
  {
    //Temporary positions and velocities used to evolve
    tmp_pos = p[i].position();
    tmp_vel = p[i].velocity();

    //K_1
    arma::vec kr_1 = p[i].velocity();
    arma::vec kv_1 = total_force(i) / p[i].mass();

    //K_2
    p[i].position() += dt2*kr_1;
    p[i].velocity() += dt2*kv_1;
    arma::vec kr_2 = p[i].velocity();
    arma::vec kv_2 = total_force(i) / p[i].mass();

    //K_3
    p[i].position() += dt2*kr_2;
    p[i].velocity() += dt2*kv_2;
    arma::vec kr_3 = p[i].velocity();
    arma::vec kv_3 = total_force(i) / p[i].mass();

    //K_4
    p[i].position() += dt*kr_3;
    p[i].velocity() += dt*kv_3;
    arma::vec kr_4 = p[i].velocity();
    arma::vec kv_4 = total_force(i) / p[i].mass();

    //evolve
    p[i].position() = tmp_pos + dt6 * (kr_1 + 2*kr_2 + 2*kr_3 + kr_4);
    p[i].velocity() = tmp_vel + dt6 * (kv_1 + 2*kv_2 + 2*kv_3 + kv_4);
  }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
  int n = p.size();
  for (int i = 0; i < n; i++)
  {
    p[i].velocity() += p[i].velocity()[i]*dt;
    p[i].position() += p[i].velocity()[i]*dt;
  }
}
