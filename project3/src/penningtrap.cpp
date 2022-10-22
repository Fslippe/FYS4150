// Definitions for the functions in the PenningTrap class

#include "penningtrap.hpp"
#include "particle.hpp"

const double k_e = 1.38935333*std::pow(10, 5); 

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool interaction_in, bool time_dependency_in)
{
  B_0 = B0_in;
  V_0 = V0_in;
  d = d_in;
  interaction = interaction_in;
  time_dependency = time_dependency_in;
  V_0_d2 = V_0/(d*d);
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
  p.push_back(p_in);
}

// Add n random particles to the trap
void PenningTrap::add_n_random_particles(int n)
{
  double q = 1;
  double m = 40.077; // Calcium ion atomic mass
  arma::arma_rng::set_seed(1);
  for (int i = 0; i < n; i++)
  {
    arma::vec r = arma::vec(3).randn() * 0.1 * d;  // random initial position
    arma::vec v = arma::vec(3).randn() * 0.1 * d;  // random initial velocity
    Particle random_particle(q, m, r, v);
    add_particle(random_particle);
  }
  
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
  arma::vec E_ext;
  if (arma::norm(r) > d)
  {
    E_ext = arma::vec({0, 0, 0}); // sets external electric field to 0 in regions outside the trap
  }
  else
  {
     //Analytic gradient:
    double E_x = -r[0]*V_0_d2;
    double E_y = -r[1]*V_0_d2;
    double E_z = 2*r[2]*V_0_d2;
    E_ext = arma::vec({-E_x, -E_y, -E_z});
  }
  return E_ext;
}
//Method to set amplitude f and angular frequency omega_v
void PenningTrap::set_amplitude_and_frquency(double amplitude_in, double frequency_in)
{
  f = amplitude_in;
  omega_v = frequency_in;
}

// External electric field at point r=(x,y,z) with time dependence
arma::vec PenningTrap::external_E_field(arma::vec r, double t)
{
  arma::vec E_ext;
  if (arma::norm(r) > d)
  {
    E_ext = arma::vec({0, 0, 0}); // sets external electric field to 0 in regions outside the trap
  }
  else
  {
    double V_0_time_dep = V_0 * (1 + f * std::cos(omega_v * t));
     //Analytic gradient:
    double E_x = -r[0]*V_0_time_dep /(d*d);
    double E_y = -r[1]*V_0_time_dep / (d*d);
    double E_z = 2*r[2]*V_0_d2;
    E_ext = arma::vec({-E_x, -E_y, -E_z});
  }
  return E_ext;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
  arma::vec B_ext;
  if (arma::norm(r) > d)
  {
    B_ext =  arma::vec({0, 0, 0}); // sets external magnetic field to 0 in regions outside the trap
  }
  else
  {
    B_ext = arma::vec({0, 0, B_0});
  }
  return B_ext;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
  arma::vec r_diff = p[i].position() - p[j].position();

  arma::vec E_i = k_e * p[j].charge() * r_diff /
  std::pow(std::sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]), 3);
  arma::vec F_int = p[i].charge() * E_i;
  
  return F_int;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i)  
{
  arma::vec E;
  if (time_dependency == true)
  {
    E = external_E_field(p[i].position(), time);
  }
  else
  {
    E = external_E_field(p[i].position());
  }
  
  arma::vec B = external_B_field(p[i].position());
  arma::vec F_ext = p[i].charge() * E + p[i].charge()*arma::cross( p[i].velocity(), B);
  return F_ext;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
  int n = p.size();
  arma::vec F_p = arma::zeros(3);
  for (int j = 0; j < n; j++)
  {
    if (i != j)
    {
      F_p += force_particle(i, j);
    }
  }
  return F_p;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{
  arma::vec F_tot;

  if (interaction)
  {
    F_tot = total_force_particles(i) + total_force_external(i);
  }
  else
  {
    F_tot = total_force_external(i);
  }
  
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
  arma::vec kr_1;
  arma::vec kv_1;
  arma::vec kr_2;
  arma::vec kv_2;
  arma::vec kr_3;
  arma::vec kv_3;
  arma::vec kr_4;
  arma::vec kv_4;

  for (int i = 0; i < n; i++)
  {
    //Temporary positions and velocities used to evolve
    tmp_pos = p[i].position();
    tmp_vel = p[i].velocity();

    //K_1
    kr_1 = p[i].v;
    kv_1 = total_force(i) / p[i].mass();

    //K_2
    p[i].r += dt2*kr_1;
    p[i].v += dt2*kv_1;
    kr_2 = p[i].v;
    kv_2 = total_force(i) / p[i].mass();

    //K_3
    p[i].r += dt2*kr_2;
    p[i].v += dt2*kv_2;
    kr_3 = p[i].v;
    kv_3 = total_force(i) / p[i].mass();

    //K_4
    p[i].r += dt*kr_3;
    p[i].v += dt*kv_3;
    kr_4 = p[i].v;
    kv_4 = total_force(i) / p[i].mass();

    //evolve
    p[i].r = tmp_pos + dt6 * (kr_1 + 2*kr_2 + 2*kr_3 + kr_4);
    p[i].v = tmp_vel + dt6 * (kv_1 + 2*kv_2 + 2*kv_3 + kv_4);
  }
  time += dt;
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
  int n = p.size();
  arma::vec tmp_vel;
  arma::vec v; 
  arma::vec r; 

  for (int i = 0; i < n; i++)
  {
    tmp_vel = p[i].velocity();
    p[i].v += total_force(i) / p[i].mass()*dt;
    p[i].r += tmp_vel*dt;
  }
  time += dt;
}

arma::mat PenningTrap::analytic(arma::vec t)
{
  arma::vec x;
  arma::vec y;
  arma::vec z;
  int n = t.size();
  arma::mat r = arma::mat(n, 3);

  if (p.size() == 1)
  {
    double omega_0 = p[0].charge() * B_0 / p[0].mass();
    double omega_z = std::sqrt( (2 * p[0].charge() * V_0) / (p[0].mass() * d*d) );
    double omega_plus = 0.5 * (omega_0 + std::sqrt(omega_0 * omega_0 - 2* omega_z * omega_z));
    double omega_minus =  0.5 * (omega_0 - std::sqrt(omega_0 * omega_0 - 2* omega_z * omega_z));
  
    double x_0 = p[0].position()[0];
    double y_0 = p[0].position()[1];
    double z_0 = p[0].position()[2];

    // The x and z component of v are 0. only using y component.
    double v_0 = p[0].velocity()[1];
   

    double A_plus = (v_0 + omega_minus * x_0) / (omega_minus - omega_plus);
    double A_minus = -(v_0 + omega_plus * x_0) / (omega_minus - omega_plus);

    x = A_plus * arma::cos(omega_plus *t) + A_minus * arma::cos(omega_minus *t);
    y = - A_plus * arma::sin(omega_plus *t) - A_minus * arma::sin(omega_minus *t);
    z = z_0 * arma::cos(omega_z * t);
  }
  else
  {
    std::cout << "Too many particles! \n" << "The analytic solution is only designed for a single particle." << "\n";
  } 
  r.col(0) = x; r.col(1) = y; r.col(2) = z;
  return r;
}

//Returns number of particles still inside trap
int PenningTrap::particles_left_in_trap()
{
  int counter = 0;
  for (int i = 0; i < p.size(); i++)
  {
   if (arma::norm(p[i].position()) < d)
   {
    counter += 1;
   }
  }
  return counter;
}

//Creates a matrix with fraction of particles left in trap for different frequencies and writes it to a .dat file.
void PenningTrap::parcticles_left_for_omega_v(double dt, int N, double omega_min, double omega_max, double omega_step)
{
 int points = (omega_max - omega_min) /omega_step;
  arma::mat frac_p_left = arma::mat(3,points);
  int j = 0; int k = 0;
  for (double f : {0.1, 0.4, 0.7})
  {
    for (double omega_v : arma::linspace(omega_min, omega_max, points))
    {
      set_amplitude_and_frquency(f, omega_v);
      for (int i = 0; i < N; i++)
      {
        evolve_RK4(dt);
      }
      frac_p_left(j,k) = particles_left_in_trap()/p.size(); // fraction of particles left in trap
      k += 1;
    }
    j += 1;
  }
  frac_p_left.save("data/frac_p_left.dat");
}
