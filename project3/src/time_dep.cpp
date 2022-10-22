//g++ time_dep.cpp penningtrap.cpp particle.cpp -o main -larmadillo

//Args: Number of timepoints, Time, Number of particles, interaction, omega_min, omega_max, omega_step 

#include "particle.hpp"
#include "penningtrap.hpp"

int main(int argc, char** argv)
{
  // Command line arguments
   if (argc != 8)
  {
    std::cout << "Expected 7 commandline arguments: " <<
    "\nint N (Timepoints)\ndouble T (Time)\nint n (number of particles)" 
    << "\nbool interaction (particle interaction true or false) \nstring method (Euler or RK4)\n";
    exit(1);
  }

  int N = atoi(argv[1]); //number of timepoints 
  double T = atoi(argv[2]); //Time
  int n = atoi(argv[3]); //number of particles
  double omega_min_in = atoi(argv[5]); // min frequency
  double omega_max_in = atoi(argv[6]); // max frequency
  double omega_step_in = atoi(argv[7]); // steps size for omega_v
  bool time_dependency = true; 
  bool interaction;
  std::string method = std::string(argv[5]);


  if (std::string(argv[4]) == "true")
  {
    interaction = true;
  }
  else if (std::string(argv[4]) == "false")
  {
    interaction = false;
  }
  else
  {
   std::cout << "argument 5 interaction has to be either true or false \n\nExiting...\n";
   exit(1);
  }
  
  // Variables
  double B0 = 96.5;
  double V0 = 2.41 * std::pow(10, 6);
  double d = 500;
  double m_ca = 40.077; // atomic mass unit
  double q_ca = 1; // elementary charge
  double dt = T / N; //Time interval
  arma::cube r; // to save nummerical postions 
  arma::cube v; // to save nummerical velocities 

  PenningTrap pt = PenningTrap(B0, V0, d, interaction, time_dependency); // Initialize PenningTrap

  if (n == 1 || n == 2)
  {
    arma::vec r0 = arma::vec({20, 0., 20});
    arma::vec v0 = arma::vec({0.,25,0.});
    Particle p_in = Particle(q_ca, m_ca, r0, v0);
    pt.add_particle(p_in);
  }
  
  if (n == 2)
  {
    arma::vec r0_2 = arma::vec({25, 25., 0});
    arma::vec v0_2 = arma::vec({0.,40, 5});
    Particle p_in_2 = Particle(q_ca, m_ca, r0_2, v0_2);
    pt.add_particle(p_in_2);
  }

  if (n > 2)  
  {
    pt.add_n_random_particles(n);
  }

// Checks what fratcion of particles are left in trap for amplitudes(0.1, 0.4, 0.7) different frquencies.
pt.parcticles_left_for_omega_v(dt, N, omega_min_in, omega_max_in, omega_step_in);


  return 0;
}