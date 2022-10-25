// Command to compile and link code:
// g++ main.cpp penningtrap.cpp particle.cpp -o main -larmadillo -O2
// See github repo for more information on how to run:
// https://github.com/Fslippe/FYS4150/tree/main/project3

#include "particle.hpp"
#include "penningtrap.hpp"

// Args: Number of timepoints, Time, Number of particles, interaction, method, time_dependency, amplitude, frequency
int main(int argc, char** argv)
{
  // Command line arguments
   if (argc != 9)
  {
    std::cout << "Expected 6 commandline arguments: " <<
    "\nint N (Timepoints)\ndouble T (Time)\nint n (number of particles)" 
    << "\nbool interaction (particle interaction true or false) \nstring method (Euler or RK4)"
    << "\nbool time_dependency (true or false)\nAmplitude \nfrequency (MHz)\n";
    exit(1);
  }

  int N = atoi(argv[1]); //number of timepoints 
  double T = atoi(argv[2]); //Time
  int n = atoi(argv[3]); //number of particles
  bool interaction;
  bool time_dependency;
  std::string method = std::string(argv[5]);

  // commandline arguments to bool statements
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
  
    if (std::string(argv[6]) == "true")
  {
    time_dependency = true;
  }
  else if (std::string(argv[6]) == "false")
  {
    time_dependency = false;
  }
  else
  {
   std::cout << "argument 7 time dependency has to be either true or false \n\nExiting...\n";
   exit(1);
  }

  // Variables
  double B0 = 96.5; // (u/microsecond/e)
  double V0 = 2.41 * std::pow(10, 6); //(u (micrometer)^2 /(microsecond)^2 / e)
  double d = 500; // micrometer
  double m_ca = 40.077; // atomic mass unit
  double q_ca = 1; // elementary charge
  double dt = T / N; //Time interval
  arma::cube r; // to save nummerical postions 
  arma::cube v; // to save nummerical velocities 

  PenningTrap pt = PenningTrap(B0, V0, d, interaction, time_dependency); // Initialize PenningTrap
  // initializing time dependent V_0
  if (time_dependency)
  {
    double f = atof(argv[7]);
    double omega = atof(argv[8]);
    pt.set_amplitude_and_frquency(f, omega);
  }
  
  // Chosen initial conditions for 1 and 2 particle simulations
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

  // For more than 2 particles, random inital positions and velocities are given
  if (n > 2)  
  {
    pt.add_n_random_particles(n);
  }

  pt.n = n; // initializing numer of particles to the Penning Trap

  // For analytic solution 
  if (std::string(method) == "Analytic")
  {
    arma::vec t = arma::linspace(0, T, N+1);
    arma::mat r_a = arma::mat(pt.analytic(t));
    r_a.save("data/r_a.dat");
    exit(0);
  }

  else
  {
    // Arrays to save r and v for nummerical solutions
    r = arma::cube(n, N+1,  3);
    v = arma::cube(n, N+1,  3);
  }
  
  // Adding initial conditions to arrays for saving
  for (int i = 0; i < n; i++)
  {
    r(arma::span(i),arma::span(0),arma::span::all) = pt.p[i].r;
    v(arma::span(i),arma::span(0),arma::span::all) = pt.p[i].v;
  }

  // Step forward in time depending on method
  // Using Euler
  if (std::string(method) == "Euler")
  {
    clock_t start = clock();

    for (int j = 1; j < N+1; j++)
    {
      pt.evolve_forward_Euler(dt);
      for (int i = 0; i < n; i++)
      {
        for (int k = 0; k < 3; k++)
        {
          r(i,j,k) = pt.p[i].r[k];
          v(i,j,k) = pt.p[i].v[k];
        }
      }
    }
    clock_t end = clock();
    double timeused = 1.*(end-start)/CLOCKS_PER_SEC;
    std::cout << "timeused = " << timeused << " seconds " << "\n";
  }     

  // Using RK4
  else if (std::string(method) == "RK4")
  {
    clock_t start = clock();
    for (int j = 1; j < N+1; j++)
    {
      pt.evolve_RK4(dt);
      for (int i = 0; i < n; i++)
      {
        for (int k = 0; k < 3; k++)
        {
          r(i,j,k) = pt.p[i].r[k];
          v(i,j,k) = pt.p[i].v[k];
        }
      }
    }
    clock_t end = clock();
    double timeused = 1.*(end-start)/CLOCKS_PER_SEC;
    std::cout << "timeused = " << timeused << " seconds " << "\n";
  }

  // Save r and v to dat files with names depending on method
  if (std::string(method) == "Euler")
    {
      r.save("data/r_Euler.dat");
      v.save("data/v_Euler.dat");
    }
  else if (std::string(method) == "RK4")
    {
      r.save("data/r_RK4.dat");
      v.save("data/v_RK4.dat");
    }

  return 0;
}
