//g++ main.cpp penningtrap.cpp particle.cpp -o main -larmadillo

#include "particle.hpp"
#include "penningtrap.hpp"


int main(int argc, char** argv)
{
  if (argc != 6)
  {
    std::cout << "Expected 6 commandline arguments: " <<
    "int N (Timepoints) double T (Time) int n (number of particles)" 
    << "bool interaction (particle interaction true or false) string method (Euler or RK4)\n\nExiting...\n";
    exit(1);
  }

  int N = atoi(argv[1]); //number of timepoints 
  double T = atoi(argv[2]); //Time
  int n = atoi(argv[3]); //number of particles
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


  std::cout << "\nN:" << N << "\nT:" << T << "\nn:" << n << "\nInteraction:" <<interaction << "\n"; 
  
  double B0 = 96.5;
  double V0 = 2.41 * std::pow(10, 6);
  double d = 500;
  double m_ca = 40.077; // atomic mass unit
  double q_ca = 1; // elementary charge
  double dt = T / N; //Time interval

  PenningTrap pt = PenningTrap(B0, V0, d, interaction);

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
    
  }
  
  arma::cube r = arma::cube(n, N,  3);
  arma::cube v = arma::cube(n, N,  3);

  if (method == "Analytic")
  {
    arma::vec t = arma::linspace(0, T, N);
    arma::mat r_a = arma::mat(pt.analytic(t));
    r_a.save("data/r_a.dat");
  }
  
  for (int i = 0; i < n; i++)
  {
    r(arma::span(i),arma::span(0),arma::span::all) = pt.p[i].r;
    v(arma::span(i),arma::span(0),arma::span::all) = pt.p[i].v;
  }

  for (int j = 1; j < N; j++)
  {
    if (method == "Euler")
    {
      pt.evolve_forward_Euler(dt);
    }
    else if (method == "RK4")
    {
      pt.evolve_RK4(dt);
    }
    for (int i = 0; i < n; i++)
    {
      r(arma::span(i),arma::span(j),arma::span::all) = pt.p[i].r;
      v(arma::span(i),arma::span(j),arma::span::all) = pt.p[i].v;
    }
  }

  if (method == "Euler")
    {
      r.save("data/r_Euler.dat");
      v.save("data/v_Euler.dat");
    }
    else if (method == "RK4")
    {
      r.save("data/r_RK4.dat");
      v.save("data/v_RK4.dat");
    }

  return 0;
}
