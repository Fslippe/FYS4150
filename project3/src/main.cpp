//g++ main.cpp penningtrap.cpp particle.cpp -o main -larmadillo

#include "particle.hpp"
#include "penningtrap.hpp"
#include <stdlib.h>


int main(int argc, char** argv)
{
  if (argc != 5)
  {
    std::cout << "Expected 5 commandline arguments: " <<
    "int N (Timepoints) double T (Time) int n (number of particles)" 
    << "bool interaction (particle interaction true or false)]\n\nExiting...\n";
    exit(1);
  }

  int N = atoi(argv[1]); //number of timepoints 
  double T = atoi(argv[2]); //Time
  int n = atoi(argv[3]); //number of particles
  bool interaction;
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
  arma::vec t = arma::linspace(0, T, N);

  PenningTrap pt_RK4 = PenningTrap(B0, V0, d, interaction);
  PenningTrap pt_E = PenningTrap(B0, V0, d, interaction);

  arma::vec r0 = arma::vec({20, 0., 20});
  arma::vec v0 = arma::vec({0.,25,0.});
  Particle p_in = Particle(q_ca, m_ca, r0, v0);
  pt_RK4.add_particle(p_in);
  pt_E.add_particle(p_in);

  arma::vec r0_2 = arma::vec({25, 25., 0});
  arma::vec v0_2 = arma::vec({0.,40, 5});
  Particle p_in_2 = Particle(q_ca, m_ca, r0_2, v0_2);
  pt_RK4.add_particle(p_in_2);
  pt_E.add_particle(p_in_2);

  arma::cube r_e = arma::cube(n, N,  3);
  arma::cube r_rk4 = arma::cube(n, N,  3);
  arma::cube v_e = arma::cube(n, N,  3);
  arma::cube v_rk4 = arma::cube(n, N,  3);

  //arma::mat r_a = arma::mat(pt_E.analytic(t));


  for (int i = 0; i < n; i++)
  {
    r_e(arma::span(i),arma::span(0),arma::span::all) = pt_E.p[i].r;
    v_e(arma::span(i),arma::span(0),arma::span::all) = pt_E.p[i].v;
    r_rk4(arma::span(i),arma::span(0),arma::span::all) = pt_RK4.p[i].r;
    v_rk4(arma::span(i),arma::span(0),arma::span::all) = pt_RK4.p[i].v;
  }
  
  for (int j = 1; j < N; j++)
  {
    pt_RK4.evolve_RK4(dt);
    pt_E.evolve_forward_Euler(dt);
    for (int i = 0; i < n; i++)
    {
      r_e(arma::span(i),arma::span(j),arma::span::all) = pt_E.p[i].r;
      v_e(arma::span(i),arma::span(j),arma::span::all) = pt_E.p[i].v;
      r_rk4(arma::span(i),arma::span(j),arma::span::all) = pt_RK4.p[i].r;
      v_rk4(arma::span(i),arma::span(j),arma::span::all) = pt_RK4.p[i].v;
    }
  }
  
  
  //r_a.save("data/r_a.dat");
  r_e.save("data/r_Euler.dat");
  v_e.save("data/v_Euler.dat");
  r_rk4.save("data/r_RK4.dat");
  v_rk4.save("data/v_RK4.dat");




  return 0;
}
