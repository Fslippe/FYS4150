//g++ main.cpp penningtrap.cpp particle.cpp -o main -larmadillo

#include "particle.hpp"
#include "penningtrap.hpp"

int main()
{
  double B0 = 96.5;
  double V0 = 2.41 * std::pow(10, 6);

  double d = 500;
  double m_ca = 40.077; // atomic mass unit
  double q_ca = 1; // elementary charge
  double dt = 0.0005; //Time interval
  double T = 50; //Time
  int N = T / dt; //number of timepoints 
  int n = 1; //number of particles
  bool interaction = false;

  arma::vec t = arma::linspace(0, T, N);
  std::cout << t << "\n";
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
  //pt_RK4.add_particle(p_in_2);
  //pt_E.add_particle(p_in_2);

  double c = p_in.charge();
  std::cout << c;
  arma::cube r_e = arma::cube(n, N,  3);
  arma::cube r_rk4 = arma::cube(n, N,  3);
  arma::cube v_e = arma::cube(n, N,  3);
  arma::cube v_rk4 = arma::cube(n, N,  3);

  arma::mat r_a = arma::mat(pt_E.analytic(t));
  std::cout << "EJKFHIL";
  for (int i = 0; i < n; i++)
  {
    r_e(arma::span(i),arma::span(0),arma::span::all) = pt_E.p[i].r;
    v_e(arma::span(i),arma::span(0),arma::span::all) = pt_E.p[i].v;
    r_rk4(arma::span(i),arma::span(0),arma::span::all) = pt_RK4.p[i].r;
    v_rk4(arma::span(i),arma::span(0),arma::span::all) = pt_RK4.p[i].v;

    for (int j = 1; j < N; j++)
    {
      pt_E.evolve_forward_Euler(dt);
      r_e(arma::span(i),arma::span(j),arma::span::all) = pt_E.p[i].r;
      v_e(arma::span(i),arma::span(j),arma::span::all) = pt_E.p[i].v;
      pt_RK4.evolve_RK4(dt);
      r_rk4(arma::span(i),arma::span(j),arma::span::all) = pt_RK4.p[i].r;
      v_rk4(arma::span(i),arma::span(j),arma::span::all) = pt_RK4.p[i].v;
    }
  }
  
  r_a.save("data/r_a.dat");
  r_e.save("data/r_Euler.dat");
  v_e.save("data/v_Euler.dat");
  r_rk4.save("data/r_RK4.dat");
  v_rk4.save("data/v_RK4.dat");




  return 0;
}
