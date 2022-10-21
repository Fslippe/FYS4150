//g++ main.cpp penningtrap.cpp particle.cpp -o main -larmadillo

#include "particle.hpp"
#include "penningtrap.hpp"

int main()
{
  double B0 = 96.5;
  double V0 = 2.41 * std::pow(10, 7);

  double d = 500;
  double m_ca = 40.077; //atomic mass unit
  double q_ca = 1; //elementary charge
  double dt = 0.0005; 
  int T = 1;
  int N = T / dt;
  int n = 2; //number of particles
  PenningTrap pt_RK4 = PenningTrap(B0, V0, d);
  PenningTrap pt_E = PenningTrap(B0, V0, d);


  arma::vec r0 = arma::vec({20, 0., 20});
  arma::vec v0 = arma::vec({0.,25,0.});
  Particle p_in = Particle(q_ca, m_ca, r0, v0);
  pt_RK4.add_particle(p_in);
  pt_E.add_particle(p_in);


  arma::vec r0_2 = arma::vec({10, 0., 10});
  arma::vec v0_2 = arma::vec({0.,25,5.});
  Particle p_in_2 = Particle(q_ca, m_ca, r0_2, v0_2);
  pt_RK4.add_particle(p_in_2);
  pt_E.add_particle(p_in_2);

  arma::cube r_e = arma::cube(n, N,  3);
  arma::cube r_rk4 = arma::cube(n, N,  3);
  arma::cube v_e = arma::cube(n, N,  3);
  arma::cube v_rk4 = arma::cube(n, N,  3);
  std::cout << r_e.tube(0,2);

  for (int j = 0; j < n; j++)
  {
    for (int i = 0; i < N; i++)
  {
    std::cout << i << j << "\n";
    pt_E.evolve_forward_Euler(dt);
    r_e.subcube(j, i, 0, j, i, 2) = pt_E.p[i].r;
    v_e.subcube(j, i, 0, j, i, 2) = pt_E.p[i].v;
    pt_RK4.evolve_RK4(dt);
    r_rk4.subcube(j, i, 0, j, i, 2) = pt_RK4.p[i].r;
    v_rk4.tube(j, i) = pt_RK4.p[i].v;
  }
  }
  
  
  
  r_e.save("data/r_Euler.dat");
  v_e.save("data/v_Euler.dat");
  r_rk4.save("data/r_RK4.dat");
  v_rk4.save("data/v_RK4.dat");




  return 0;
}
