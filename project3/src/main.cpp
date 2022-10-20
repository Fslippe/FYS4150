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
  arma::vec r0 = arma::vec({20, 0., 20});
  arma::vec v0 = arma::vec({0.,25,0.});
  Particle p_in = Particle(q_ca, m_ca, r0, v0);
  arma::vec r0_2 = arma::vec({10, 0., 10});
  arma::vec v0_2 = arma::vec({0.,25,5.});
  Particle p_in_2 = Particle(q_ca, m_ca, r0_2, v0_2);
  PenningTrap pt = PenningTrap(B0, V0, d);
  pt.add_particle(p_in);
  pt.add_particle(p_in_2);

  std::cout << pt.p[0].info();
  std::cout << pt.p[1].info();

  pt.evolve_forward_Euler(5.);
  std::cout << pt.p[0].info();
  std::cout << pt.p[1].info();




  return 0;
}
