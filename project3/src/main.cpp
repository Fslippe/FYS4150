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
  Particle p = Particle(q_ca, m_ca, r0, v0);
  std::cout << d << "\n\n\n";
  PenningTrap pt = PenningTrap(B0, V0, d);
  pt.add_particle(p);
  std::cout << p.info();

  pt.evolve_forward_Euler(5.);



  return 0;
}
