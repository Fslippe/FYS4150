// Definitions for the functions in the Particle class

#include "particle.hpp"

// Constructor
Particle::Particle(int charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in)
{
  q = charge_in;
  m = mass_in;
  r = position_in;
  v = velocity_in;
}

// Method that returns the charge
int Particle::charge()
{
  return q;
}

// Method that returns the mass
double Particle::mass()
{
  return m;
}

// Method that returns the position
arma::vec Particle::position()
{
  return r;
}

// Method that returns the velocity
arma::vec Particle::velocity()
{
  return v;
}

// Method that returns a string with info in the form "charge mass position velocity"
std::string Particle::info()
{
  std::ostringstream r_print,v_print;
  r_print << r;
  v_print << v;
  std::string info_string = std::to_string(q) + " " + std::to_string(m) + " " + r_print.str() + " " + v_print.str();
  return info_string;
}
