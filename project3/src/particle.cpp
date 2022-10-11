// Definitions for the functions in the particle class

#include "particle.hpp"

// Constructor
particle::particle(int charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in)
{
  charge_ = charge_in;
  mass_ = mass_in;
  position_ = position_in;
  velocity_ = velocity_in;
}

// Method that returns the charge
int particle::charge()
{
  return charge_;
}

// Method that returns the mass
double particle::mass()
{
  return mass_;
}

// Method that returns the position
arma::vec particle::position()
{
  return position_;
}

// Method that returns the velocity
arma::vec particle::velocity()
{
  return velocity_;
}

// Method that returns a string with info in the form "charge mass position velocity"
std::string particle::info()
{
  std::ostringstream r,v;
  r << position_;
  v << velocity_;
  std::string info_string = std::to_string(charge_) + " " + std::to_string(mass_) + " " + r.str() + " " + v.str();
  return info_string;
}