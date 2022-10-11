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
int particle::position()
{
  return position_;
}

// Method that returns the velocity
int particle::velocity()
{
  return velocity_;
}

// Method that returns a string with info in the form "charge mass position velocity"
std::string particle::info()
{
  std::string info_string = std::to_string(charge_) + " " + std::to_string(mass_) + " " + std::to_string(position_) + " " + std::to_string(velocity_);
  return info_string;
}