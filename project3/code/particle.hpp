// Particle class

#ifndef __particle_hpp__  
#define __particle_hpp__

#include <string>
#include <armadillo>

class Particle
{
  public:
    // Member variables
    int q;
    double m;
    arma::vec r;
    arma::vec v;

    // Constructor
    Particle(int charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in);

    // Method that returns the charge
    int charge();

    // Method that returns the mass
    double mass();

    // Method that returns the position
    arma::vec position();

    // Method that returns the velocity
    arma::vec velocity();

    // Method that returns a string with info in the form "charge mass position velocity"
    std::string info();

};

#endif