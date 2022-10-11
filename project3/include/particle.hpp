// particle class

#ifndef __particle_hpp__  
#define __particle_hpp__

#include <string>
#include <armadillo>

class particle
{
  private:

    int charge_;
    double mass_;
    arma::vec position_;
    arma::vec velocity_;

  public:

    // Constructor
    particle(int charge_in, double mass_in, arma::vec position_in, aram::vec velocity_in);

    // Method that returns the charge
    int charge();

    // Method that returns the mass
    double mass();

    // Method that returns the position
    double position();

    // Method that returns the velocity
    double velocity();

    // Method that returns a string with info in the form "charge mass position velocity"
    std::string info();

};

#endif