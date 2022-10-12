// PenningTrap class

#ifndef __penningtrap_hpp__  
#define __tenningtrap_hpp__

#include <string>
#include <armadillo>

class PenningTrap
{
  public:
    // Member variables
    double B_0; // magnetic field strength
    double V_0; // applied potetial
    double d; // characteristic dimension
    std::vector<Particle> particle_objects; // conntains Particle objects in the Penning trap

    // Constructor
    PenningTrap(double B_0_in, double V_0_in, double d_in, std::vector<Particle> particle_objects_in);


    //Method that adds a particle to the Penning trap
    void add_p(Particle p_in);

    // Method that returns the ecternal electric field
    arma::vec ext_el_field();

    // Method that returns the ecternal magnetic field
    arma::vec ext_mag_field();

    // Method that returns the force due to the interaction among particles
    arma::vec p_interaction_force();

};

#endif