// To be added to PenningTrap class

// Analytic solution for time evolution of particle motion
// Takes in a vektor t with time steps and outputs three vektors (x, y, z)
std::vec PenningTrap::analytic(std::vec t)
{
  if (p.size() == 1)
  {
    double omega_0 = p[0].charge() * B_0 / p[0].mass();
    double omega_z = std::sqrt( (2 * p[0].charge() * V_0) / (p[0].mass() * d*d) );
    double omega_plus = 0.5 * (omega_0 + std::sqrt(omega_0 * omega_0 - 2* omega_z * omega_z));
    double omega_minus =  0.5 * (omega_0 - std::sqrt(omega_0 * omega_0 - 2* omega_z * omega_z));
  
    double x_0 = p[0].position()[0];
    double y_0 = p[0].position()[1];
    double z_0 = p[0].position()[2];

    // The x and z component of v are 0. only using y component.
    double v_0 = p[0].velocity()[1];
   

    double A_plus = (v_0 + omega_minus * x_0) / (omega_minus - omega_plus);
    double A_minus = (v_0 + omega_plus * x_0) / (omega_minus - omega_plus);

    std::vec x = A_plus * std::cos(omega_plus *t) + A_minus * std::cos(omega_minus *t);
    std::vec y = - A_plus * std::sin(omega_plus *t) - A_minus * std::sin(omega_minus *t);
    std::vec z = z_0 * std::cos(omega_z * t);
  }
  else
  {
    std::cout << 'Too many particles! \n The analytic solution is only designed for a singel particle.' << "\n";
  }
  
    return std::vec(x, y, z)
}