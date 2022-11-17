#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

//Decalarations
// Expectation values
double E(double beta, int J ,double Z);

double E_pow2(double beta, int J, double Z);

// print this val
double e(double beta, int J, double Z);

double M(double beta, int J, double Z);

double M_pow2(double beta, int J, double Z);

//print this val
double m(double beta, int J, double Z);
// print this val
double C_v(double beta, int J, double Z, int T);

// print this val
double X(double beta, int J, double Z, int T);

int main()
{
    int J = 1;
    int T = 1*J;

    // inverse temperature
    double beta = 1 / T;

    //partitition function
    double Z = 2*exp(beta*8*J) + 2*exp(-beta*8*J) + 12;
    //double Z = 4*(cosh(9*beta*J) + 3);
    cout << "Z: "<<Z;

    double expect_e = e(beta, J, Z);
    double expect_m = m(beta, J, Z);
    double expect_C_v = C_v(beta, J, Z, T);
    double expect_X = X(beta, J, Z, T);
    cout << "\ne: " << expect_e;
    cout << "\nm: " << expect_m;
    cout << "\nC_v: " << expect_C_v;
    cout << "\nX: " << expect_X;

    
    return 0;
}

// Expectation values
double E(double beta, int J ,double Z)
{
    return (J / Z) * (16*exp(- beta * 8 * J) - 16*exp(beta * 8 * J));
}

double E_pow2(double beta, int J, double Z)
{
    return (J*J / Z) * (128*exp(- beta * 8 * J) + 128*exp(beta * 8 * J));
}

// print this val
double e(double beta, int J, double Z)
{
    return E(beta, J, Z) / 4;
}

double M(double beta, int J, double Z)
{
    return (8*exp(beta * 8 * J) + 16)/ Z;
}

double M_pow2(double beta, int J, double Z)
{
    return (32*exp(beta * 8 * J) + 32)/ Z;
}

//print this val
double m(double beta, int J, double Z)
{
    return M(beta, J, Z) / 4;
}
// print this val
double C_v(double beta, int J, double Z, int T)
{
    return (E_pow2(beta , J, Z) - ((E(beta, J, Z)*E(beta, J, Z)))) / (4 * T);
}

// print this val
double X(double beta, int J, double Z, int T)
{
    return (M_pow2(beta , J, Z) - (M(beta, J, Z)*M(beta, J, Z))) / (4 * T);
    // where does the last k_b go?
}

