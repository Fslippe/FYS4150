// g++ main.cpp ising_model.cpp -o -larmadillo

#include "ising_model.hpp"

int main()
{
    IsingModel IM = IsingModel(2);
    IM.init_lattice(1);
    int a = IM.test();
    std::cout << a; //.MC_sample(1);
    return 0;
}