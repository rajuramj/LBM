#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "Lattice.hpp"
#include <memory>   //for shared pointer

class Simulation{

private:
    //src and dest points to a lattice object, from/to where the information is to be read and written resp;
    std::shared_ptr<Lattice> src;   // Similar to Lattice *ptr
    std::shared_ptr<Lattice> dest;

public:
    Simulation(const size_t&, const size_t&);

    // Prints lattice contents
    void printLattice();

};

#endif
