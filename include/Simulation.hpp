#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "Lattice.hpp"
#include <memory>       //for shared pointer

class Simulation{

private:
    //src and dest points to a lattice object, from/to where the information is to be read and written resp;
    std::shared_ptr<Lattice> src;   // Similar to Lattice *src
    std::shared_ptr<Lattice> dest;

    size_t numCellsX;  // This includes ghost cells
    size_t numCellsY;

    // Arrays to denote the  vectors in x and y directions.
    static constexpr int dir_x[] = {0, 0, 0, -1, 1, 1, -1, -1, 1};
    static constexpr int dir_y[] = {0, 1, -1, 0, 0, 1, 1, -1, -1};

public:
    Simulation(const size_t&, const size_t&);

    // Prints lattice contents
    void printLattice();

    // Sets the periodic BC's in East and West directions
    void setPeriodicBCs();

    // Sets the reflecting BC's in North and South directions
    void setNoSlipBCs();

    // Stream and collide are coded in one function to implement loop fusion
    void stream_Collide();

    // Perform all simulation steps of LBM
    void runSimulation();

};

#endif
