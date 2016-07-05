#include "Simulation.hpp"

Simulation::Simulation(const size_t& dim_x, const size_t& dim_y){

    //Allocate memory for lattice object pointed by src,
    // Similar to src = new Lattice(dim_x + 2, dim_y + 2);
    this->src = std::make_shared<Lattice>(dim_x + 2, dim_y + 2);
    this->dest = std::make_shared<Lattice>(dim_x + 2, dim_y + 2);

    // Init src lattice with weights
    src->init();
}

void Simulation::printLattice(){

    std::cout << "Contents of src lattice\n";
    this->src->display();

    std::cout << "Contents of dest lattice\n";
    this->dest->display();
}
