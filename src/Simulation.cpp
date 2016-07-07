#include "Simulation.hpp"

constexpr int Simulation::dir_x[];
constexpr int Simulation::dir_y[];

Simulation::Simulation(const size_t& dim_x, const size_t& dim_y){

    this->numCellsX = dim_x + 2;
    this->numCellsY = dim_y + 2;

//    std::cout<<"dim_x & dim_y in Simulation class "<< std::endl;
//    std::cout << "dim_x :" << dim_x << std::endl;
//    std::cout << "dim_y :" << dim_y << std::endl;
//    std::cout<<" \n "<< std::endl;

        std::cout<<"=============  numCellsX & numCellsY in Simulation class =========  "<< std::endl;
        std::cout << "numCellsX :" << numCellsX << std::endl;
        std::cout << "numCellsY :" << numCellsY << std::endl;

    //Allocate memory for lattice object pointed by src,
    // Similar to src = new Lattice(dim_x + 2, dim_y + 2);
    this->src = std::make_shared<Lattice>(this->numCellsX, this->numCellsY);
    this->dest = std::make_shared<Lattice>(this->numCellsX, this->numCellsY);

    // Init src lattice with weights
    src->init();
}

void Simulation::printLattice(){

    std::cout << "Contents of src lattice\n";
    this->src->display();

//    std::cout << "Contents of dest lattice\n";
//    this->dest->display();
}


void Simulation::setPeriodicBCs(){

    const size_t i_left_src = numCellsX - 2;
    const size_t i_left_dest = 0U; // left ghost cell
    const size_t i_right_src = 1U;
    const size_t i_right_dest = numCellsX - 1; // right ghost cell

    // Iterate over all non ghost cells rows
    for(size_t j=1; j< numCellsY - 1; ++j){

        // Iterate over all direction inside a cell
        for(size_t q=0; q< NUM_DIR; ++q){

            //filling left ghost cell
            (*src)(i_left_dest,j,q) = (*src)(i_left_src,j,q);

            //filling right ghost cell
            (*src)(i_right_dest,j,q) = (*src)(i_right_src,j,q);
        }

    }

    std::cout << "Period BC's set successfully\n";
}

void Simulation::setNoSlipBCs(){

    //Down ghost layer
    size_t j_src = 1U;
    size_t j_dest = 0U;

    // The for loop are written seperately to utilise cache lines effectively
    for(size_t i=1; i< numCellsX - 1; ++i) {

        (*src)(i-1, j_dest, NE) =  (*src)(i, j_src, SW);
        (*src)(i, j_dest, N) =  (*src)(i, j_src, S);
        (*src)(i+1, j_dest, NW) =  (*src)(i, j_src, SE);
    }

    //Top Ghost layer
    j_src = numCellsY - 2;
    j_dest = numCellsY - 1;

    for(size_t i=1; i< numCellsX - 1; ++i) {

        (*src)(i-1, j_dest, SE) =  (*src)(i, j_src, NW);
        (*src)(i, j_dest, S) =  (*src)(i, j_src, N);
        (*src)(i+1, j_dest, SW) =  (*src)(i, j_src, NE);
    }

    std::cout << "Reflective BC's set successfully\n";

}



