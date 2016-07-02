#include "Lattice.hpp"

//Initialise the static map
// C, N W, N, N E, W, E, SW, S, SE

Lattice::dirMap.insert(Pair ("C", 0));
Lattice::dirMap.insert(Pair ("N", 1));
Lattice::dirMap.insert(Pair ("S", 2));
Lattice::dirMap.insert(Pair ("W", 3));
Lattice::dirMap.insert(Pair ("E", 4));
Lattice::dirMap.insert(Pair ("NE", 5));
Lattice::dirMap.insert(Pair ("NW", 6));
Lattice::dirMap.insert(Pair ("SW", 7));
Lattice::dirMap.insert(Pair ("SE", 8));


Lattice::Lattice(const size_t& dim_x, const size_t& dim_y){

    std::cout << "c'tr of Lattice" << std::endl;

    this->numCellsX = dim_x;
    this->numCellsY = dim_y;
    this->data_.resize(dim_x * dim_y * NUM_DIR);
}

//real& Lattice::operator() (const size_t& i, const size_t& j, const std::string& dir)
real& Lattice::operator() (const size_t& i, const size_t& j) {

    assert(i <numCellsX &&  j <numCellsY);
    return NUM_DIR*(j*numCellsX + i) + dirMap[dir];
}

//const real& Lattice::operator() (const size_t& i, const size_t& j, const std::string& dir)
const real& Lattice::operator() (const size_t& i, const size_t& j) const {

    assert(i <numCellsX &&  j <numCellsY);
    return NUM_DIR*(j*numCellsX + i) + dirMap[dir];
}

void Lattice::display() {

    for (auto const &f_q : this->data_)
        {
            std::cout << f_q ;
        }

    std::cout << std::endl;
}



