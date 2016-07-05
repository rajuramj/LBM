#include "Lattice.hpp"

//Initialise the static map
// C, N W, N, N E, W, E, SW, S, SE

//Map Lattice::dirMap ={{"C", 0}, {"N", 1}, {"S", 2}, {"W", 3}, {"E", 4}, {"NE", 5}, {"NW", 6}, {"SW", 7}, {"SE", 8} };


constexpr real Lattice::weights[] ;//= {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

Lattice::Lattice(const size_t& dim_x, const size_t& dim_y){

    std::cout << "c'tr of Lattice" << std::endl;

    this->numCellsX = dim_x;
    this->numCellsY = dim_y;
    this->data_.resize(dim_x * dim_y * NUM_DIR);
}


// Initialise the lattice with weights.
void Lattice::init() {

    for(size_t i=0; i< data_.size(); ++i) {
    this->data_[i] = weights[i % NUM_DIR];
    }
}




real& Lattice::operator() (const size_t& i, const size_t& j, const Direction& dir){

    assert(i <numCellsX &&  j <numCellsY);
    std::cout << "NON const version operator() called\n";
    return this->data_[NUM_DIR*(j*numCellsX + i) + dir];
}

const real& Lattice::operator() (const size_t& i, const size_t& j, const Direction& dir) const{
    assert(i <numCellsX &&  j <numCellsY);
    std::cout << "Const version operator() called\n";

    return this->data_[NUM_DIR*(j*numCellsX + i)+ dir];
}

void Lattice::display() const {

    for (auto const &f_q : this->data_)
        {
            std::cout << f_q << "\t" ;
        }

    std::cout << std::endl;
}



