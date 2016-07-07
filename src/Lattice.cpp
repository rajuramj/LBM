#include "Lattice.hpp"

//Initialise the static map
// C, N W, N, N E, W, E, SW, S, SE

//Map Lattice::dirMap ={{"C", 0}, {"N", 1}, {"S", 2}, {"W", 3}, {"E", 4}, {"NE", 5}, {"NW", 6}, {"SW", 7}, {"SE", 8} };


constexpr real Lattice::weights[] ;//= {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

Lattice::Lattice(const size_t& dim_x, const size_t& dim_y){

    std::cout << "c'tr of Lattice" << std::endl;
    std::cout<<" \n "<< std::endl;

    this->numCellsX = dim_x;
    this->numCellsY = dim_y;
    this->data_.resize(dim_x * dim_y * NUM_DIR);
}


// Initialise the lattice with weights.
void Lattice::init() {

    for(size_t i=0; i< data_.size(); ++i) {
    this->data_[i] = weights[i % NUM_DIR] + i/NUM_DIR;
    }
}



real& Lattice::operator() (const size_t& i, const size_t& j, const Direction& dir){

   assert(i>=0 && j>=0  && i <numCellsX &&  j <numCellsY);
    //std::cout << "NON const version operator() called\n";
    return this->data_[NUM_DIR*(j*numCellsX + i) + dir];
}

real& Lattice::operator() (const size_t& i, const size_t& j, const size_t& k){

    assert(i>=0 && j>=0 && k>=0 && i <numCellsX &&  j <numCellsY && k <NUM_DIR);
    //std::cout << "NON const version operator() called\n";
    return this->data_[NUM_DIR*(j*numCellsX + i) + k];
}

const real& Lattice::operator() (const size_t& i, const size_t& j, const Direction& dir) const{

    assert(i>=0 && j>=0 && i <numCellsX &&  j <numCellsY);
    //std::cout << "Const version operator() called\n";

    return this->data_[NUM_DIR*(j*numCellsX + i)+ dir];
}

const real& Lattice::operator() (const size_t& i, const size_t& j, const size_t& k) const{

    assert(i>=0 && j>=0 && k>=0 && i <numCellsX &&  j <numCellsY && k <NUM_DIR);
    std::cout << "Const version operator() called\n";

    return this->data_[NUM_DIR*(j*numCellsX + i)+ k];
}

void Lattice::display() const {

//    size_t counter =0;
//    for (auto const &f_q : this->data_)
//        {
//            std::cout << f_q << "\t" ;

//            counter = (counter + 1) % NUM_DIR;
//            if(counter ==0)
//                std::cout << std::endl;
//        }

    for(int j=numCellsY - 1; j >=0; --j){
    //size_t j=0;
        std::cout << std::endl;

        for(size_t i = 0; i < numCellsX; ++i){
            std::cout << (*this)(i,j, NW) << "\t";
            std::cout << (*this)(i,j, N) << "\t";
            std::cout << (*this)(i,j, NE) << "\t";
            std::cout << "\t";
        }
        std::cout << std::endl;

        for(size_t i = 0; i < numCellsX; ++i){
            std::cout << (*this)(i,j, W) << "\t";
            std::cout << (*this)(i,j, C) << "\t";
            std::cout << (*this)(i,j, E) << "\t";
            std::cout << "\t";
        }
        std::cout << std::endl;

        for(size_t i = 0; i < numCellsX; ++i){
            std::cout << (*this)(i,j, SW) << "\t";
            std::cout << (*this)(i,j, S) << "\t";
            std::cout << (*this)(i,j, SE) << "\t";
            std::cout << "\t";
        }
        std::cout << std::endl;
    }


}



