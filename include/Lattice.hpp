#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Type.hpp"
#include <vector>
#include <map>
#include <iostream>
#include <cassert>

#define NUM_DIR 9

class Lattice{

private:
    size_t numCellsX;  // This includes ghost cells
    size_t numCellsY;

    //vector to store probability density function(f_q) values.
    std::vector<real> data_;

    // init lattice weights
    static constexpr real weights[] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

public:
    //Constructor
    Lattice(const size_t&, const size_t&);

    // Sets the initails f_q's
    void init();

    //Non const version, used for assigning
     real& operator() (const size_t&, const size_t&, const Direction&);
     real& operator() (const size_t&, const size_t&, const size_t&);

    //Const version, used for accessing const array object, Safe(returns const reference)
    // This operator is never used in this assignment
    const real& operator() (const size_t&, const size_t&, const Direction&) const;
    const real& operator() (const size_t&, const size_t&, const size_t&) const;

    void display() const;
};



#endif
