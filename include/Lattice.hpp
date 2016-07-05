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
    size_t numCellsX;
    size_t numCellsY;

    //vector to store probability density function values.
    std::vector<real> data_;

public:
    // enum data for directions
    enum Direction {
        C = 0,
        N = 1,
        S = 2,
        W = 3,
        E = 4,
        NE = 5,
        NW = 6,
        SW = 7,
        SE = 8
    };

    // Maps the direction to index
    static std::map<std::string, size_t> dirMap; //= {{"C", 0}, {"N", 1}, {"S", 2}, {"W", 3}, {"E", 4}, {"NE", 5}, {"NW", 6}, {"SW", 7}, {"SE", 8} };

    //static  void fillMap();
    //Constructor
    Lattice(const size_t&, const size_t&);

    //Non const version, used for assigning
    //real& operator() (const size_t&, const size_t&);
    real& operator() (const size_t&, const size_t&, const Direction&);

    //Const version, used for accessing const array object, Safe(returns const reference)
    //const real& operator() (const size_t&, const size_t&) const;
    const real& operator() (const size_t&, const size_t&, const Direction&) const;

    void display() const;
};

#endif
