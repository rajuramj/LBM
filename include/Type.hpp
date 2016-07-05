#ifndef TYPE_HPP
#define TYPE_HPP
#include <utility> // for std::pair
#include <string>
#include <map>

// Lets rename double as real
typedef double real;
//typedef std::pair<std::string, size_t> Pair;
//typedef std::map<std::string, size_t> Map;

// Enum to define the mapping from direction to index in std::vector data_.
typedef enum  {
    C = 0,
    N = 1,
    S = 2,
    W = 3,
    E = 4,
    NE = 5,
    NW = 6,
    SW = 7,
    SE = 8
} Direction;

//// Enum to define the weights of the stencil which is used in init() function of Lattice class.
//typedef enum  {
//    w_C = 4.0/9.0,
//    w_N = 1.0/9.0,
//    w_S = w_N,
//    w_W = w_N,
//    w_E = w_N,
//    w_NE = 1.0/36.0,
//    w_NW = w_NE,
//    w_SW = w_NE,
//    w_SE = w_NE
//} Stencil;

#endif
