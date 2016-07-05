#ifndef TYPE_HPP
#define TYPE_HPP
#include <utility> // for std::pair
#include <string>
#include <map>

// Lets rename double as real in the assignment
typedef double real;
typedef std::pair<std::string, size_t> Pair;
typedef std::map<std::string, size_t> Map;

//enum Direction {
//    C = 0,
//    N = 1,
//    S = 2,
//    W = 3,
//    E = 4,
//    NE = 5,
//    NW = 6,
//    SW = 7,
//    SE = 8
//};

#endif
