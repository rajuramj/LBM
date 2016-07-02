#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include "Type.hpp"

/*using std::cout;
using std::cin;
using std::endl;
using std::ostream;
*/

class Vector
{
    public:
    std::vector<double> vec_;


    Vector(const size_t&, real=0);              // Constructor
    Vector();
    ~Vector();                                 // Destructor

    Vector operator+ (const Vector& );
    Vector operator- (const Vector& );
    Vector operator- ();                 // v = -v
    real operator* (const Vector& );    // Dot product of 2 vectors
    Vector operator* (const real& );    // Multiplying vector with a scalar
    void operator= (const Vector& );
    void display();
} ;

#endif
