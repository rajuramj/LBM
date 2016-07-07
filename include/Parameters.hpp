#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

#include "Type.hpp"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>     /* atoi */
#include <string>     // std::string, std::stoi


using namespace std;

extern real nx, ny;
extern real latticeVisc, latticeAcc;
extern real relaxRate; //relaxation rate

class Parameters
{
public:
    Parameters(std::string&);
    void calcDomDim();      // To calculate the dimensions of the Domain
//    int getScenario(const Scenario&);

//    int scene;
//    std::string scene_;
//    Scenario& scenario;

    int scene_;

    inline real convVisc(real visc, real dx, real dt);
    inline real convAcc(real acc, real dx, real dt);

private:

    real length_, width_, dia_, centerX_, centerY_;
    real viscosity, simTime, acceleration;
    size_t cylinderResolution;
    real dx, dt;        // cell width and Timestep

};


inline real Parameters::convVisc(real visc, real dx, real dt)
{
    real dx_inv = 1.0 / (dx*dx) ;
    latticeVisc = visc * dt * dx_inv;
    return latticeVisc;
}

inline real Parameters::convAcc(real acc, real dx, real dt)
{
    latticeAcc = (acc * dt * dt) /dx;
    return latticeAcc;
}




#endif // SIMULATIONPARAMETERS_HPP
