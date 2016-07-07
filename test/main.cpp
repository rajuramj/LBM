#include "Lattice.hpp"
#include "Parameters.hpp"
#include "Simulation.hpp"

real nx, ny;
real latticeVisc, latticeAcc;
real relaxRate; //relaxation rate


int main(int argc, char** argv)
{	
    if(argc < 2)
        std::cerr<<"Insufficient number of input parameters"<<std::endl;

    // The program will terminate after printing error message.
    assert(argc ==2);

    std::string s1("scenario1");
    std::string s2("scenario2");

    if(s1.compare(argv[1])  && s2.compare(argv[1])) {

        std::cerr << "argv[1] must be scenario1 or scenario2!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout<<" \n "<< std::endl;
    std::cout << "Data successfully read!" << std::endl;

//    Lattice l1(2,2);
//    l1(0, 0, C) = 1.5;
//    l1(0, 0, N) = 2.5;
//    l1(0, 0, NE) = 2.455;

//    l1.display();

//    std::cout << "l1(0,1,NE) is : " << l1(0,0,NE) << std::endl;

    std::string arg = argv[1];
    Parameters param(arg);
    param.calcDomDim();
    std::cout<< "Param Converted !"<< std::endl;

    Simulation sim(3,3);
    sim.printLattice();
    sim.setNoSlipBCs();
    sim.printLattice();

    return 0;
}

