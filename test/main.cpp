#include "Lattice.hpp"

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

    std::cout << "Data successfully read!" << std::endl;

//    //Lattice::fillMap();
//     Lattice l1(2,2);
//    //l1(std::string("N"),0, 0) = 2.5;
//    l1(0, 0) = 1.2;
//    l1(1,0) = 2.4;

    Lattice l1(2,2);
    l1(0, 0, l1.C) = 1.5;
    l1(0, 0, l1.N) = 2.5;
    l1(0, 0, l1.NE) = 2.455;

    l1.display();

    std::cout << "l1(0,1,l1.NE) is : " << l1(0,0,l1.NE) << std::endl;

    return 0;
}

