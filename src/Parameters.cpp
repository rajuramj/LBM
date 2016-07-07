#include "Parameters.hpp"


Parameters::Parameters(std::string& temp) : length_(0.06), width_(0.02), dia_(0.005), centerX_(0.02), centerY_(0.008)
{

    std::cout << "c'tr of Parameter " << std::endl;

    //   scene_ = temp;
    if (temp == "scenario1") { scene_ =1; }
    else if(temp == "scenario2") {scene_ = 2;}
    else { std::cout<< "No matching scenarios found !" <<std::endl; }

}

//int getScenario( const Scenario& scene_)
//{
////int scene = std::stoi("scenario");
//    int scene = scene_ ;
//    return scene;
//}


void Parameters::calcDomDim()
{
//    std::string scenario;
//    int i1 = std::stoi("scenario1");
//    int i2 = std::stoi("scenario2");

//    int scene = Scenario& scene_;
//    int scene = getScenario(Scenario& scene_);



    switch (scene_)
    {

    case 1: // scenario1
        viscosity = 1e-06;           // m^2/s
        simTime = 3;                // seconds
        acceleration =  0.01;       // m/s^2
        cylinderResolution = 30;    // no. of cells

        std::cout<<" \n " << std::endl;
        std::cout<< "*****Param of Scenario 1 ******"<< std::endl;

        std::cout<< "Diameter :" << dia_ << std::endl;
        std::cout<< "Cylinder Resolution :" << cylinderResolution << std::endl;

        std::cout<< "length_ :" << length_ << std::endl;
        std::cout<< "width_ :" << width_ << std::endl;
        std::cout<<" \n "<< std::endl;

        dx = dia_ / cylinderResolution;
        //dt = 1e-11;
        dt = 1e-4;

        std::cout<< "dx :" << dx << std::endl;
        std::cout<< "dt :" << dt << std::endl;

        nx = length_ / dx;
        ny = width_ / dx;

        latticeVisc = convVisc(viscosity, dx, dt);
        latticeAcc = convAcc(acceleration, dx, dt);

        relaxRate = 1 / ((3*latticeVisc) + 0.5);

        std::cout<<" \n " << std::endl;
        std::cout<< "*****Param after conversion ******"<< std::endl;

        std::cout<< "No. of cells in X :" << nx << std::endl;
        std::cout<< "No. of cells in Y :" << ny << std::endl;

        std::cout<< "latticeVisc :" << latticeVisc << std::endl;
        std::cout<< "latticeAcc :" << latticeAcc << std::endl;
        std::cout<< "relaxRate :" << relaxRate << std::endl;
        std::cout<<" \n "<< std::endl;

        break;

    case 2: //scenario2
        viscosity = 1e-06;           // m^2/s
        simTime = 5;                // seconds
        acceleration =  0.016;       // m/s^2
        cylinderResolution = 60;    // no. of cells

        std::cout<<" \n " << std::endl;
        std::cout<< "*****Param of Scenario 2 ******"<< std::endl;

        std::cout<< "Diameter :" << dia_ << std::endl;
        std::cout<< "Cylinder Resolution :" << cylinderResolution << std::endl;

        std::cout<< "length_ :" << length_ << std::endl;
        std::cout<< "width_ :" << width_ << std::endl;
        std::cout<<" \n "<< std::endl;

        dx = dia_ / cylinderResolution;
        dt = 1e-4;

        std::cout<< "dx :" << dx << std::endl;
        std::cout<< "dt :" << dt << std::endl;

        nx = length_ / dx;
        ny = width_ / dx;

        latticeVisc = convVisc(viscosity, dx, dt);
        latticeAcc = convAcc(acceleration, dx, dt);

        relaxRate = 1 / ((3*latticeVisc) + 0.5);

        std::cout<<" \n " << std::endl;
        std::cout<< "*****Param after conversion ******"<< std::endl;

        std::cout<< "No. of cells in X :" << nx << std::endl;
        std::cout<< "No. of cells in Y :" << ny << std::endl;

        std::cout<< "latticeVisc :" << latticeVisc << std::endl;
        std::cout<< "latticeAcc :" << latticeAcc << std::endl;
        std::cout<< "relaxRate :" << relaxRate << std::endl;
        std::cout<<" \n "<< std::endl;

        break;

     default:
        std::cerr << "Please choose scenario1 OR scenario2\n";

    }


}
