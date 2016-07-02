/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
 *  Address: EPFL-STI-LIN Station 9, 1015 Lausanne, Switzerland
 *  E-mail: orestis.malaspinas@epfl.ch
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#include "olb2D.h"
#include "olb2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

#define NSDESCRIPTOR ForcedD2Q9Descriptor
#define TDESCRIPTOR AdvectionDiffusionD2Q5Descriptor

const int maxIter  = 1000000000;
const int saveIter = 2500;

void iniGeometry( BlockStructure2D<T, NSDESCRIPTOR>& NSlattice,
                  BlockStructure2D<T, TDESCRIPTOR>& ADlattice,
                  ForcedBGKdynamics<T, NSDESCRIPTOR> &bulkDynamics,
                  Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                  OnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>& NSboundaryCondition,
                  OnLatticeAdvectionDiffusionBoundaryCondition2D<T,TDESCRIPTOR>& TboundaryCondition,
                  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter
                )
{
  OstreamManager cout(std::cout,"iniGeometry");

  typedef advectionDiffusionLbHelpers<T,TDESCRIPTOR> TlbH;
  int nx = ADlattice.getNx();
  int ny = ADlattice.getNy();

  double Tomega  = converter.getOmegaT();
  double NSomega = converter.getOmegaNS();
  cout << "defining dynamics" << endl;

  ADlattice.defineDynamics(0,nx-1, 0,ny-1, &advectionDiffusionBulkDynamics);
  NSlattice.defineDynamics(0,nx-1, 0,ny-1, &bulkDynamics);

  NSboundaryCondition.addVelocityBoundary1P(0,nx-1,ny-1,ny-1, NSomega);
  NSboundaryCondition.addVelocityBoundary1N(0,nx-1,   0,   0, NSomega);

  TboundaryCondition.addTemperatureBoundary1P(0,nx-1,ny-1,ny-1, Tomega);
  TboundaryCondition.addTemperatureBoundary1N(0,nx-1,   0,   0, Tomega);

  for (int iX=0; iX<nx; ++iX) {
    for (int iY=0; iY<ny; ++iY) {
      T u[2] = {0.,0.};
      T rho = (T)1;
      NSlattice.get(iX,iY).defineRhoU(rho, u);
      NSlattice.get(iX,iY).iniEquilibrium(rho, u);
      T force[2] = {T(), T()};
      NSlattice.get(iX,iY).defineExternalField (
        NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
        NSDESCRIPTOR<T>::ExternalField::sizeOfForce,
        force );
    }
  }

  T Tcold = converter.getTcold();
  T Thot  = converter.getThot();

  for (int iX=0; iX<nx; ++iX) {
    for (int iY=0; iY<ny; ++iY) {
      T u[2] = {0.,0.};
      T temperature =Tcold;
      if (iY==0)
        temperature = Thot;

      ADlattice.get(iX,iY).defineRho(temperature);
      T tEq[TDESCRIPTOR<T>::q];
      for (int iPop = 0; iPop < TDESCRIPTOR<T>::q; ++iPop) {
        tEq[iPop] = TlbH::equilibrium(iPop,temperature,u);
      }
      ADlattice.get(iX,iY).definePopulations(tEq);
      ADlattice.get(iX,iY).defineExternalField (
        TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
        TDESCRIPTOR<T>::ExternalField::sizeOfVelocity,
        u );
    }
  }

  ADlattice.get(nx/2,1).defineRho(Thot+Thot/(T)5);

  NSlattice.initialize();
  ADlattice.initialize();
}

void couplingBetweenNavierStokesAndAdvectionDiffusion(
  BlockStructure2D<T, TDESCRIPTOR>& ADlattice,
  BlockStructure2D<T,NSDESCRIPTOR>& NSlattice,
  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter)
{

  int nx = NSlattice.getNx();
  int ny = NSlattice.getNy();

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<SpatiallyExtendedObject2D* > partnerForNavierStokes;

  partnerForNavierStokes.push_back(&ADlattice);
  partnerForNavierStokes.push_back(&NSlattice);

  std::vector<T> dir;
  dir.push_back(T());
  dir.push_back((T)1);

  NavierStokesAdvectionDiffusionCouplingGenerator2D<T,NSDESCRIPTOR>
  coupling(0,nx-1,0,ny-1, converter.getGravity(),
           converter.getT0(), converter.getDeltaTemperature(), dir);

  NSlattice.addLatticeCoupling(coupling,partnerForNavierStokes);
}

void plotStatistics(BlockStructure2D<T, NSDESCRIPTOR>&    NSlattice,
                    BlockStructure2D<T, TDESCRIPTOR>&    ADlattice, int iT)
{
  OstreamManager cout(std::cout,"plotStatistics");

  cout << "Writing Gif..." << std::endl;
  ImageWriter<T> imageCreator("leeloo.map");

  imageCreator.writeScaledGif(createFileName("t", iT, 7),
                              ADlattice.getDataAnalysis().getPressure(),
                              600,300);
  imageCreator.writeScaledGif(createFileName("uy", iT, 6),
                              NSlattice.getDataAnalysis().getVelocityNorm(),
                              600,300);
}

T computeNusselt(BlockStructure2D<T, NSDESCRIPTOR>& NSlattice,
                 BlockStructure2D<T,TDESCRIPTOR> &ADlattice,
                 AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter)
{
  int nx = NSlattice.getNx();
  int ny = NSlattice.getNy();

  T u_T = T();
  for (int iX = 0; iX < nx; ++iX)
  {
    for (int iY = 0; iY < ny; ++iY)
    {
      const T uy =
        NSlattice.getDataAnalysis().getVelocity().extractComponent(1).get(iX,iY);
      u_T += uy * ADlattice.get(iX,iY).computeRho();
    }
  }
  T nusselt = (T)1 + u_T*converter.getDeltaX() /
              (converter.getKappa() *
               (converter.getThot()-converter.getTcold()));

  return nusselt;
}

int main(int argc, char *argv[])
{
  OstreamManager clout(std::cout,"main");

  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  if (argc != 5)
  {
    clout << "Error : Wrong parameters specified." << endl;
    clout << "1 : Rayleigh." << endl;
    clout << "2 : Prandtl." << endl;
    clout << "3 : N." << endl;
    clout << "4 : Delta t." << endl;
    exit(1);
  }

  const double Ra    = atof(argv[1]);
  const double Pr    = atof(argv[2]);
  const int N        = atoi(argv[3]);
  const double dt    = atof(argv[4]);


  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> converter(
    Ra,   // Ra
    Pr,   // Pr
    0.0,  // Tcold
    1.0,  // Thot
    N,    // N
    dt,   // dt
    2.0,  // lx
    1.0   // ly
  );
  writeLogFile<T,NSDESCRIPTOR,TDESCRIPTOR>(converter,"2D rayleigh-benard");

  const double Raprova = converter.getN() * converter.getN() *
                         converter.getN() * converter.getDeltaTemperature() *
                         converter.getGravity() / (converter.getNu()*converter.getKappa());

  const double Prprova = converter.getNu() / converter.getKappa();

  clout << Raprova << " " << Prprova << endl;
  clout << converter.getOmegaNS() << " " << converter.getOmegaT() << endl;

  int nx = converter.getNx();
  int ny = converter.getNy();


#ifndef PARALLEL_MODE_MPI  // sequential program execution
  MultiDataDistribution2D distribution = createRegularDataDistribution( nx, ny, 1, 1, 1 );
#else                      // parallel program execution
  MultiDataDistribution2D distribution = createRegularDataDistribution( nx, ny );
#endif
  MultiBlockLattice2D<T, TDESCRIPTOR> ADlattice(distribution);
  MultiBlockLattice2D<T, NSDESCRIPTOR> NSlattice(distribution);

  OnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>* NSboundaryCondition =
    createLocalBoundaryCondition2D(NSlattice);

  OnLatticeAdvectionDiffusionBoundaryCondition2D<T,TDESCRIPTOR>* TboundaryCondition =
    createAdvectionDiffusionBoundaryCondition2D(ADlattice);

  ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getOmegaNS(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>());


  AdvectionDiffusionRLBdynamics<T, TDESCRIPTOR> TbulkDynamics (
    converter.getOmegaT(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>()
  );

  couplingBetweenNavierStokesAndAdvectionDiffusion(ADlattice, NSlattice, converter);

  iniGeometry(NSlattice, ADlattice,
              NSbulkDynamics, TbulkDynamics,
              *NSboundaryCondition, *TboundaryCondition,
              converter );

  util::ValueTracer<T> converge(0.01,(T)ny,1.0e-5);
  int iT = 0;
  clout << "starting simulation" << endl;
  for (iT=0; iT<maxIter; ++iT)
  {

    if (converge.hasConverged())
    {
      clout << "FINAL CONVERGENCE!!! " << endl;
      break;
    }

    if (iT%saveIter==0)
    {
      clout << "Writing stats at time " << iT << "." << std::endl;
      clout << ADlattice.getStatistics().getAverageEnergy() << endl;
      plotStatistics(NSlattice, ADlattice, iT );
    }

    ADlattice.collideAndStream(true);
    NSlattice.collideAndStream(true);

    converge.takeValue(ADlattice.getStatistics().getAverageEnergy(),true);

    NSlattice.executeCoupling();
    ADlattice.executeCoupling();
  }

  clout << "Time " << iT << "." << std::endl;
  clout << "Nusselt = " << computeNusselt(NSlattice, ADlattice, converter) << endl;

  delete NSboundaryCondition;
  delete TboundaryCondition;
}
