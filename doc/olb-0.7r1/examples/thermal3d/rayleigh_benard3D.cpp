/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2008 Jonas Latt, Orestis Malaspina, Andrea Parmigiani
 *  E-mail: Jonas.Latt@cui.unige.ch
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

#include "olb3D.h"
#include "olb3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#define TDESCRIPTOR AdvectionDiffusionD3Q7Descriptor
#define NSDESCRIPTOR ForcedD3Q19Descriptor

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

const int maxIter  = 1000000;
const int saveIter = 2000;

void iniGeometry( BlockStructure3D<T, NSDESCRIPTOR>& NSlattice,
                  BlockStructure3D<T, TDESCRIPTOR>& ADlattice,
                  ForcedBGKdynamics<T, NSDESCRIPTOR> &bulkDynamics,
                  Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                  OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>& NSboundaryCondition,
                  OnLatticeAdvectionDiffusionBoundaryCondition3D<T,TDESCRIPTOR>& TboundaryCondition,
                  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter
                )
{
  OstreamManager cout(std::cout,"iniGeometry");
  typedef advectionDiffusionLbHelpers<T,TDESCRIPTOR> TlbH;

  int nx = ADlattice.getNx();
  int ny = ADlattice.getNy();
  int nz = ADlattice.getNz();

  double Tomega  = converter.getOmegaT();
  double NSomega = converter.getOmegaNS();

  cout << "defining dynamics" << endl;

  ADlattice.defineDynamics(0,nx-1, 0,ny-1, 0,nz-1, &advectionDiffusionBulkDynamics);

  NSlattice.defineDynamics(0,nx-1, 0,ny-1, 0,nz-1, &bulkDynamics);

  NSboundaryCondition.addVelocityBoundary2P(0,nx-1,0,ny-1,nz-1,nz-1, NSomega);
  NSboundaryCondition.addVelocityBoundary2N(0,nx-1,0,ny-1,0,   0,    NSomega);


  TboundaryCondition.addTemperatureBoundary2P(0,nx-1,0,ny-1,nz-1,nz-1, Tomega);
  TboundaryCondition.addTemperatureBoundary2N(0,nx-1,0,ny-1,0,   0,    Tomega);

  for (int iX=0; iX<nx; ++iX)
  {
    for (int iY=0; iY<ny; ++iY)
    {
      for (int iZ=0; iZ<nz; ++iZ)
      {
        T u[3] = {0.,0.,0.};
        T rho = (T)1;
        NSlattice.get(iX,iY,iZ).defineRhoU(rho, u);
        NSlattice.get(iX,iY,iZ).iniEquilibrium(rho, u);
        T force[3] = {T(), T(), T()};
        NSlattice.get(iX,iY,iZ).defineExternalField (
          NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
          NSDESCRIPTOR<T>::ExternalField::sizeOfForce,
          force );
      }
    }
  }

  T Tcold = converter.getTcold();
  T Thot  = converter.getThot();

  for (int iX=0; iX<nx; ++iX)
  {
    for (int iY=0; iY<ny; ++iY)
    {
      for (int iZ=0; iZ<nz; ++iZ)
      {
        T u[3] = {0.,0.,0.};
        T temperature = Tcold;
        if (iZ == 0)
          temperature = Thot;


        ADlattice.get(iX,iY,iZ).defineRho(temperature);
        T tEq[TDESCRIPTOR<T>::q];
        for (int iPop = 0; iPop < TDESCRIPTOR<T>::q; ++iPop) {
          tEq[iPop] = TlbH::equilibrium(iPop,temperature,u);
        }
        ADlattice.get(iX,iY,iZ).definePopulations(tEq);
        ADlattice.get(iX,iY,iZ).defineExternalField (
          TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
          TDESCRIPTOR<T>::ExternalField::sizeOfVelocity,
          u );
      }
    }
  }

  T u[3] = {0.,0.,0.};
  T temperature = Thot+Thot/(T)5;
  ADlattice.get(nx/2,ny/2,1).defineRho(temperature);
  T tEq[TDESCRIPTOR<T>::q];
  for (int iPop = 0; iPop < TDESCRIPTOR<T>::q; ++iPop)
  {
    tEq[iPop] = TlbH::equilibrium(iPop,temperature,u);
  }
  ADlattice.get(nx/2,ny/2,1).definePopulations(tEq);
  ADlattice.get(nx/2,ny/2,1).defineExternalField (
    TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
    TDESCRIPTOR<T>::ExternalField::sizeOfVelocity,
    u );

  NSlattice.initialize();
  ADlattice.initialize();
}

void couplingBetweenNavierStokesAndAdvectionDiffusion(
  BlockStructure3D<T, TDESCRIPTOR>& ADlattice,
  BlockStructure3D<T,NSDESCRIPTOR>& NSlattice,
  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter)
{

  int nx = NSlattice.getNx();
  int ny = NSlattice.getNy();
  int nz = NSlattice.getNz();

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<SpatiallyExtendedObject3D*> partnerForNavierStokes;

  partnerForNavierStokes.push_back(&ADlattice);

  std::vector<T> dir;  // the vector of the direction of the force
  dir.push_back(T());
  dir.push_back(T());
  dir.push_back((T)1);

  NavierStokesAdvectionDiffusionCouplingGenerator3D<T,NSDESCRIPTOR>
  tCoupling(0,nx-1,0,ny-1,0,nz-1, converter.getGravity(),
            converter.getT0(), converter.getDeltaTemperature(), dir);

  NSlattice.addLatticeCoupling(tCoupling,partnerForNavierStokes);
}

T computeNusselt(BlockStructure3D<T, NSDESCRIPTOR>& NSlattice,
                 BlockStructure3D<T,TDESCRIPTOR> &ADlattice,
                 AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter)
{
  int nx = ADlattice.getNx();
  int ny = ADlattice.getNy();
  int nz = ADlattice.getNz();
  T u_T = T();
  for (int iX = 0; iX < nx; ++iX)
  {
    for (int iY = 0; iY < ny; ++iY)
    {
      for (int iZ = 0; iZ < nz; ++iZ)
      {
        const T uz = NSlattice.getDataAnalysis().getVelocity().extractComponent(2).get(iX,iY,iZ);
        u_T += uz * ADlattice.get(iX,iY,iZ).computeRho();
      }
    }
  }
  T nusselt = (T)1 + u_T /
              ((T)(nx-1)*(T)(ny-1)*converter.getKappa() *
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
  const double dt = atof(argv[4]);


  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> converter(
    Ra,  // Ra
    Pr,  // Pr
    0.0, // Tcold
    1.0, // Thot
    N,   // N
    dt,  // dt
    2.0, // lx
    1.0, // ly
    1.0  // lz
  );
  writeLogFile<T,NSDESCRIPTOR,TDESCRIPTOR>(converter,"3D rayleigh-benard");

  const double Raprova = converter.getN() * converter.getN() *
                         converter.getN() * converter.getDeltaTemperature() *
                         converter.getGravity() / (converter.getNu()*converter.getKappa());

  const double Prprova = converter.getNu() / converter.getKappa();

  clout << Raprova << " " << Prprova << endl;
  clout << converter.getOmegaNS() << " " << converter.getOmegaT() << endl;

  int nx = converter.getNx();
  int ny = converter.getNy();
  int nz = converter.getNz();

#ifndef PARALLEL_MODE_MPI  // sequential program execution
  BlockLattice3D<T, TDESCRIPTOR> ADlattice(nx, ny, nz );
  BlockLattice3D<T, NSDESCRIPTOR> NSlattice(nx, ny, nz );
#else                      // parallel program execution
  MultiDataDistribution3D distribution =
    createRegularDataDistribution( nx, ny, nz );
  MultiBlockLattice3D<T, TDESCRIPTOR> ADlattice(distribution);
  MultiBlockLattice3D<T, NSDESCRIPTOR> NSlattice(distribution);
#endif

  OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>* NSboundaryCondition =
    createLocalBoundaryCondition3D(NSlattice);

  OnLatticeAdvectionDiffusionBoundaryCondition3D<T,TDESCRIPTOR>* TboundaryCondition =
    createAdvectionDiffusionBoundaryCondition3D<T,TDESCRIPTOR,
    AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> >(ADlattice);

  ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getOmegaNS(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>());


  AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> TbulkDynamics (
    converter.getOmegaT(),
    instances::getBulkMomenta<T,TDESCRIPTOR>()
  );

  couplingBetweenNavierStokesAndAdvectionDiffusion(ADlattice, NSlattice, converter);

  iniGeometry(NSlattice, ADlattice,
              NSbulkDynamics, TbulkDynamics,
              *NSboundaryCondition, *TboundaryCondition,
              converter );

  util::ValueTracer<T> converge(0.01,(T)N,1.0e-5);
  int iT = 0;
  clout << "starting simulation" << endl;
  for (iT=0; iT<maxIter; ++iT)
  {
    ADlattice.collideAndStream(true);
    NSlattice.collideAndStream(true);

    converge.takeValue(NSlattice.getStatistics().getAverageEnergy(),true);

    NSlattice.executeCoupling();
    ADlattice.executeCoupling();
    if (converge.hasConverged())
    {
      clout << "Simulation converged." << endl;
      break;
    }
  }

  clout << "Nusselt = " << computeNusselt(NSlattice, ADlattice, converter) << endl;

  delete NSboundaryCondition;
  delete TboundaryCondition;
}

