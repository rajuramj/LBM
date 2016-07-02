/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
 *  Address: EPFL, STI-LIN Station 9, 1015 Lausanne, Switzerland
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

#include "olb3D.h"
#include "olb3D.hh"
#include <cstdlib>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR ForcedShanChenD3Q19Descriptor


void iniGeometry( BlockStructure3D<T, DESCRIPTOR>& latticeOne,
                  BlockStructure3D<T, DESCRIPTOR>& latticeTwo,
                  Dynamics<T, DESCRIPTOR>& bulkDynamics1,
                  Dynamics<T, DESCRIPTOR>& bulkDynamics2,
                  Dynamics<T, DESCRIPTOR>& bounceBackRho0,
                  Dynamics<T, DESCRIPTOR>& bounceBackRho1,
                  T force
                )
{
  // The setup is: periodicity along horizontal direction, bounce-back on top
  // and bottom. The upper half is initially filled with fluid 1 + random noise,
  // and the lower half with fluid 2. Only fluid 1 experiences a forces,
  // directed downwards.
  T noise      = 1.e-2;

  int nx = latticeOne.getNx();
  int ny = latticeOne.getNy();
  int nz = latticeOne.getNz();

  latticeOne.defineDynamics(0,nx-1, 0,ny-1, 0,nz-1, &bulkDynamics1);
  latticeTwo.defineDynamics(0,nx-1, 0,ny-1, 0,nz-1, &bulkDynamics2);

  latticeOne.defineDynamics(0,nx-1, 0,ny-1, 0,0, &bounceBackRho0);
  latticeTwo.defineDynamics(0,nx-1, 0,ny-1, 0,0, &bounceBackRho1);
  latticeOne.defineDynamics(0,nx-1, 0,ny-1, nz-1,nz-1, &bounceBackRho1);
  latticeTwo.defineDynamics(0,nx-1, 0,ny-1, nz-1,nz-1, &bounceBackRho0);

  for (int iX=0; iX<nx; ++iX) {
    for (int iY=0; iY<ny; ++iY) {
      for (int iZ=0; iZ<nz; ++iZ) {
        T f[3] = {0.,0.,-force*(1.+(double)random()/(double)RAND_MAX * noise)};
        T zeroVelocity[3] = {0.,0.,0.};
        T noForce[3] = {0.,0.,0.};
        T rho1 = (T)1;
        T rho2 = (T)1;
        if (iZ>nz/2) {
          rho2 = 0.;
          //rho1 += (double)random()/(double)RAND_MAX * noise;
        }
        else {
          rho1 = 0.;
        }

        Cell<T,DESCRIPTOR>& cell1 = latticeOne.get(iX,iY,iZ);
        cell1.defineRhoU(rho1, zeroVelocity);
        cell1.iniEquilibrium(rho1, zeroVelocity);
        cell1.defineExternalField (
          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
          DESCRIPTOR<T>::ExternalField::sizeOfForce, f );

        Cell<T,DESCRIPTOR>& cell2 = latticeTwo.get(iX,iY,iZ);
        cell2.defineRhoU(rho2, zeroVelocity);
        cell2.iniEquilibrium(rho2, zeroVelocity);
        cell2.defineExternalField (
          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
          DESCRIPTOR<T>::ExternalField::sizeOfForce, noForce );
      }
    }
  }

  latticeOne.initialize();
  latticeTwo.initialize();
}

void plotStatistics(BlockStructure3D<T, DESCRIPTOR>& latticeTwo,
                    BlockStructure3D<T, DESCRIPTOR>& latticeOne, int iT)
{
  OstreamManager cout(std::cout,"plotStatistics");

  cout << "Writing Gif..." << std::endl;
  int nx = latticeOne.getNx();
  ImageWriter<T> imageCreator("leeloo.map");

  imageCreator.writeScaledGif(createFileName("p", iT, 6),
                              latticeOne.getDataAnalysis().getPressure().sliceX(nx/2),
                              400,400);
}

void produceVTK(BlockStructure3D<T, DESCRIPTOR>& latticeOne, int iT) {
  OstreamManager cout(std::cout,"produceVTK");

  cout << "Writing VTK" << std::endl;

  VtkImageOutput3D<T> vtkOut(createFileName("vtk", iT, 7), (T)1);
  vtkOut.writeData<T,float>(latticeOne.getDataAnalysis().getPressure(),
                            "p", (T)1);
}

int main(int argc, char *argv[])
{
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");

  const T omega1 = 1.0;
  const T omega2 = 1.0;
  const int nx   = 75;
  const int ny   = 75;
  const int nz   = 75;
  const T G      = 0.1;
  T force        = 2.e-3;

  const int maxIter  = 1000000;
  const int saveIter = 100;
  const int vtkIter  = 200;
  const int statIter = 10;


#ifndef PARALLEL_MODE_MPI  // sequential program execution
  MultiDataDistribution3D distribution = createRegularDataDistribution( nx, ny, nz, 1, 1, 1, 1 );
#else                      // parallel program execution
  MultiDataDistribution3D distribution = createRegularDataDistribution( nx, ny, nz );
#endif
  // Use a multi-block lattice even in sequential program execution. In this
  //   way, it is easier to have periodic boundaries, because one can
  //   access periodic data from the lattice envelope
  MultiBlockLattice3D<T, DESCRIPTOR> latticeOne(distribution);
  MultiBlockLattice3D<T, DESCRIPTOR> latticeTwo(distribution);

  BGKdynamics<T, DESCRIPTOR> bulkDynamics1 (
    omega1, instances::getExternalVelocityMomenta<T,DESCRIPTOR>() );
  BGKdynamics<T, DESCRIPTOR> bulkDynamics2 (
    omega2, instances::getExternalVelocityMomenta<T,DESCRIPTOR>() );
  // A bounce-back node with fictitious density 1,
  //   which is experienced by the partner fluid
  BounceBack<T, DESCRIPTOR> bounceBackRho1( 1. );
  // A bounce-back node with fictitious density 0,
  //   which is experienced by the partner fluid
  BounceBack<T, DESCRIPTOR> bounceBackRho0( 0. );

  vector<SpatiallyExtendedObject3D* > partnerForOne;
  partnerForOne.push_back(&latticeTwo);

  ForcedShanChenCouplingGenerator3D<T,DESCRIPTOR> coupling(0,nx-1,0,ny-1, 1,nz-2,G);
  latticeOne.addLatticeCoupling(coupling, partnerForOne);

  iniGeometry(latticeOne, latticeTwo, bulkDynamics1, bulkDynamics2,
              bounceBackRho0, bounceBackRho1, force);

  int iT = 0;
  clout << "starting simulation..." << endl;
  for (iT=0; iT<maxIter; ++iT) {
    if (iT%statIter==0 && iT > 0) {
      clout << "averageRhoFluidOne="   << latticeOne.getStatistics().getAverageRho();
      clout << "; averageRhoFluidTwo=" << latticeTwo.getStatistics().getAverageRho() << std::endl;
    }
    if (iT%saveIter==0) {
      plotStatistics(latticeTwo, latticeOne, iT);
    }
    if (iT%vtkIter==0) {
      produceVTK(latticeOne, iT);
    }

    latticeOne.collideAndStream(true);
    latticeTwo.collideAndStream(true);

    latticeOne.executeCoupling();
    latticeTwo.executeCoupling();
  }

}

