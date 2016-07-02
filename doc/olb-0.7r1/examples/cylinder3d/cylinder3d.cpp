/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011 Mathias J. Krause, Thomas Henn
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
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

/* cylinder3d.cpp:
 * This example examines a steady flow past a cylinder placed in a
 * channel. At the inlet, a Poiseuille profile is imposed on the velocity,
 * where as the outlet implements an outflow condition: grad_n u = 0 and
 * rho = 1. It also shows the usage of the STL-reader and explains how
 * to set boundary conditions automatically.
 */

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
#include "olb3D.hh"     // Include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor

/// Parameters for the simulation setup
const T uMax            = 0.05;    /// uMax for the converter
const int iTMax         = 10000;   /// no of time steps for smooth start-up


void iniGeometry(SuperLattice3D<T, D3Q19Descriptor>& lattice,
                 LBconverter<T> const& converter, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& bc,
                 BlockGeometry3D& blockGeometry) {

  const int nx = blockGeometry.getNx();
  const int ny = blockGeometry.getNy();
  const int nz = blockGeometry.getNz();
  const T omega = converter.getOmega();

  /// Defining the statistics
  BlockGeometryStatistics3D bloGeSta(&blockGeometry);
  if (singleton::mpi().isMainProcessor()) {
    bloGeSta.countVoxel();
  }
  /// Material=0 -->do nothing
  lattice.defineDynamics(&bloGeSta,
                         &instances::getNoDynamics<T, DESCRIPTOR>(), 0);

  /// Material=1 -->bulk dynamics
  lattice.defineDynamics(&bloGeSta, &bulkDynamics, 1);

  /// Material=2 -->bounce back
  lattice.defineDynamics(&bloGeSta, &bulkDynamics, 2);
  bc.addVelocityBoundary(&bloGeSta, omega, 2);

  /// Material=3 -->bulk dynamics (inflow)
  lattice.defineDynamics(&bloGeSta, &bulkDynamics, 3);

  /// Material=4 -->bulk dynamics (outflow)
  lattice.defineDynamics(&bloGeSta, &bulkDynamics, 4);

  /// Setting of the boundary conditions
  bc.addVelocityBoundary(&bloGeSta, omega, 3);
  bc.addPressureBoundary(&bloGeSta, omega, 4);

  for (int iX = 0; iX < nx; ++iX) {
    for (int iY = 0; iY < ny; ++iY) {
      for (int iZ = 0; iZ < nz; ++iZ) {
        T vel[] = { T(), T(), T() };
        lattice.defineRhoU(iX, iX, iY, iY, iZ, iZ, (T) 1, vel);
        lattice.iniEquilibrium(iX, iX, iY, iY, iZ, iZ, (T) 1, vel);
      }
    }
  }

  /// Computation of the poiseuille inflow
  T H = ny - 1;
  for (int iY = 0; iY < ny; ++iY) {
    for (int iZ = 0; iZ < nz; ++iZ) {
      T u = converter.getLatticeU() * 9 / 4 * 16 * iY / H * iZ / H * (1
            - iY / H) * (1 - iZ / H);
      T vel[] = { u, T(), T() };
      lattice.defineRhoU(0, 0, iY, iY, iZ, iZ, (T) 1, vel);
      lattice.iniEquilibrium(0, 0, iY, iY, iZ, iZ, (T) 1, vel);
    }
  }

  /// Make the lattice ready for simulation
  lattice.initialize();
  return;
}

/// Generates a slowly increasing sinuidal inflow for the first iTMax timesteps
void reIniGeometry(SuperLattice3D<T, D3Q19Descriptor>& lattice,
                   LBconverter<T> const& converter, int iT,
                   Dynamics<T, D3Q19Descriptor>& bulkDynamics,
                   BlockGeometry3D& blockGeometry) {

  OstreamManager clout(std::cout,"reIniGeometry");

  if (iT <= iTMax) {

    const int ny = blockGeometry.getNy();
    const int nz = blockGeometry.getNz();

    T frac = (sin(-3.1416 / 2.0 + (T) iT / (T) iTMax * 3.1416) + 1.0) / 2.0;

    /// Computation of the poiseuille inflow
    T H = ny - 1;
    for (int iY = 0; iY < ny; ++iY) {
      for (int iZ = 0; iZ < nz; ++iZ) {
        T u = converter.getLatticeU() * 9 / 4 * 16 * iY / H * iZ / H
              * (1 - iY / H) * (1 - iZ / H) * frac;
        T vel[] = { u, T(), T() };
        lattice.defineU(0, 0, iY, iY, iZ, iZ, vel);
      }
    }

    if (iT % 50 == 0)
      clout << "step=" << iT << "; scalingFactor=" << frac << std::endl;
  }
}

/// Computes the pressure drop between the voxels before and after the cylinder
void computeResults(SuperLattice3D<T, D3Q19Descriptor>& lattice,
                    LBconverter<T> const& converter, int iter,
                    Dynamics<T, D3Q19Descriptor>& bulkDynamics,
                    BlockGeometry3D& blockGeometry) {

  OstreamManager clout(std::cout,"computeResults");
  Cell<T, D3Q19Descriptor> cell;

  T iLu, iLl, iRu, iRl, interpIl, interpIr;

  /// Interpolation of 4 inner nodes
  if (lattice.get(45, 20, 21, cell)) {
    T rho;
    T vel[] = { T(), T(), T() };
    cell.defineDynamics(&bulkDynamics);
    cell.computeRhoU(rho, vel);
    rho = (rho - 1.) / 3.;
    iLu = rho * 1 * 0.2 / uMax * 0.2 / uMax;
  }

  if (lattice.get(45, 20, 20, cell)) {
    T rho;
    T vel[] = { T(), T(), T() };
    cell.defineDynamics(&bulkDynamics);
    cell.computeRhoU(rho, vel);
    rho = (rho - 1.) / 3.;
    iLl = rho * 1 * 0.2 / uMax * 0.2 / uMax;
  }
  interpIl = (iLu + iLl) / 2;

  if (lattice.get(56, 20, 21, cell)) {
    T rho;
    T vel[] = { T(), T(), T() };
    cell.defineDynamics(&bulkDynamics);
    cell.computeRhoU(rho, vel);
    rho = (rho - 1.) / 3.;
    iRu = rho * 1 * 0.2 / uMax * 0.2 / uMax;
  }

  if (lattice.get(56, 20, 20, cell)) {
    T rho;
    T vel[] = { T(), T(), T() };
    cell.defineDynamics(&bulkDynamics);
    cell.computeRhoU(rho, vel);
    rho = (rho - 1.) / 3.;
    iRl = rho * 1 * 0.2 / uMax * 0.2 / uMax;
  }
  interpIr = (iRu + iRl) / 2;
  clout << "deltaP=" << interpIl - interpIr << endl;
}

/// Writing the vtk files
void writeVTK(SuperLattice3D<T, D3Q19Descriptor>& sLattice,
              LBconverter<T> const& converter, int iter) {
  vector<const ScalarFieldBase3D<T>*> scalar;
  vector<const TensorFieldBase3D<T, 3>*> tensor;

  for (int iC = 0; iC < sLattice.get_load().size(); iC++) {
    const TensorFieldBase3D<T, 3>* velocity =  &sLattice.get_lattice(iC).getDataAnalysis().getVelocity();
    const ScalarFieldBase3D<T>*    pressure = &sLattice.get_lattice(iC).getDataAnalysis().getPressure();
    scalar.push_back(pressure);
    tensor.push_back(velocity);
  }
  CuboidVTKout3D<T>::writeFlowField(createFileName("vtkCylinder", iter, 6),
                                    "Pressure", scalar, "Velocity", tensor, sLattice.get_cGeometry(),
                                    sLattice.get_load(), converter.getDeltaT());
}

int main(int argc, char* argv[]) {

  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  int TaskID;
#ifdef PARALLEL_MODE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&TaskID);
#endif

  /// Instantiation of the STLreader class
  STLreader<T> stlreader("cylinder.stl");

  /// Instantiation of a empty blockGeometry
  BlockGeometry3D blockGeometry;

  /// Sets the material number outside the surface
  stlreader.setOuterMaterialNo(2);

  /// Reads the stl-file and stores it in blockGeometry
  stlreader.read(blockGeometry, 0u, 250);

  if (singleton::mpi().isMainProcessor()) {
    clout << "Nx=" << blockGeometry.getNx() << "; Ny="  << blockGeometry.getNy() << "; Nz=" << blockGeometry.getNz() << endl;
  }

  /// Set material number for inflow
  for (int iY = 1; iY < blockGeometry.getNy() - 1; iY++) {
    for (int iZ = 1; iZ < blockGeometry.getNz() - 1; iZ++) {
      blockGeometry.setMaterial(0, iY, iZ, 3);
    }
  }
  /// Set material number for outflow
  for (int iY = 1; iY < blockGeometry.getNy() - 1; iY++) {
    for (int iZ = 1; iZ < blockGeometry.getNz() - 1; iZ++) {
      blockGeometry.setMaterial(251, iY, iZ, 4);
    }
  }

  /// Removes all not needed boundary voxels outside the surface
  blockGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  blockGeometry.innerClean();
  blockGeometry.checkForErrors();

  /// Writes the blockGeometry as vti file for visualization
  blockGeometry.writeVti("geometry");

  LBconverter<T> converter(
    (int) 3,                               // dim
    (T)   1.,                              // latticeL_
    (T)   uMax,                            // latticeU_
    (T)   1./2.,                           // charNu_
    (T)   10.                              // charL_ = 1,
  );
  writeLogFile(converter, "cylinder3d");

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cGeometry(blockGeometry.getPositionX(), blockGeometry.getPositionY(), blockGeometry.getPositionZ(), blockGeometry.getSpacing(), blockGeometry.getNx(), blockGeometry.getNy(), blockGeometry.getNz(), singleton::mpi().getSize() );
#else
  CuboidGeometry3D<T> cGeometry(blockGeometry.getPositionX(), blockGeometry.getPositionY(), blockGeometry.getPositionZ(), blockGeometry.getSpacing(), blockGeometry.getNx(), blockGeometry.getNy(), blockGeometry.getNz(), 7 );
#endif

  SuperLattice3D<T, DESCRIPTOR> sLattice(cGeometry, 1);

  BGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>());
  sOnLatticeBoundaryCondition3D<T, D3Q19Descriptor> sBoundaryCondition(sLattice);
  createInterpBoundaryCondition3D<T, D3Q19Descriptor> (sBoundaryCondition);

  iniGeometry(sLattice, converter, bulkDynamics, sBoundaryCondition, blockGeometry);

  for (int iT = 0; iT <= 20000; ++iT) {

    if (iT % 50 == 0 && singleton::mpi().isMainProcessor()) {
      sLattice.getStatistics().print(iT,iT*converter.getDeltaT());
    }

    if (iT % 50 == 0 && iT > 0) {
      computeResults(sLattice, converter, iT, bulkDynamics, blockGeometry);
    }

    if (iT % 1000 == 0 && iT > 0) {
      cout << "Saving VTK file ..." << endl;
      writeVTK(sLattice, converter, iT);
    }
    sLattice.collideAndStream();

    reIniGeometry(sLattice, converter, iT, bulkDynamics, blockGeometry);
  }
}
