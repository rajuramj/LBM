/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
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

/** \file
 * The description of a 3D super lattice -- generic implementation.
 */

#ifndef SUPER_LATTICE_3D_HH
#define SUPER_LATTICE_3D_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "core/blockLattice3D.h"
#include "core/cell.h"
#include "communicator3D.h"
#include "cuboidGeometry3D.h"
#include "core/loadBalancer.h"
#include "core/postProcessing.h"
#include "superLattice3D.h"

namespace olb {

////////////////////// Class SuperLattice3D /////////////////////////


template<typename T, template<typename U> class Lattice>
SuperLattice3D<T, Lattice>::SuperLattice3D(CuboidGeometry3D<T>& cGeometry,
    int overlapBC, loadBalancer *lb) :
  _cGeometry(cGeometry), _commStream(*this), _commBC(*this),
  _load(lb) {

  if (overlapBC >= 1) {
    _commBC_on = true;
    _overlap = overlapBC;
  } else {
    _commBC_on = false;
    _overlap = 1;
  }

  int rank = 0;
  int size = 1;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = singleton::mpi().getSize();
#endif

  if(!_load) {
    _load = new loadBalancer(rank, size, _cGeometry.get_nC(), 0);
  }
  for (int iC = 0; iC < _load->size(); iC++) {
    int nX = _cGeometry.get_cuboid(_load->glob(iC)).get_nX() + 2 * _overlap;
    int nY = _cGeometry.get_cuboid(_load->glob(iC)).get_nY() + 2 * _overlap;
    int nZ = _cGeometry.get_cuboid(_load->glob(iC)).get_nZ() + 2 * _overlap;

    BlockLattice3D<T, Lattice> tmp(0, 0, 0);
    _blockLattices.push_back(tmp);
    BlockLattice3D<T, Lattice> blockLattice(nX, nY, nZ);
    _blockLattices.back().swap(blockLattice);
  }
  for (int iC = 0; iC < _load->size(); iC++) {
    BlockLatticeView3D<T, Lattice> lattice(_blockLattices[iC], _overlap,
                                           _blockLattices[iC].getNx() - _overlap - 1, _overlap,
                                           _blockLattices[iC].getNy() - _overlap - 1, _overlap,
                                           _blockLattices[iC].getNz() - _overlap - 1);
    _lattices.push_back(lattice);
  }
  _statistics = new LatticeStatistics<T> ;
  _statistics_on = true;

  _commStream.init_nh();
  _commStream.add_cells(1);
  _commStream.init();

  if (_commBC_on)
    _commBC.init_nh();
}

template<typename T, template<typename U> class Lattice>
bool SuperLattice3D<T, Lattice>::set(T iX, T iY, T iZ,
                                     Cell<T, Lattice> const& cell) {

  bool found = false;
  int locX, locY, locZ;
  for (int iC = 0; iC < _load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkPoint(iX, iY, iZ, locX,
        locY, locZ, _overlap)) {
      _blockLattices[iC].get(locX, locY, locZ) = cell;
      found = true;
    }
  }
  return found;
}

template<typename T, template<typename U> class Lattice>
bool SuperLattice3D<T, Lattice>::get(T iX, T iY, T iZ, Cell<T, Lattice>& cell) const {

  int locX, locY, locZ;
  bool found = false;
  int foundIC = 0;

  for (int iC = 0; iC < _cGeometry.get_nC(); iC++) {
    if (_cGeometry.get_cuboid(iC).checkPoint(iX, iY, iZ, locX, locY, locZ)) {
      found = true;
      foundIC = iC;
    }
  }

#ifdef PARALLEL_MODE_MPI
  const int sizeOfCell = Lattice<T>::q + Lattice<T>::ExternalField::numScalars;
  T* cellData = new T[sizeOfCell];

  if (found) {
    if (_load->rank(foundIC)==singleton::mpi().getRank()) {
      _lattices[_load->loc(foundIC)].get(locX,locY,locZ).serialize(cellData);
    }
    singleton::mpi().bCast(cellData, sizeOfCell, _load->rank(foundIC));
    cell.unSerialize(cellData);
    delete [] cellData;
  }
#else
  if (found) {
    cell = _lattices[_load->loc(foundIC)].get(locX, locY, locZ);
  }
#endif

  return found;
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::initialize() {

  if (_commBC_on) {
    _commBC.init();
  }
  for (int iC = 0; iC < _load->size(); iC++) {
    //AENDERN VON INI in BLOCKLATTICEVIEW!!!!
    //_lattices[iC].initialize();
    _lattices[iC].postProcess();
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::defineDynamics(T x0, T x1, T y0, T y1, T z0,
    T z1, Dynamics<T, Lattice>* dynamics) {

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < _load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      _blockLattices[iC].defineDynamics(locX0, locX1, locY0, locY1,
                                        locZ0, locZ1, dynamics);
    }
  }
}

template<typename T, template<typename U> class Lattice> void SuperLattice3D<T,
         Lattice>::defineDynamics(BlockGeometryStatistics3D* blockGeoSta, T x0,
                                  T x1, T y0, T y1, T z0, T z1, Dynamics<T, Lattice>* dynamics,
int material) {
  int locX, locY, locZ;
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {

        if (blockGeoSta->getBlockGeometry()->getMaterial(iX, iY, iZ)
        == material) {
          for (int iC = 0; iC < _load->size(); iC++) {
            if (_cGeometry.get_cuboid(_load->glob(iC)).checkPoint(
            iX, iY, iZ, locX, locY, locZ, _overlap)) {
              _blockLattices[iC].defineDynamics(locX, locX, locY,
                                                locY, locZ, locZ, dynamics);
            }
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::defineDynamics(
  BlockGeometryStatistics3D* blockGeoSta, Dynamics<T, Lattice>* dynamics,
  int material) {

  defineDynamics(blockGeoSta,
                 -_overlap, blockGeoSta->getBlockGeometry()->getNx() - 1+_overlap,
                 -_overlap, blockGeoSta->getBlockGeometry()->getNy() - 1+_overlap,
                 -_overlap, blockGeoSta->getBlockGeometry()->getNz() - 1+_overlap,
                 dynamics, material);

}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::defineRhoU(T x0, T x1, T y0, T y1, T z0, T z1,
    T rho, const T u[Lattice<T>::d]) {
  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < _load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      for (int iX = locX0; iX <= locX1; iX++) {
        for (int iY = locY0; iY <= locY1; iY++) {
          for (int iZ = locZ0; iZ <= locZ1; iZ++) {
            _blockLattices[iC].get(iX, iY, iZ).defineRhoU(rho, u);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T,Lattice>::defineRho (
  T x0, T x1, T y0, T y1, T z0, T z1, T rho ) {

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC=0; iC<_load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters( x0, x1, y0, y1, z0, z1,
        locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      for (int iX=locX0; iX<=locX1; iX++) {
        for (int iY=locY0; iY<=locY1; iY++) {
          for (int iZ=locZ0; iZ<=locZ1; iZ++) {
            _blockLattices[iC].get(iX,iY,iZ).defineRho(rho);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T,Lattice>::defineU (
  T x0, T x1, T y0, T y1, T z0, T z1, const T u[Lattice<T>::d] ) {

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC=0; iC<_load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters( x0, x1, y0, y1, z0, z1,
        locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      for (int iX=locX0; iX<=locX1; iX++) {
        for (int iY=locY0; iY<=locY1; iY++) {
          for (int iZ=locZ0; iZ<=locZ1; iZ++) {
            _blockLattices[iC].get(iX,iY,iZ).defineU(u);
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T,Lattice>::defineExternalField (
  T x0, T x1, T y0, T y1, T z0, T z1, int fieldBeginsAt, int sizeOfField,
  T* field ) {

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC=0; iC<_load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters( x0, x1, y0, y1, z0,
        z1, locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      for (int iX=locX0; iX<=locX1; iX++) {
        for (int iY=locY0; iY<=locY1; iY++) {
          for (int iZ=locZ0; iZ<=locZ1; iZ++) {
            _blockLattices[iC].get(iX,iY,iZ).defineExternalField (
              fieldBeginsAt, sizeOfField, field);
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::iniEquilibrium(T x0, T x1, T y0, T y1, T z0,
    T z1, T rho, const T u[Lattice<T>::d]) {
  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < _load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      for (int iX = locX0; iX <= locX1; iX++) {
        for (int iY = locY0; iY <= locY1; iY++) {
          for (int iZ = locZ0; iZ <= locZ1; iZ++) {
            _blockLattices[iC].get(iX, iY, iZ).iniEquilibrium(rho,
                u);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::collide() {

  for (int iC = 0; iC < _load->size(); iC++) {
    _lattices[iC].collide();
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::collide(T x0, T x1, T y0, T y1, T z0, T z1) {

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < _load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1)) {
      _lattices[iC].collide(locX0, locX1, locY0, locY1, locZ0, locZ1);
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::stream() {

  _commStream.send();
  _commStream.receive();
  _commStream.write();

  for (int iC = 0; iC < _load->size(); iC++) {
    _blockLattices[iC].stream(_overlap - 1, _blockLattices[iC].getNx()
                              - _overlap, _overlap - 1,
                              _blockLattices[iC].getNy() - _overlap, _overlap - 1,
                              _blockLattices[iC].getNz() - _overlap);
  }
  if (_commBC_on) {
    _commBC.send();
    _commBC.receive();
    _commBC.write();
  }

  for (int iC = 0; iC < _load->size(); iC++) {
    _lattices[iC].postProcess();
  }
  if (_statistics_on)
    reset_statistics();
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::stream(T x0, T x1, T y0, T y1, T z0, T z1) {

  _commStream.send();
  _commStream.receive();
  _commStream.write();

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < _load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      _blockLattices[iC].stream(locX0, locX1, locY0, locY1, locZ0, locZ1);
    }
  }
  if (_commBC_on) {
    _commBC.send();
    _commBC.receive();
    _commBC.write();
  }

  for (int iC = 0; iC < _load->size(); iC++) {
    _lattices[iC].postProcess();
  }
  if (_statistics_on)
    reset_statistics();
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::collideAndStream() {

  for (int iC = 0; iC < _load->size(); iC++) {
    int x1 = _lattices[iC].getNx() - 1;
    int y1 = _lattices[iC].getNy() - 1;
    int z1 = _lattices[iC].getNz() - 1;

    _lattices[iC].collide(0, x1, 0, 0, 0, z1);
    _lattices[iC].collide(0, x1, y1, y1, 0, z1);
    _lattices[iC].collide(0, 0, 1, y1 - 1, 0, z1);
    _lattices[iC].collide(x1, x1, 1, y1 - 1, 0, z1);
    _lattices[iC].collide(1, x1 - 1, 1, y1 - 1, 0, 0);
    _lattices[iC].collide(1, x1 - 1, 1, y1 - 1, z1, z1);
  }

  for (int iC = 0; iC < _load->size(); iC++) {
    _blockLattices[iC].bulkCollideAndStream(_overlap + 1,
                                            _blockLattices[iC].getNx() - _overlap - 2, _overlap + 1,
                                            _blockLattices[iC].getNy() - _overlap - 2, _overlap + 1,
                                            _blockLattices[iC].getNz() - _overlap - 2);
  }

  _commStream.send();
  _commStream.receive();
  _commStream.write();

  for (int iC = 0; iC < _load->size(); iC++) {
    int x1 = _blockLattices[iC].getNx() - 1;
    int y1 = _blockLattices[iC].getNy() - 1;
    int z1 = _blockLattices[iC].getNz() - 1;

    _blockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1, _overlap - 1, x1
                                      - _overlap + 1, _overlap - 1, _overlap, _overlap - 1, z1
                                      - _overlap + 1);
    _blockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1, _overlap - 1, x1
                                      - _overlap + 1, y1 - _overlap, y1 - _overlap + 1, _overlap - 1,
                                      z1 - _overlap + 1);
    _blockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1, _overlap - 1,
                                      _overlap, _overlap + 1, y1 - _overlap - 1, _overlap - 1, z1
                                      - _overlap + 1);
    _blockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1, x1 - _overlap,
                                      x1 - _overlap + 1, _overlap + 1, y1 - _overlap - 1, _overlap
                                      - 1, z1 - _overlap + 1);
    _blockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1, _overlap + 1, x1
                                      - _overlap - 1, _overlap + 1, y1 - _overlap - 1, _overlap - 1,
                                      _overlap);
    _blockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1, _overlap + 1, x1
                                      - _overlap - 1, _overlap + 1, y1 - _overlap - 1, z1 - _overlap,
                                      z1 - _overlap + 1);
  }

  if (_commBC_on) {
    //_commBC.send();
    //_commBC.receive();
    //_commBC.write();
    _commStream.send();
    _commStream.receive();
    _commStream.write();
  }

  for (int iC = 0; iC < _load->size(); iC++) {
    _lattices[iC].postProcess();
  }
  if (_statistics_on)
    reset_statistics();
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T,Lattice>::stripeOffDensityOffset (
  int x0, int x1, int y0, int y1, int z0, int z1, T offset )
{

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < _load->size(); iC++) {
    if (_cGeometry.get_cuboid(_load->glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, _overlap)) {
      _blockLattices[iC].stripeOffDensityOffset(locX0, locX1, locY0, locY1,
          locZ0, locZ1, offset);
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T,Lattice>::stripeOffDensityOffset(T offset) {

  for (int iC = 0; iC < _load->size(); iC++) {
    _blockLattices[iC].stripeOffDensityOffset(offset);
  }
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T>& SuperLattice3D<T, Lattice>::getStatistics() {

  return *_statistics;
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T> const&
SuperLattice3D<T, Lattice>::getStatistics() const {

  return *_statistics;
}

template<typename T, template<typename U> class Lattice>
void SuperLattice3D<T, Lattice>::reset_statistics() {

  T weight;
  T sum_weight = 0;
  T average_rho = 0;
  T average_energy = 0;
  T maxU = 0;
  T delta = 0;

  getStatistics().reset();

  for (int iC = 0; iC < _load->size(); iC++) {
    delta = _cGeometry.get_cuboid(_load->glob(iC)).get_delta();
    weight = _blockLattices[iC].getStatistics().getNumCells() * delta
             * delta * delta;
    sum_weight += weight;
    average_rho += _blockLattices[iC].getStatistics().getAverageRho()
                   * weight;
    average_energy += _blockLattices[iC].getStatistics().getAverageEnergy()
                      * weight;
    if (maxU < _blockLattices[iC].getStatistics().getMaxU() / delta) {
      maxU = _blockLattices[iC].getStatistics().getMaxU() / delta;
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(sum_weight, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_rho, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_energy, MPI_SUM);
  singleton::mpi().reduceAndBcast(maxU, MPI_MAX);
#endif

  average_rho = average_rho / sum_weight;
  average_energy = average_energy / sum_weight;

  getStatistics().reset(average_rho, average_energy, maxU, (int) sum_weight);
  for (int iC = 0; iC < _load->size(); iC++) {
    delta = _cGeometry.get_cuboid(_load->glob(iC)).get_delta();
    _blockLattices[iC].getStatistics().reset(average_rho, average_energy,
        maxU * delta, (int) sum_weight);
  }
}

} // namespace olb

#endif
