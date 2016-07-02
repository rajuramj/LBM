/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2011 Jonas Latt, Mathias J. Krause,
 *  Jonas Kratzke
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

#ifndef UNITS_HH
#define UNITS_HH

#include "units.h"

namespace olb {

template<typename T>
void writeLogFile(LBunits<T> const& converter,
                  std::string const& title)
{
  std::string fullName = singleton::directories().getLogOutDir() + "olbLog.dat";
  olb_ofstream ofile(fullName.c_str());
  ofile << title << "\n\n";
  ofile << "Velocity in lattice units: u="         << converter.getLatticeU()   << "\n";
  ofile << "Reynolds number:           Re="        << converter.getRe()         << "\n";
  ofile << "Lattice viscosity:         latticeNu=" << converter.getLatticeNu()  << "\n";
  ofile << "Lattice resolution:        N="         << converter.getResolution() << "\n";
  ofile << "Extent of the system:      lx="        << converter.getLx()         << "\n";
  ofile << "Extent of the system:      ly="        << converter.getLy()         << "\n";
  ofile << "Extent of the system:      lz="        << converter.getLz()         << "\n";
  ofile << "Grid spacing deltaX:       dx="        << converter.getDeltaX()     << "\n";
  ofile << "Time step deltaT:          dt="        << converter.getDeltaT()     << "\n";
  ofile << "Relaxation time:           tau="       << converter.getTau()        << "\n";
  ofile << "Relaxation frequency:      omega="     << converter.getOmega()      << "\n";
}

template<typename T>
void writeLogFile(LBconverter<T> const& converter,
                  std::string const& title)
{
  std::string fullName = singleton::directories().getLogOutDir() + title + ".dat";
  olb_ofstream ofile(fullName.c_str());
  ofile << "LBconverter information\n\n";
  ofile << "characteristical values\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "Dimension(d):                     dim="              << converter.getDim()              << "\n";
  ofile << "Characteristical length(m):       charL="            << converter.getCharL()            << "\n";
  ofile << "Characteristical speed(m/s):      charU="            << converter.getCharU()            << "\n";
  ofile << "Characteristical time(s):         charT="            << converter.getCharTime()         << "\n";
  ofile << "Density factor(kg/m^d):           charRho="          << converter.getCharRho()          << "\n";
  ofile << "Characterestical mass(kg):        charMass="         << converter.getCharMass()         << "\n";
  ofile << "Characterestical force(N):        charForce="        << converter.getCharForce()        << "\n";
  ofile << "Characterestical pressure(Pa):    charPressure="     << converter.getCharPressure()     << "\n";
  ofile << "Pressure level(Pa):               pressureLevel="    << converter.getPressureLevel()    << "\n";
  ofile << "Phys. kinematic viscosity(m^2/s): charNu="           << converter.getCharNu()           << "\n";
  ofile << "======================================================================\n\n";
  ofile << "lattice values\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "DeltaX:                           deltaX="           << converter.getDeltaX()           << "\n";
  ofile << "Lattice velocity:                 latticeU="         << converter.getLatticeU()         << "\n";
  ofile << "DeltaT:                           deltaT="           << converter.getDeltaT()           << "\n";
  ofile << "Reynolds number:                  Re="               << converter.getRe()               << "\n";
  ofile << "DimlessNu:                        dNu="              << converter.getDimlessNu()        << "\n";
  ofile << "Viscosity for computation:        latticeNu="        << converter.getLatticeNu()        << "\n";
  ofile << "Relaxation time:                  tau="              << converter.getTau()              << "\n";
  ofile << "Relaxation frequency:             omega="            << converter.getOmega()            << "\n";
  ofile << "======================================================================\n\n";
  ofile << "conversion factors\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "latticeL(m):                      latticeL="         << converter.getLatticeL()         << "\n";
  ofile << "Time step (s):                    physTime="         << converter.physTime()            << "\n";
  ofile << "Velocity factor(m/s):             physVelocity="     << converter.physVelocity()        << "\n";
  ofile << "FlowRate factor(m^d/s):           physFlowRate="     << converter.physFlowRate()        << "\n";
  ofile << "Mass factor(kg):                  physMass="         << converter.physMass()            << "\n";
  ofile << "Force factor(N):                  physForce="        << converter.physForce()           << "\n";
  ofile << "Force factor massless(N/kg):      physMasslessForce="<< converter.physMasslessForce()   << "\n";
  ofile << "Pressure factor(Pa):              physPressure="     << converter.physPressure(4)       << "\n";
  ofile << "latticePressure:                  latticeP="         << converter.latticePressure()     << "\n";

}

template<typename T>
void write(LBconverter<T> const& converter)
{
  OstreamManager clout(std::cout,"ConvLog");

  clout << "LBconverter information" << std::endl;
  clout << "characteristical values" << std::endl;
  clout << "Dimension(d):                     dim="              << converter.getDim()              << std::endl;
  clout << "Characteristical length(m):       charL="            << converter.getCharL()            << std::endl;
  clout << "Characteristical speed(m/s):      charU="            << converter.getCharU()            << std::endl;
  clout << "Characteristical time(s):         charT="            << converter.getCharTime()         << std::endl;
  clout << "Density factor(kg/m^d):           charRho="          << converter.getCharRho()          << std::endl;
  clout << "Characterestical mass(kg):        charMass="         << converter.getCharMass()         << std::endl;
  clout << "Characterestical force(N):        charForce="        << converter.getCharForce()        << std::endl;
  clout << "Characterestical pressure(Pa):    charPressure="     << converter.getCharPressure()     << std::endl;
  clout << "Pressure level(Pa):               pressureLevel="    << converter.getPressureLevel()    << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s): charNu="           << converter.getCharNu()           << std::endl;

  clout << "lattice values" << std::endl;
  clout << "DeltaX:                           deltaX="           << converter.getDeltaX()           << std::endl;
  clout << "Lattice velocity:                 latticeU="         << converter.getLatticeU()         << std::endl;
  clout << "DeltaT:                           deltaT="           << converter.getDeltaT()           << std::endl;
  clout << "Reynolds number:                  Re="               << converter.getRe()               << std::endl;
  clout << "DimlessNu:                        dNu="              << converter.getDimlessNu()        << std::endl;
  clout << "Viscosity for computation:        latticeNu="        << converter.getLatticeNu()        << std::endl;
  clout << "Relaxation time:                  tau="              << converter.getTau()              << std::endl;
  clout << "Relaxation frequency:             omega="            << converter.getOmega()            << std::endl;


  clout << "conversion factors" << std::endl;
  clout << "latticeL(m):                      latticeL="         << converter.getLatticeL()         << std::endl;
  clout << "Time step (s):                    physTime="         << converter.physTime()            << std::endl;
  clout << "Velocity factor(m/s):             physVelocity="     << converter.physVelocity()        << std::endl;
  clout << "FlowRate factor(m^d/s):           physFlowRate="     << converter.physFlowRate()        << std::endl;
  clout << "Mass factor(kg):                  physMass="         << converter.physMass()            << std::endl;
  clout << "Force factor(N):                  physForce="        << converter.physForce()           << std::endl;
  clout << "Force factor massless(N/kg):      physMasslessForce="<< converter.physMasslessForce()   << std::endl;
  clout << "Pressure factor(Pa):              physPressure="     << converter.physPressure(4)       << std::endl;
  clout << "latticePressure:                  latticeP="         << converter.latticePressure()     << std::endl;

}

template<typename T>
LBconverter<T>* createLBconverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createLBconverter");

  int dim;
  T latticeL;
  T deltaX;
  int N;
  T latticeU;
  T charNu;
  T Re;
  T charL = T(1);
  T charU = T(1);
  T charRho = T(1);
  T pressureLevel = T(0);

  // params[parameter].read(value) sets the value or returns false if the parameter can not be found

  if (!params["dim"].read<int>(dim)) {
    clout << "Error: Cannot read parameter: dim" << std::endl;
    exit (1);
  }
  if (!params["DiscParam"]["latticeU"].read(latticeU)) {
    clout << "Error: Cannot read parameter: latticeU"
          << std::endl;
    exit (1);
  }
  if (!params["PhysParam"]["charL"].read(charL)) {
    clout << "Parameter charL not found. Set default: charL = 1."
          << std::endl;
  }
  if (!params["DiscParam"]["latticeL"].read(latticeL)) {
    clout << "Parameter latticeL not found. Use deltaX instead."
          << std::endl;
    if (!params["DiscParam"]["deltaX"].read(deltaX)) {
      clout << "Parameter deltaX not found. Use resolution instead."
            << std::endl;
      if (!params["DiscParam"]["resolution"].read(N)) {
        clout << "Error: Cannot read any of the parameters: "
              << "latticeL, deltaX, resolution"
              << std::endl;
        exit (1);
      }
      deltaX = (T) 1 / (T) N;
    }
    latticeL = deltaX * charL;
  }
  if (!params["PhysParam"]["charU"].read(charU)) {
    clout << "Parameter charU not found. Set default: charU = 1."
          << std::endl;
  }
  if (!params["PhysParam"]["charRho"].read(charRho)) {
    clout << "Parameter charRho not found. Set default: charRho = 1."
          << std::endl;
  }
  if (!params["PhysParam"]["charPressure"].read(pressureLevel)) {
    clout << "Parameter charPressure not found. Set default: charPressure = 0."
          << std::endl;
  }
  if (!params["PhysParam"]["charNu"].read(charNu)) {
    if (!params["PhysParam"]["Re"].read(Re)) {
      clout << "Error: Cannot read neither charNu nor Re"
            << std::endl;
      exit (1);
    }
    clout << "Parameter charNu not found. Use Re instead."
          << std::endl;
    charNu = charL*charU / Re;
  }

  return new LBconverter<T>(dim, latticeL, latticeU, charNu, charL, charU,
                            charRho, pressureLevel);
}

}  // namespace olb

#endif
