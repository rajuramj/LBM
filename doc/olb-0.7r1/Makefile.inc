#  This file is part of the OpenLB library
#
#  Copyright (C) 2007 Mathias Krause
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public 
#  License along with this program; if not, write to the Free 
#  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.

###########################################################################
###########################################################################
## DEFINITIONS TO BE CHANGED

CXX             := g++
#CXX             := icpc -D__aligned__=ignored
#CXX             := mpiCC

OPTIM           := -O3 -Wall
DEBUG           := -g -DOLB_DEBUG

CXXFLAGS        := $(OPTIM)
#CXXFLAGS        := $(DEBUG)

LDFLAGS         :=

PARALLEL_MODE   := OFF
#PARALLEL_MODE   := MPI
#PARALLEL_MODE   := OMP
#PARALLEL_MODE   := HYBRID

MPIFLAGS        :=
OMPFLAGS        := -fopenmp


BUILDTYPE       := precompiled
#BUILDTYPE       := generic


###########################################################################
## conditional settings

ifeq ($(BUILDTYPE), precompiled)
   CXXFLAGS := -DOLB_PRECOMPILED $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), MPI)
   CXXFLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS) $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), OMP)
   CXXFLAGS := -DPARALLEL_MODE_OMP $(OMPFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), HYBRID)
   CXXFLAGS := -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

###########################################################################
## defines shell

SHELL           := /bin/sh

###########################################################################
## dependencies, object, library directory and library name

DEPENDDIR       := build/$(BUILDTYPE)/dep
OBJDIR          := build/$(BUILDTYPE)/obj
LIBDIR          := build/$(BUILDTYPE)/lib
LIB             := olb

###########################################################################
## search directories

SUBDIRS         := src/core \
                   src/boundaryConditions \
                   src/complexDynamics \
                   src/complexGrids/mpiManager \
                   src/complexGrids/cuboidStructure \
                   src/complexGrids/multiBlockStructure \
                   src/utilities \
                   src/io \
                   src/external/tinyxml

EXAMPLEDIRS     := examples/bstep2d \
                   examples/bstep3d \
                   examples/cavity2d \
                   examples/cavity3d \
                   examples/cavity3d_cuboid_mpi \
                   examples/cylinder2d \
                   examples/cylinder3d \
                   examples/forcedPoiseuille2d \
                   examples/mrt2d \
                   examples/multiComponent2d \
                   examples/multiComponent3d \
                   examples/poiseuille2d \
                   examples/thermal2d \
                   examples/thermal3d

INCLUDEDIRS     := src \
		   src/external \
                   src/external/cvmlcpp \
		   src/external/boost \
		   src/external/kaffpa
			

BUILDTYPEDIRS   := build/precompiled \
                   build/generic

SUBDIRSH        := $(foreach d,$(SUBDIRS),$(ROOT)/$(d))
IDIR            := $(foreach d,$(INCLUDEDIRS),-I$(ROOT)/$(d))


###########################################################################
###########################################################################
