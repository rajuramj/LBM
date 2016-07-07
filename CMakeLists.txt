cmake_minimum_required(VERSION 2.8.9)

# Defining the project name
project (SiWiR2_LBM)

# setting the compiler flags 
#set(CMAKE_CXX_FLAGS "-Wall -Winline -Wshadow -pedantic")
set(CMAKE_CXX_FLAGS "-Wall -Winline -pedantic")

add_definitions(-std=c++11)

# Bringing in the include directories
include_directories(include)

#Adding the sources using the set command
set(SOURCES test/main.cpp src/Lattice.cpp src/Parameters.cpp src/Simulation.cpp) 

#Setting the executable file
add_executable(lbm ${SOURCES})



