cmake_minimum_required(VERSION 3.6.2)
project(doutzen)

#find_package(MPI REQUIRED)

#Bring the headers, such as Student.h into the project
#include_directories(${MPI_INCLUDE_PATH})
#include_directories(include /User/aaram/local/include /usr/lib/openmpi/include)
set(EIGEN_INCLUDE_PATH $ENV{HOME}/local/include)
include_directories(include ${EIGEN_INCLUDE_PATH} ${MPI_INCLUDE_PATH})
set(FFTW_LIBRARY_PATH /opt/local/lib)
link_directories(${FFTW_LIBRARY_PATH})

#Can manually add the sources using the set command as follows:
#set(SOURCES src/main.cpp src/parameter.cpp src/lattice.cpp src/basis.cpp)

#However, the file(GLOB...) allows for wildcard additions:

#SET(CMAKE_CXX_COMPILER "mpic++.openmpi")
#add_definitions(${CMAKE_CXX_COMPILER},${GCC_COVERAGE_COMPILE_FLAGS})
SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -O3 -std=c++17")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

#file(GLOB SOURCES src/main.cpp src/diagram.cpp src/element.cpp src/graph.cpp src/green.cpp src/lattice.cpp src/parameters.cpp)
#add_executable(doutzen ${SOURCES})
#target_link_libraries(doutzen ${MPI_LIBRARIES} fftw3)

file(GLOB SOURCES src/exactCorrelator.cpp)
add_executable(RSJJ_exact ${SOURCES})

file(GLOB SOURCES src/selfconsistency_NCA.cpp src/dyson.cpp)
add_executable(RSJJ_NCA ${SOURCES})

file(GLOB SOURCES src/selfconsistency_OCA.cpp src/dyson.cpp)
add_executable(RSJJ_OCA ${SOURCES})

file(GLOB SOURCES src/selfconsistency_TOA.cpp src/dyson.cpp)
add_executable(RSJJ_TOA ${SOURCES})

#file(GLOB SOURCES src/selfconsistency_Ladder.cpp src/dyson.cpp)
#add_executable(doutzen_Ladder ${SOURCES})
#
#file(GLOB SOURCES src/selfconsistency_TCA.cpp src/dyson.cpp)
#add_executable(doutzen_TCA ${SOURCES})
#
#file(GLOB SOURCES src/selfconsistency_Q2.cpp src/dyson.cpp)
#add_executable(doutzen_Q2 ${SOURCES})
#
#file(GLOB SOURCES src/diagMCvertex.cpp src/dyson.cpp src/diagram.cpp src/configuration.cpp src/element.cpp src/green.cpp src/parameters.cpp src/graph.cpp)
#add_executable(doutzen_diagMC ${SOURCES})
#target_link_libraries(doutzen_diagMC ${MPI_LIBRARIES} fftw3)
#
#file(GLOB SOURCES src/singleModeWaveguide_diagMC.cpp src/dyson.cpp src/diagram.cpp src/configuration.cpp src/element.cpp src/green.cpp src/parameters.cpp src/graph.cpp)
#add_executable(doutzen_smw_diagMC ${SOURCES})
#target_link_libraries(doutzen_smw_diagMC ${MPI_LIBRARIES} fftw3)
