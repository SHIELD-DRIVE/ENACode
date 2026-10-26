cmake_minimum_required(VERSION 3.16)
project(ENACode CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(ENA_WITH_HDF5 "Build HDF5 input/output (needed for checkpoint dumps)" ON)
option(ENA_WITH_OPENMP "Build with OpenMP threads (milestone 5); results are identical for any thread count" ON)
option(ENA_WITH_MPI "Build with MPI (milestone 6): time-dependent maps over several ranks; results are identical for any layout" ON)

add_library(ena
  src/config.cpp
  src/ascii.cpp
  src/inputs.cpp
  src/streamlines.cpp
  src/regions.cpp
  src/checkpoints.cpp
  src/physics.cpp
  src/stream_processing.cpp
  src/flux.cpp
  src/flux_mi.cpp
  src/flux_esi.cpp
  src/pipeline.cpp
  src/time_dependent.cpp
  src/comm.cpp
)
target_include_directories(ena PUBLIC include)
target_compile_options(ena PRIVATE -Wall -Wextra)
# No fused multiply-add: keeps results reproducible on every CPU.
target_compile_options(ena PUBLIC -ffp-contract=off)
# Intel compilers (icpx, icpc) default to a fast floating-point model that reorders sums and
# uses approximate maths functions; ask for the standard-conforming one, as GCC and Clang use.
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  target_compile_options(ena PUBLIC -fp-model=precise)
endif()

if(ENA_WITH_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(ena PUBLIC OpenMP::OpenMP_CXX)
  else()
    message(WARNING "OpenMP not found: building single-threaded")
  endif()
endif()

if(ENA_WITH_MPI)
  find_package(MPI COMPONENTS CXX)
  if(MPI_CXX_FOUND)
    target_link_libraries(ena PUBLIC MPI::MPI_CXX)
    target_compile_definitions(ena PUBLIC ENA_HAVE_MPI=1)
  else()
    message(WARNING "MPI not found: building without MPI (one rank; threads still work)")
  endif()
endif()

if(ENA_WITH_HDF5)
  find_package(HDF5 COMPONENTS C)
  if(HDF5_FOUND)
    target_sources(ena PRIVATE src/h5io.cpp)
    target_include_directories(ena PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(ena PUBLIC ${HDF5_C_LIBRARIES})
    target_compile_definitions(ena PUBLIC ENA_HAVE_HDF5=1 ${HDF5_DEFINITIONS})
  else()
    message(WARNING "HDF5 not found: checkpoint dumps and the converter are disabled")
  endif()
endif()

add_executable(ena_load apps/ena_load.cpp)
target_link_libraries(ena_load PRIVATE ena)
add_executable(ena_run apps/ena_run.cpp)
target_link_libraries(ena_run PRIVATE ena)
# Step 1 of the Tecplot input preparation (1_fix_neutrals.mcr) in C++: reads the neutral solution
# and writes the fixed neutrals as a Tecplot binary file (see tools/tecplot_prep.py).
add_executable(ena_fixneutrals apps/ena_fixneutrals.cpp)
target_link_libraries(ena_fixneutrals PRIVATE ena)
# The streamlines of the input preparation, traced reproducibly in C++ instead of by Tecplot
# (see tools/tecplot_prep.py and configs/macros/single_ion/3_calculate_files_trace.mcr).
add_executable(ena_trace apps/ena_trace.cpp)
target_link_libraries(ena_trace PRIVATE ena)

if(HDF5_FOUND)
  add_executable(ena_convert apps/ena_convert.cpp)
  target_link_libraries(ena_convert PRIVATE ena)
  add_executable(test_physics tests/test_physics.cpp)
  target_link_libraries(test_physics PRIVATE ena)
  add_executable(test_physics_mi tests/test_physics_mi.cpp)
  target_link_libraries(test_physics_mi PRIVATE ena)
endif()

enable_testing()
add_executable(test_m2 tests/test_m2.cpp)
target_link_libraries(test_m2 PRIVATE ena)
add_test(NAME milestone2 COMMAND test_m2 ${CMAKE_CURRENT_BINARY_DIR}/test_tmp)
