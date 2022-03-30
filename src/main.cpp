#include <assert.h>
#include <iomanip>
#include <iostream>
#include <vector>
//#ifdef USEMPI
#include <mpi.h>
//#endif
#include "arguments.hpp"
#include "solver.hpp"
 
int main(int argc, char *argv[]) {

  int rank=0;
  int numproc=1;
  int ndims = 1;
//#ifdef USEMPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//#endif

  // parse command line arguments

  auto resolution = convertTo<int>(2, 32, argc, argv);
  auto iterations = convertTo<int>(3, 1000, argc, argv);
  std::string str_ndims = argv[1];
  
  if (str_ndims == "2D") {
    ndims = 1;
  };
  
  // std::cout << "numproc=" << numproc << std::endl;
  // std::cout << "resolution=" << resolution << std::endl;
  // std::cout << "iterations=" << iterations << std::endl;

  assert(resolution > 0);
  assert(iterations > 0);

  solve(resolution, iterations,rank,numproc, ndims);

//#ifdef USEMPI
  MPI_Finalize();
//#endif

  return 0;
}
