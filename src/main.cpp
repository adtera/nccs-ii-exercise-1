#include <assert.h>
#include <iomanip>
#include <iostream>
#include <vector>
#ifdef USEMPI
#include <mpi.h>
#endif
#include "arguments.hpp"
#include "solver.hpp"
 
int main(int argc, char *argv[]) {

  int rank=0;
  int numproc=1;
#ifdef USEMPI
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  // parse command line arguments
  auto resolution = convertTo<int>(1, 32, argc, argv);
  auto iterations = convertTo<int>(2, 1000, argc, argv);

  std::cout << "numproc=" << numproc << std::endl;
  std::cout << "resolution=" << resolution << std::endl;
  std::cout << "iterations=" << iterations << std::endl;

  assert(resolution > 0);
  assert(iterations > 0);

  solve(resolution, iterations,rank,numproc);

#ifdef USEMPI
  MPI_Finalize();
#endif

  return 0;
}
