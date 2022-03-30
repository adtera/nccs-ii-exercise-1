#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
//#ifdef USEMPI
#include <mpi.h>
//#endif
#include <assert.h>

template <typename Type> class MatrixView {
private:
  std::vector<Type> &v;
  MatrixView(const MatrixView &);
  MatrixView &operator=(const MatrixView &);

public:
  const size_t N, M;
  MatrixView(std::vector<Type> &v, size_t N, size_t M) : v(v), N(N), M(M) {
    assert(v.size() / N == M);
  }
  Type &set(size_t i, size_t j) { return v[i + N * j]; }
  const Type &get(size_t i, size_t j) { return v[i + N * j]; }
  Type &set(size_t n) { return v[n]; }
  const Type &get(size_t n) { return v[n]; }
};

double ParticularSolution(double x, double y) {
  return sin(2 * M_PI * x) * sinh(2 * M_PI * y);
}

double NormL2(const std::vector<double> &v) {
  double norm = 0;
  for (const auto &value : v) {
    norm += value * value;
  }
  return sqrt(norm);
}

double NormInf(const std::vector<double> &v) {
  double max = std::numeric_limits<double>::lowest();
  for (const auto &value : v) {
    max = std::fabs(value) > max ? std::fabs(value) : max;
  }
  return max;
}

struct Stencil {
  Stencil(double h)
      : C(4.0 / (h * h) + 4 * M_PI * M_PI), N(-1.0 / (h * h)),
        S(-1.0 / (h * h)), W(-1.0 / (h * h)), E(-1.0 / (h * h)) {}
  const double C, N, S, W, E;
};

enum Cell { UNKNOWN = 0, DIR = 1, DOWN = 2, UP = 3, LEFT = 4, RIGHT = 5 };

void solve(size_t resolution, size_t iterations, int mpi_rank,
           int mpi_numproc, int ndims) {

// Solver Template stops here, MPI subdomain solve starts 


//  #ifdef USEMPI


  int myrank,numprocs;
  // double h = 1.0 / (NY - 1);
  double h = 1.0 / (resolution - 1);
  
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

  // Build Cartesian Grid
//  std::cout << "Init Processor Rank " << myrank << std::endl;
  int dims[2] = {0,0};
  if (ndims == 1) {
    dims[1] = 1;
//    std::cout << "This is dims " << dims << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Dims_create(numprocs, ndims, dims);
  MPI_Comm GRID_COMM;  
  int bcs[2] = {0,0};
  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &GRID_COMM);
  int coords[2] = {};
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Cart_coords(GRID_COMM, myrank, ndims, coords);
    


  int width_x, width_y, last_width_x, last_width_y;
  int modulo_x, modulo_y;
  size_t NY=0; 
  size_t NX=0;
  bool is_subdomain = coords[0] + 1 != dims[0] && coords[1] + 1 != dims[1];
  bool is_x_last_subdomain = coords[0] + 1 == dims[0] && coords[1] + 1 != dims[1];
  bool is_y_last_subdomain = coords[0] + 1 != dims[0] && coords[1] + 1 == dims[1];
  bool is_xy_last_subdomain = coords[0] + 1 == dims[0] && coords[1] + 1 == dims[1];

  width_x = (resolution - 2)/dims[0];
  width_y = (resolution - 2)/dims[1];
    
  modulo_x = (resolution - 2) % dims[0];
  modulo_y = (resolution - 2) % dims[1];
  
  last_width_x = width_x + modulo_x;
  last_width_y = width_y + modulo_y;
  

  if (myrank == 0) {
 	  std::cout << "numprocs = " << numprocs << std::endl;
    std::cout << "dims=(" << dims[0] << "," << dims[1] <<  ")" << std::endl;
    std::cout << "coords=(" << coords[0] << "," << coords[1] << ")" << std::endl;
  }
  

  if (is_subdomain) {
    NX = width_x+2;
    NY = width_y+2;
  }
   else if (is_x_last_subdomain) {
    NX = last_width_x+2;
    NY = width_y+2;
  }
   else if (is_y_last_subdomain) {
    NX = width_x+2;
    NY = last_width_y+2;
  }
   else if (is_xy_last_subdomain) {
    NX = last_width_x+2;
    NY = last_width_y+2;
  }


  std::vector<int> domain(NX * NY, Cell::UNKNOWN);
  MatrixView<int> domainView(domain, NX, NY);
  for (size_t i = 1; i != NX - 1; ++i) {
    domainView.set(i, 0) = Cell::DOWN;
    if (coords[1] + 1 == dims[1]) {
      domainView.set(i,NY - 1) = Cell::DIR;
    } else { 
      domainView.set(i, NY - 1) = Cell::UP;
    };
  }
  for (size_t j = 1; j != NY - 1; ++j) {
    domainView.set(0, j) = Cell::LEFT;
    domainView.set(NX - 1, j) = Cell::RIGHT;
  }

//  std::cout << "Init Done Processor Rank " << myrank << std::endl;

	MPI_Barrier(GRID_COMM);
	std::cout << "coords: " << coords[0] << "," << coords[1] << " rank: " << myrank << " size: " << NX << "," << NY << std::endl;

// each worker has its subdomain incl. ghost layers

  const auto stencil = Stencil(h);

  // referenceSolution, needed to compute error 
  std::vector<double> referenceSolution(NX * NY, 0);
  MatrixView<double> referenceSolutionView(referenceSolution, NX, NY);
  for (size_t j = 0; j != NY; ++j) {
    for (size_t i = 0; i != NX; ++i) {
      referenceSolutionView.set(i, j) = ParticularSolution(((NX-2)*coords[0] + i) * h, ((NY-2)*coords[1] + j) * h);
    }
  }


  // right hand side for subdomain
  std::vector<double> rightHandSide(NX * NY, 0);
  MatrixView<double> rightHandSideView(rightHandSide, NX, NY);
  for (size_t j = 0; j != NY; ++j) {
    for (size_t i = 0; i != NX; ++i) {
      rightHandSideView.set(i, j) =
          ParticularSolution(((NX-2)*coords[0] + i) * h,((NY-2)*coords[1] + j) * h) * 4 * M_PI * M_PI;
    }
  }



  auto SolverJacobi = [](std::vector<double> &sol, std::vector<double> &sol2,
                         std::vector<double> &rhs, const Stencil &stencil,
                         size_t NX, size_t NY) {
    MatrixView<double> solView(sol, NX, NY);
    MatrixView<double> sol2View(sol2, NX, NY);
    MatrixView<double> rhsView(rhs, NX, NY);

    for (size_t j = 1; j != NY - 1; ++j) {
      for (size_t i = 1; i != NX - 1; ++i) {
        sol2View.set(i, j) =
            1.0 / stencil.C *
            (rhsView.set(i, j) - (solView.get(i + 1, j) * stencil.E +
                                  solView.get(i - 1, j) * stencil.W +
                                  solView.get(i, j + 1) * stencil.S +
                                  solView.get(i, j - 1) * stencil.N));
      }
    }
    sol.swap(sol2);
  };

  auto ComputeResidual = [](std::vector<double> &sol, std::vector<double> &rhs,
                            const Stencil &stencil, size_t resolutionX, size_t resolutionY) {
    MatrixView<double> solView(sol, resolutionX, resolutionY);
    MatrixView<double> rhsView(rhs, resolutionX, resolutionY);

    std::vector<double> residual(resolutionX * resolutionY, 0);
    MatrixView<double> residualView(residual, resolutionX, resolutionY);
    for (size_t j = 1; j != resolutionY - 1; ++j) {
      for (size_t i = 1; i != resolutionX - 1; ++i) {
        residualView.set(i, j) =
            rhsView.get(i, j) -
            (solView.get(i, j) * stencil.C + solView.get(i + 1, j) * stencil.E +
             solView.get(i - 1, j) * stencil.W +
             solView.get(i, j - 1) * stencil.S +
             solView.get(i, j + 1) * stencil.N);
      }
    }
    return residual;
  };
  auto ComputeError = [](std::vector<double> &sol,
                         std::vector<double> &reference, size_t NX, size_t NY) {
    MatrixView<double> solView(sol, NX, NY);
    MatrixView<double> referenceView(reference, NX, NY);

    std::vector<double> error(NX * NY, 0);
    MatrixView<double> errorView(error, NX, NY);

    for (size_t j = 1; j != NY - 1; ++j) {
      for (size_t i = 1; i != NX - 1; ++i) {
        errorView.set(i, j) = referenceView.get(i, j) - solView.get(i, j);
      }
    }
    return error;
  
  };



  // solution approximation starting with boundary initialized to dirichlet
  // conditions, else 0
  
  std::vector<double> solution(NX * NY, 0);
  MatrixView<double> solutionView(solution, NX, NY);
  for (size_t j = 0; j != NY; ++j) {
    for (size_t i = 0; i != NX; ++i) {
      if (domainView.get(i, j) == Cell::DIR) {
        solutionView.set(i, j) = ParticularSolution(((NX-2)*coords[0] + i) * h, ((NY-2)*coords[1] + j) * h);
      }
    }
  };

  std::vector <double> solution2 = solution;
  std::vector <double> UP_SEND(NX, 0);
  std::vector <double> DOWN_SEND(NX, 0);
  std::vector <double> LEFT_SEND(NY, 0);
  std::vector <double> RIGHT_SEND(NY, 0);
  std::vector <double> UP_RECV(NX, 0);
  std::vector <double> DOWN_RECV(NX, 0);
  std::vector <double> LEFT_RECV(NY, 0);
  std::vector <double> RIGHT_RECV(NY, 0);

  int up_rank;
  int down_rank;
  int left_rank; 
  int right_rank;

  MPI_Cart_shift(GRID_COMM, 0, 1, &left_rank, &right_rank);

  MPI_Cart_shift(GRID_COMM, 1, 1, &down_rank, &up_rank);


//  std::cout << "solve LSE using stencil jacobi" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  for (size_t iter = 0; iter <= iterations; ++iter) {
    SolverJacobi(solution, solution2, rightHandSide, stencil, NX, NY);
    
    for (size_t j = 0; j != NY; ++j) {
      for (size_t i = 0; i != NX; ++i) {
        if (domainView.get(i, j) == Cell::UP && dims[1] != 1) {
          UP_SEND[i] = solutionView.get(i, j-1);
        } else if (domainView.get(i, j) == Cell::DOWN && dims[1] != 1){
          DOWN_SEND[i] = solutionView.get(i, j+1);
        } else if (domainView.get(i, j) == Cell::LEFT) {
          LEFT_SEND[j] = solutionView.get(i+1,j);
        } else if (domainView.get(i, j) == Cell::RIGHT) {
          RIGHT_SEND[j] = solutionView.get(i-1,j);
        };
      }
    }

    if (dims[1] != 1) {
      MPI_Send(UP_SEND.data(), NX, MPI_DOUBLE, up_rank, 1, GRID_COMM);
      MPI_Recv(DOWN_RECV.data(), NX, MPI_DOUBLE, down_rank, 1, GRID_COMM, MPI_STATUS_IGNORE);
    
      MPI_Send(DOWN_SEND.data(), NX, MPI_DOUBLE, down_rank, 2, GRID_COMM);
      MPI_Recv(UP_RECV.data(), NX, MPI_DOUBLE, up_rank, 2, GRID_COMM, MPI_STATUS_IGNORE);
    };
    MPI_Send(LEFT_SEND.data(), NY, MPI_DOUBLE, left_rank, 3, GRID_COMM);
    MPI_Recv(RIGHT_RECV.data(), NY, MPI_DOUBLE, right_rank, 3, GRID_COMM, MPI_STATUS_IGNORE);
    
    MPI_Send(RIGHT_SEND.data(), NY, MPI_DOUBLE, right_rank, 4, GRID_COMM);
    MPI_Recv(LEFT_RECV.data(), NY, MPI_DOUBLE, left_rank, 4, GRID_COMM, MPI_STATUS_IGNORE);


    for (size_t j = 0; j != NY; ++j) {
      for (size_t i = 0; i != NX; ++i) {
        if (domainView.get(i, j) == Cell::UP && dims[1] != 1){
          solutionView.set(i, j) = UP_RECV[i];
        } else if (domainView.get(i, j) == Cell::DOWN && dims[1] != 1){
          solutionView.set(i, j) = DOWN_RECV[i];
        } else if (domainView.get(i,j) == Cell::LEFT) {
          solutionView.set(i,j) = LEFT_RECV[j];
        } else if (domainView.get(i,j) == Cell::RIGHT) {
          solutionView.set(i,j) = RIGHT_RECV[j];
        };
      }
    }
  }


{
  auto stop = std::chrono::high_resolution_clock::now();
  auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
  auto residual = ComputeResidual(solution, rightHandSide, stencil, NX, NY);
  auto residualNorm = NormL2(residual);
  auto residualMax = NormInf(residual);
  auto error = ComputeError(solution, referenceSolution, NX, NY);
  auto errorNorm = NormL2(error);
  auto errorMax = NormInf(error);

  double seconds_sum;
  double residualNorm_sum;
  double residualMax_max;
  double errorNorm_sum;
  double errorMax_max;

  MPI_Reduce(&seconds, &seconds_sum, 1, MPI_DOUBLE, MPI_SUM, 0, GRID_COMM);
  MPI_Reduce(&residualNorm, &residualNorm_sum, 1, MPI_DOUBLE, MPI_SUM, 0, GRID_COMM);
  MPI_Reduce(&residualMax, &residualMax_max, 1, MPI_DOUBLE, MPI_MAX, 0, GRID_COMM);
  MPI_Reduce(&errorNorm, &errorNorm_sum, 1, MPI_DOUBLE, MPI_SUM, 0, GRID_COMM);
  MPI_Reduce(&errorMax, &errorMax_max, 1, MPI_DOUBLE, MPI_MAX, 0, GRID_COMM);

  if (myrank == 0) {
  std::cout << std::scientific << "|total runtime|= " << seconds << " seconds" << std::endl;
  std::cout << std::scientific << "|average runtime per process|= " << seconds/numprocs << " seconds per processor" << std::endl;
  std::cout << std::scientific << "|residual|=" << residualNorm << std::endl;
  std::cout << std::scientific << "|residualMax|=" << residualMax << std::endl;
  std::cout << std::scientific << "|error|=" << errorNorm << std::endl;
  std::cout << std::scientific << "|errorMax|=" << errorMax << std::endl;
  }
//  #endif
};

}

