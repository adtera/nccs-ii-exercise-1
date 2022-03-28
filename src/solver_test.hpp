#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#ifdef USEMPI
#include <mpi.h>
#endif
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

enum Cell { UNKNOWN = 0, DIR = 1, NEU = 2, ROB = 0 };

void solve(size_t resolution, size_t iterations, int mpi_rank,
           int mpi_numproc) {

#ifdef USEMPI
  std::cout << "Rank=" << mpi_rank << std::endl;
  int ndims = 2;
  int dims[2] = {};
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Dims_create(mpi_numproc, ndims, dims);
  std::cout << "dims=(" << dims[0] << ", " <  < dims[1] << ")" << std::endl;
  MPI_Comm topo_com;
  int bcs[2] = {0, 0};
  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &topo_com);
  int coords[2] = {};
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Cart_coords(topo_com, mpi_rank, ndims, coords);
  std::cout << "coords=(" << coords[0] << ", " << coords[1] << ")" << std::endl;
#endif

  size_t NY = resolution;
  size_t NX = resolution;
  double h = 1.0 / (NY - 1);

  const auto stencil = Stencil(h);

  // domain cell types
  std::vector<int> domain(NX * NY, Cell::UNKNOWN);
  MatrixView<int> domainView(domain, NX, NY);
  for (size_t i = 0; i != NX; ++i) {
    domainView.set(i, 0) = Cell::DIR;
    domainView.set(i, NY - 1) = Cell::DIR;
  }
  for (size_t j = 0; j != NY; ++j) {
    domainView.set(0, j) = Cell::DIR;
    domainView.set(NX - 1, j) = Cell::DIR;
  }

  // referenceSolution
  std::vector<double> referenceSolution(NX * NY, 0);
  MatrixView<double> referenceSolutionView(referenceSolution, NX, NY);
  for (size_t j = 0; j != NY; ++j) {
    for (size_t i = 0; i != NX; ++i) {
      referenceSolutionView.set(i, j) = ParticularSolution(i * h, j * h);
    }
  }

  // right hand side
  std::vector<double> rightHandSide(NX * NY, 0);
  MatrixView<double> rightHandSideView(rightHandSide, NX, NY);
  for (size_t j = 0; j != NY; ++j) {
    for (size_t i = 0; i != NX; ++i) {
      rightHandSideView.set(i, j) =
          ParticularSolution(i * h, j * h) * 4 * M_PI * M_PI;
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
                            const Stencil &stencil, size_t NX, size_t NY) {
    MatrixView<double> solView(sol, NX, NY);
    MatrixView<double> rhsView(rhs, NX, NY);

    std::vector<double> residual(NX * NY, 0);
    MatrixView<double> residualView(residual, NX, NY);
    for (size_t j = 1; j != NY - 1; ++j) {
      for (size_t i = 1; i != NX - 1; ++i) {
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
      if (domainView.get(i, j) == Cell::DIR)
        // particular soolution for each subdomain
        solutionView.set(i, j) = ParticularSolution(i * h, j * h);
    }
  } 
  std::vector<double> solution2 = solution;
  std::cout << "solve LSE using stencil jacobi" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  for (size_t iter = 0; iter <= iterations; ++iter) {
      // exchange ghost layers between workers
    SolverJacobi(solution, solution2, rightHandSide, stencil, NX, NY);
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto seconds =
      std::chrono::duration_cast<std::chrono::duration<double>>(stop - start)
          .count();
  std::cout << std::scientific << "runtime=" << seconds << std::endl;

  {
    auto residual = ComputeResidual(solution, rightHandSide, stencil, NX, NY);
    auto residualNorm = NormL2(residual);
    std::cout << std::scientific << "|residual|=" << residualNorm << std::endl;
    auto residualMax = NormInf(residual);
    std::cout << std::scientific << "|residualMax|=" << residualMax
              << std::endl;
    auto error = ComputeError(solution, referenceSolution, NX, NY);
    auto errorNorm = NormL2(error);
    std::cout << std::scientific << "|error|=" << errorNorm << std::endl;
    auto errorMax = NormInf(error);
    std::cout << std::scientific << "|errorMax|=" << errorMax << std::endl;
  }
}
