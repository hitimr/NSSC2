#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <assert.h>


#ifdef USEMPI
	#include <mpi.h>
#endif

#include "logging.hpp"
#include "splitting.hpp"

// globals according to command line
extern int g_n_processes;
extern int g_dim;
extern size_t g_resolution;
extern int g_iterations;
extern int g_my_rank;


template <typename Type> class MatrixView 
{
private:	
	MatrixView(const MatrixView &);
	std::vector<Type> &v;
	MatrixView &operator=(const MatrixView &);

public:
	const size_t N, M;
	MatrixView(std::vector<Type> &v, size_t N, size_t M) : v(v), N(N), M(M) 
	{
		assert(v.size() / N == M);
	}
	Type &set(size_t i, size_t j) { return v[i + N * j]; }
	const Type &get(size_t i, size_t j) { return v[i + N * j]; }
	Type &set(size_t n) { return v[n]; }
	const Type &get(size_t n) { return v[n]; }
};

FP_TYPE ParticularSolution(FP_TYPE x, FP_TYPE y) 
{
	return sin(2 * M_PI * x) * sinh(2 * M_PI * y);
}

FP_TYPE NormL2(const std::vector<FP_TYPE> &v) 
{
	FP_TYPE norm = 0;
	for (const auto &value : v) {
		norm += value * value;
	}
	return sqrt(norm);
}	

FP_TYPE NormInf(const std::vector<FP_TYPE> &v) 
{
	FP_TYPE max = std::numeric_limits<FP_TYPE>::lowest();
	for (const auto &value : v) {
		max = std::fabs(value) > max ? std::fabs(value) : max;
	}
	return max;
}

struct Stencil 
{
	Stencil(FP_TYPE h)
			: C(4.0 / (h * h) + 4 * M_PI * M_PI), N(-1.0 / (h * h)),
				S(-1.0 / (h * h)), W(-1.0 / (h * h)), E(-1.0 / (h * h)) {}
	const FP_TYPE C, N, S, W, E;
};

enum Cell { UNKNOWN = 0, DIR = 1, NEU = 2, ROB = 0 };

void solve(size_t resolution, size_t iterations, int mpi_rank) 
{
	Logger log(mpi_rank);
	

	#ifdef USEMPI
		MPI_Barrier(MPI_COMM_WORLD);	
		std::cout << "Rank=" << mpi_rank << std::endl;
		
		// Create new communicator
		int ndims = 2;
		int dims[2] = {};
		int bcs[2] = {0,0};	 // logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
		int reorder = 1;			// ranking may be reordered (true) or not (false) (logical)

		MPI_Comm topo_com;
		MPI_Dims_create(g_n_processes, ndims, dims);	// https://www.mpich.org/static/docs/v3.3.x/www3/MPI_Dims_create.html 
		MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &topo_com);	// https://www.mpich.org/static/docs/v3.3/www3/MPI_Cart_create.html
		MPI_Barrier(MPI_COMM_WORLD);	

		int coords[2] = {};								 
		MPI_Cart_coords(topo_com, mpi_rank, ndims, coords);
		std::cout << "dims=(" << dims[0] << ", " << dims[1] << ")" << std::endl;
		std::cout << "coords=(" << coords[0] << ", " << coords[1] << ")" << std::endl;
	#endif	

		size_t NY = resolution;
		size_t NX = (2.0 * NY) - 1;
		FP_TYPE h = 1.0 / (NY - 1);

		auto grid_size = local_grid_size(g_my_rank);
		NX = grid_size[COORD_X];
		NY = grid_size[COORD_Y];



		const auto stencil = Stencil(h);


		// domain cell types
		std::vector<int> domain(NX * NY, Cell::UNKNOWN);
		MatrixView<int> domainView(domain, NX, NY);
		for (size_t i = 0; i != NX; ++i) 
		{
			domainView.set(i, 0) = Cell::DIR;
			domainView.set(i, NY - 1) = Cell::DIR;
		}
		for (size_t j = 0; j != NY; ++j) 
		{
			domainView.set(0, j) = Cell::DIR;
			domainView.set(NX - 1, j) = Cell::DIR;
		}

		// referenceSolution
		std::vector<FP_TYPE> referenceSolution(NX * NY, 0);
		MatrixView<FP_TYPE> referenceSolutionView(referenceSolution, NX, NY);
		for (size_t j = 0; j != NY; ++j) 
		{
			for (size_t i = 0; i != NX; ++i) 
			{
				referenceSolutionView.set(i, j) = ParticularSolution(i * h, j * h);
			}
		}

		// right hand side
		std::vector<FP_TYPE> rightHandSide(NX * NY, 0);
		MatrixView<FP_TYPE> rightHandSideView(rightHandSide, NX, NY);
		for (size_t j = 0; j != NY; ++j) 
		{
			for (size_t i = 0; i != NX; ++i) 
			{
				rightHandSideView.set(i, j) =
						ParticularSolution(i * h, j * h) * 4 * M_PI * M_PI;
			}
		}

		auto SolverJacobi = [](std::vector<FP_TYPE> &sol, std::vector<FP_TYPE> &sol2,
													std::vector<FP_TYPE> &rhs, const Stencil &stencil,
													size_t NX, size_t NY) {
			MatrixView<FP_TYPE> solView(sol, NX, NY);
			MatrixView<FP_TYPE> sol2View(sol2, NX, NY);
			MatrixView<FP_TYPE> rhsView(rhs, NX, NY);

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

		auto ComputeResidual = [](std::vector<FP_TYPE> &sol, std::vector<FP_TYPE> &rhs,
															const Stencil &stencil, size_t NX, size_t NY) {
			MatrixView<FP_TYPE> solView(sol, NX, NY);
			MatrixView<FP_TYPE> rhsView(rhs, NX, NY);

			std::vector<FP_TYPE> residual(NX * NY, 0);
			MatrixView<FP_TYPE> residualView(residual, NX, NY);
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
		auto ComputeError = [](std::vector<FP_TYPE> &sol,
													std::vector<FP_TYPE> &reference, size_t NX, size_t NY) {
			MatrixView<FP_TYPE> solView(sol, NX, NY);
			MatrixView<FP_TYPE> referenceView(reference, NX, NY);

			std::vector<FP_TYPE> error(NX * NY, 0);
			MatrixView<FP_TYPE> errorView(error, NX, NY);

			for (size_t j = 1; j != NY - 1; ++j) {
				for (size_t i = 1; i != NX - 1; ++i) {
					errorView.set(i, j) = referenceView.get(i, j) - solView.get(i, j);
				}
			}
			return error;
		};

		// solution approximation starting with boundary initialized to dirichlet
		// conditions, else 0
		std::vector<FP_TYPE> solution(NX * NY, 0);
		MatrixView<FP_TYPE> solutionView(solution, NX, NY);
		for (size_t j = 0; j != NY; ++j) {
			for (size_t i = 0; i != NX; ++i) {
				if (domainView.get(i, j) == Cell::DIR)
					solutionView.set(i, j) = ParticularSolution(i * h, j * h);
			}
		}
		std::vector<FP_TYPE> solution2 = solution;
		std::cout << "solve LSE using stencil jacobi" << std::endl;
		auto start = std::chrono::high_resolution_clock::now();
		for (size_t iter = 0; iter <= iterations; ++iter) {
			SolverJacobi(solution, solution2, rightHandSide, stencil, NX, NY);
		}

		auto stop = std::chrono::high_resolution_clock::now();
		auto seconds =
				std::chrono::duration_cast<std::chrono::duration<FP_TYPE>>(stop - start)
						.count();
		std::cout << std::scientific << "runtime=" << seconds << std::endl;


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

		log.add("n_processes", std::to_string(g_n_processes));
		log.add("runtime", std::to_string(seconds));
		log.add("error", std::to_string(errorNorm));
}
