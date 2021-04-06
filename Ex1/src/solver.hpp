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
	extern MPI_Comm g_topo_com;
#endif

#include "logging.hpp"
#include "splitting.hpp"
#include "common.hpp"



template <typename Type> class MatrixView 
{	
public:
	MatrixView(const MatrixView &);
	std::vector<Type> &v;
	MatrixView &operator=(const MatrixView &);

	const size_t N, M;
	MatrixView(std::vector<Type> &v, size_t N, size_t M) : v(v), N(N), M(M) 
	{
		assert(v.size() / N == M);
	}
	Type &set(size_t i, size_t j) { return v[i + N * j]; }
	Type &get(size_t i, size_t j) { return v[i + N * j]; }
	Type &set(size_t n) { return v[n]; }
	Type &get(size_t n) { return v[n]; }
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

void solve(size_t resolution, size_t iterations) 
{
	//Logger log(g_my_rank);
	

	#ifdef USEMPI
		MPI_Barrier(MPI_COMM_WORLD);	
		std::cout << "Rank=" << g_my_rank << std::endl;
		
		// Create new communicator
		int ndims = 2;
		int dims[2] = {};
		int bcs[2] = {0,0};	 // logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
		int reorder = 1;			// ranking may be reordered (true) or not (false) (logical)

		MPI_Dims_create(g_n_processes, ndims, dims);	// https://www.mpich.org/static/docs/v3.3.x/www3/MPI_Dims_create.html 
		MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &g_topo_com);	// https://www.mpich.org/static/docs/v3.3/www3/MPI_Cart_create.html
		MPI_Barrier(MPI_COMM_WORLD);	

		int coords[2] = {};								 
		MPI_Cart_coords(g_topo_com, g_my_rank, ndims, coords);
		std::cout << "dims=(" << dims[0] << ", " << dims[1] << ")" << std::endl;
		std::cout << "coords=(" << coords[0] << ", " << coords[1] << ")" << std::endl;


		// calculate local grid size
		vector<int> my_coords;
		if(g_dim == DIM1)
		{
			my_coords = {coords[0]};
		} 
		else 
		{
			my_coords = {coords[0], coords[1]};	// transform to vector
		}
		
		auto grid_size = local_grid_size(my_coords, true);
		size_t NX = grid_size[COORD_X];
		size_t NY = grid_size[COORD_Y];
		FP_TYPE h = 1.0 / (NY - 1);	
	#else
		// Serial
		size_t NY = resolution;
		size_t NX = (2.0 * NY) - 1;
		FP_TYPE h = 1.0 / (NY - 1);
	#endif	


		const auto stencil = Stencil(h);


		// domain cell types
		std::vector<int> domain(NX * NY, Cell::UNKNOWN);
		MatrixView<int> domainView(domain, NX, NY);

		int topProc;
		int botProc;
		MPI_Cart_shift(g_topo_com, 1, 1, &topProc, &botProc);
		auto borders = border_types(my_coords);

		for (size_t i = 0; i != NX; ++i) 
		{
			domainView.set(i, 0) = Cell::DIR;	// left
			domainView.set(i, NY - 1) = Cell::DIR;	// right
		}

		// top/bot
		for (size_t j = 0; j != NY; ++j) 
		{
			if(borders[TOP] == BORDER_DIR)
				domainView.set(0, j) = Cell::DIR;

			if(borders[BOTTOM] == BORDER_DIR)
				domainView.set(NX - 1, j) = Cell::DIR;			
		}


		// referenceSolution
		std::vector<FP_TYPE> referenceSolution(NX * NY, 0);
		MatrixView<FP_TYPE> referenceSolutionView(referenceSolution, NX, NY);
		std::vector<int> offset = to_global_grid_coords(my_coords, {0,0});

		for (size_t j = 0; j != NY; ++j) 
		{
			for (size_t i = 0; i != NX; ++i) 
			{
				referenceSolutionView.set(i, j) = ParticularSolution( i * h, j * h);
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
						ParticularSolution((i + offset[0]) * h, (j + offset[1]) * h) * 4 * M_PI * M_PI;
			}
		}

		auto SolverJacobi = [](std::vector<FP_TYPE> &sol, std::vector<FP_TYPE> &sol2,
													std::vector<FP_TYPE> &rhs, const Stencil &stencil,
													size_t NX, size_t NY) {
			MatrixView<FP_TYPE> solView(sol, NX, NY);
			MatrixView<FP_TYPE> sol2View(sol2, NX, NY);
			MatrixView<FP_TYPE> rhsView(rhs, NX, NY);

			int topProc;
			int botProc;

			MPI_Cart_shift(g_topo_com, 1, 1, &topProc, &botProc);
			MPI_Status status;
			if(topProc >= 0)
			{	// send up
				MPI_Sendrecv(
					&solView.get(1, NY-1),	// send data start
					NX-2, 	// xmit length
					MPI_FP_TYPE,	//xmit dtypte
					topProc,	// rank of receiver
					0,	// Tag
					&sol2View.get(1, NY), // desstination buffers
					NX-2,	// Receive length
					MPI_FP_TYPE, // Receive dtype
					g_my_rank, // source rank
					0, // Tag
					g_topo_com,
					&status);
			}

			if(botProc >= 0)
			{	// send down
				MPI_Sendrecv(
					&solView.get(1, 1),	// send data start
					NX-2, 	// xmit length
					MPI_FP_TYPE,	//xmit dtypte
					botProc,	// rank of receiver
					0,	// Tag
					&sol2View.get(1, 0), // desstination buffers
					NX-2,	// Receive length
					MPI_FP_TYPE, // Receive dtype
					g_my_rank, // source rank
					0, // Tag
					g_topo_com,
					&status);
			}

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
			MPI_Barrier(g_topo_com);
			sol.swap(sol2);
			MPI_Barrier(g_topo_com);
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


		// local solution
		std::vector<FP_TYPE> solution2 = solution;
		std::cout << "solve LSE using stencil jacobi" << std::endl;
		auto start = std::chrono::high_resolution_clock::now();
		for (size_t iter = 0; iter <= iterations; ++iter) 
		{
			SolverJacobi(solution, solution2, rightHandSide, stencil, NX, NY);
		}

#ifdef USEMPI
		std::vector<FP_TYPE> recv_buf;
		if(g_my_rank == MASTER)
		{
			recv_buf.resize(g_resolution * g_resolution);			
		}		
		
		
		int x_start = borders[LEFT] 	== BORDER_GHOST ? 1 : 0;
		int x_end 	= borders[RIGHT] 	== BORDER_GHOST ? grid_size[COORD_X] - 1 : grid_size[COORD_X];
		int y_start = borders[BOTTOM] 	== BORDER_GHOST ? 1 : 0;
		int y_end 	= borders[TOP] 		== BORDER_GHOST ? grid_size[COORD_Y] - 1 : grid_size[COORD_Y];

		int len = 0;
		std::vector<FP_TYPE> send_buf(solution.size());
		for(int y = y_start; y < y_end; y++)
		{
			std::vector<FP_TYPE> next_line(&solution[x_start + grid_size[COORD_Y] * y], &solution[x_end + grid_size[COORD_Y] * y]);
			send_buf.insert(send_buf.begin() + len, next_line.begin(), next_line.end());
			len += next_line.size();
		}
		send_buf.resize(len);

		MPI_Gather(
			&send_buf,		// send_data
			send_buf.size(),// send_count
			MPI_FP_TYPE,	// send_datatype
			&recv_buf,		// recv_data
			recv_buf.size(),// recv_datatype
			MPI_FP_TYPE,	// send_datatype
			MASTER,			// root (rank of the receiver)
			g_topo_com		// communicator
		);
#endif

		auto stop = std::chrono::high_resolution_clock::now();
		auto seconds = std::chrono::duration_cast<std::chrono::duration<FP_TYPE>>(stop - start).count();
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
		std::cout << "--------------solver.hpp----------------\n";
		/*
		log.add("n_processes", std::to_string(g_n_processes));
		log.add("runtime", std::to_string(seconds));
		log.add("error", std::to_string(errorNorm));
		log.add("error", std::to_string(2308945720935784));
		*/
}
