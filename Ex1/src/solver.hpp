#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <assert.h>
#include <fstream>
#include <string>


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

void print_matrixView(MatrixView<FP_TYPE> & mv, std::string fileName)
{
	ofstream out(fileName);
	for(size_t y = 0; y < mv.M; y++)
	{
		for(size_t x = 0; x < mv.N; x++)
		{
			double val = mv.get(x,y);
			out  << val << "\t";
		}
		out << endl;
	}
	out.close();
}


void print_matrix(std::vector<FP_TYPE> &v, std::string fileName, size_t NX, size_t NY)
{
	MatrixView<FP_TYPE> mv(v, NX, NY);
	print_matrixView(mv, fileName);
}

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
		
		
		// Create new communicator
		
		int ndims = g_dim;
		std::vector<int> dims;
		dims.resize(g_dim); 
		int reorder = 1;			// ranking may be reordered (true) or not (false) (logical)
		std::vector<int> bcs;	// logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
		bcs.resize(g_dim, 0);
		
		MPI_Dims_create(g_n_processes, ndims, dims.data());	// https://www.mpich.org/static/docs/v3.3.x/www3/MPI_Dims_create.html 
		MPI_Cart_create(MPI_COMM_WORLD, ndims, dims.data(), bcs.data(), reorder, &g_topo_com);	// https://www.mpich.org/static/docs/v3.3/www3/MPI_Cart_create.html
		MPI_Barrier(MPI_COMM_WORLD);	

		int coords[2] = {};								 
		MPI_Cart_coords(g_topo_com, g_my_rank, ndims, coords);		// TODO provide vector immediately
		if(g_my_rank == MASTER)
		{
			std::cout << "dims= (";
			for(int i = 0; i < g_dim; i++) std::cout << dims[i] << ",";
			std::cout << ")" << endl;
		}

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

		std::cout << "Rank=" << g_my_rank << "; coords =" << coords[0] <<std::endl;
		
		auto grid_size = local_grid_size(my_coords, true);
		size_t NX = grid_size[COORD_X];
		size_t NY = grid_size[COORD_Y];
		FP_TYPE h = 1.0 / (g_resolution - 1);	
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
		
		// right hand side
		std::vector<FP_TYPE> rightHandSide(NX * NY, 0);
		MatrixView<FP_TYPE> rightHandSideView(rightHandSide, NX, NY);
		auto coord_offset = to_global_grid_coords(my_coords, {(int) 0, (int) 0});
		for (size_t j = 0; j != NY; ++j) 
		{
			for (size_t i = 0; i != NX; ++i) 
			{
				
//if(g_my_rank == DEBUG_RANK)
	//cout << "(" << global_coords[0] << "|" << global_coords[1] << ")";

				rightHandSideView.set(i, j) =
						ParticularSolution(i * h, j * h) * 4 * M_PI * M_PI;	// TODO calculate fixed offset beforehand
			}
//if(g_my_rank == DEBUG_RANK) cout << endl;
		}


		auto SolverJacobi = [](std::vector<FP_TYPE> &sol, std::vector<FP_TYPE> &sol2,
													std::vector<FP_TYPE> &rhs, const Stencil &stencil,
													size_t NX, size_t NY) {
			MatrixView<FP_TYPE> solView(sol, NX, NY);
			MatrixView<FP_TYPE> sol2View(sol2, NX, NY);
			MatrixView<FP_TYPE> rhsView(rhs, NX, NY);

			int topProc;
			int botProc;

			//MPI_Cart_shift(g_topo_com, 0, 1, &g_my_rank, &botProc);

			MPI_Request req;
			MPI_Status status;

			botProc = 0;
			topProc = 1;

			if(g_dim == DIM1)
			{				
				if(g_my_rank == 1) MPI_Isend(&solView.get(1, 1), 		NX-2, MPI_FP_TYPE, botProc, 0, g_topo_com, &req);	// send down
				if(g_my_rank == 0) MPI_Isend(&solView.get(1, NY-2), 	NX-2, MPI_FP_TYPE, topProc, 0, g_topo_com, &req);	// send up
				
				if(g_my_rank == 1) MPI_Recv(&solView.get(1, 0), 		NX-2, MPI_FP_TYPE, botProc, 0, g_topo_com, &status);	// receivce from bot
				if(g_my_rank == 0) MPI_Recv(&solView.get(1, NY-1), 		NX-2, MPI_FP_TYPE, topProc, 0, g_topo_com, &status);	// receivce from top
			}
			else
			{
				// TODO 2D
			}	

	
			MPI_Barrier(g_topo_com);
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
if(g_my_rank == DEBUG_RANK)
print_matrixView(solView, "out/solution.txt");	
	
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
#ifdef USEMPI			
			MPI_Barrier(g_topo_com);
#endif
		}
string fileName = "out/solution"+std::to_string(g_my_rank)+".txt";
print_matrix(solution, fileName, NX, NY);

#ifdef USEMPI		
		
		// check if we need to remove the ghost layer before sending it
		int x_start = borders[LEFT] 	== BORDER_GHOST ? 1 : 0;
		int x_end 	= borders[RIGHT] 	== BORDER_GHOST ? grid_size[COORD_X] - 1 : grid_size[COORD_X];
		int y_start = borders[BOTTOM] 	== BORDER_GHOST ? 1 : 0;
		int y_end 	= borders[TOP] 		== BORDER_GHOST ? grid_size[COORD_Y] - 1 : grid_size[COORD_Y];

		std::vector<FP_TYPE> send_buf;
		for(int y = y_start; y < y_end; y++)
		{
			// copy each line into the send-buffer
			// since this buffer could be smaller than the solution due to ghost layers we need to keep track of its length
			std::vector<FP_TYPE> next_line(&solution[x_start + grid_size[COORD_Y] * y], &solution[x_end + grid_size[COORD_Y] * y]); // exract line without ghost layers
			send_buf.insert(send_buf.end(), next_line.begin(), next_line.end());	// add to send buffer
		}
		
		
		int send_buf_size = send_buf.size();
		cout << "rank " << g_my_rank << " send_buf_size=" << send_buf_size << endl;

		vector<int> recvcounts(g_my_rank == MASTER ? g_n_processes : 0);	// only allocate memory on master
		MPI_Gather(
			&send_buf_size,		// send_data
			1,			// send_count
			MPI_INT,	// send_datatype
			g_my_rank == MASTER ? recvcounts.data() : nullptr,		// recv_data
			g_my_rank == MASTER ? 1 : 0,		// recv_count
			MPI_INT,		// send_datatype
			MASTER,			// root (rank of the receiver)
			g_topo_com		// communicator
		);				

		int offset = 0;
		vector<int> displs(g_n_processes);
		for(size_t i = 0; i < recvcounts.size(); i++)
		{
			displs[i] = offset;
			offset += recvcounts[i] + 0;
		}

		vector<FP_TYPE> rbuf(g_my_rank == MASTER ? offset : 0);	// only allocate memory on master
cout << "Rank " << g_my_rank << "sending"  << endl;
		
		MPI_Gatherv(
			send_buf.data(),// sendbuf
			send_buf.size(),// sendcount
			MPI_FP_TYPE,	// MPI_Datatype
			g_my_rank == MASTER ? rbuf.data() : nullptr,	// recvbuf
			recvcounts.data(), 	// recvcounts
			displs.data(),		// displs
			MPI_FP_TYPE,
			MASTER,				// root (rank of the receiver)
			g_topo_com			// communicator
		);




#endif	
		if(g_my_rank == MASTER)
		{	
			// change to global coordinates
			NY = g_resolution;
			NX = (2.0 * NY) - 1;
			h = 1.0 / (NY - 1);


		
			// Assemble solution
			solution.clear();
			solution.insert(solution.begin(), rbuf.begin(), rbuf.begin() + offset);		

print_matrix(solution, "out/solution.txt", NX, NY);


			MatrixView<FP_TYPE> solution_view(solution, NX, NY);

			// right hand side
			std::vector<FP_TYPE> global_rightHandSide(NX * NY, 0);
			MatrixView<FP_TYPE>  global_rightHandSideView(global_rightHandSide, NX, NY);
			for (size_t j = 0; j != NY; ++j) 
			{
				for (size_t i = 0; i != NX; ++i) 
				{
					global_rightHandSideView.set(i, j) =
							ParticularSolution(i * h, j * h) * 4 * M_PI * M_PI;
				}
			}

						
			// referenceSolution
			std::vector<FP_TYPE> global_referenceSolution(NX * NY, 0);
			MatrixView<FP_TYPE> global_referenceSolutionView(global_referenceSolution, NX, NY);

			for (size_t j = 0; j != NY; ++j) 
			{
				for (size_t i = 0; i != NX; ++i) 
				{
					global_referenceSolutionView.set(i, j) = ParticularSolution( i * h, j * h);
				}
			}

//print_matrix(global_referenceSolution, "out/solution.txt", NX, NY);

			auto stop = std::chrono::high_resolution_clock::now();
			auto seconds = std::chrono::duration_cast<std::chrono::duration<FP_TYPE>>(stop - start).count();
			std::cout << std::scientific << "runtime=" << seconds << std::endl;

			
			auto residual = ComputeResidual(solution, global_rightHandSide, stencil, NX, NY);
			auto residualNorm = NormL2(residual);
			std::cout << std::scientific << "|residual|=" << residualNorm << std::endl;
			auto residualMax = NormInf(residual);
			std::cout << std::scientific << "|residualMax|=" << residualMax
								<< std::endl;
			auto error = ComputeError(solution, global_referenceSolution, NX, NY);
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
		else
		{
			//MPI_Send(&send_buf[0], 95, MPI_FP_TYPE, 0, 0, g_topo_com);
		}
		MPI_Barrier(g_topo_com);
}




