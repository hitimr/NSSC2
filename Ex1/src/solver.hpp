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

void print_matrixView(MatrixView<FP_TYPE> & mv)
{
	for(size_t y = 0; y < mv.M; y++)
	{
		for(size_t x = 0; x < mv.N; x++)
		{
			double val = mv.get(x,y);
			cout  << val << "\t";
		}
		cout << endl;
	}	
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
	Logger log(g_my_rank);

#ifdef USEMPI	
	// -------------------------------------------------------------------------
	// Prepare Communicator	
	// -------------------------------------------------------------------------

	// Create new communicator	
	// communicator is 2D by default. in 1D case x = 0 for all processes
	// this ensures that processes are stacked vertically in 1D
	vector<int> dims = get_topo_shape();	// shape of MPI topology
	if(dims[COORD_X] == 1)
	{
		// number of processes is prime (x = 1 for all processes)
		// -> revert to 1D
		g_dim = 1;
	}
	vector<int> periods = {false, false};	// no periodic grid
	MPI_Cart_create(MPI_COMM_WORLD, DIM2, dims.data(), periods.data(), true, &g_topo_com);	// https://www.mpich.org/static/docs/v3.3/www3/MPI_Cart_create.html
	MPI_Barrier(MPI_COMM_WORLD);	

	vector<int> coords(2); // coordinates of rank within the grid						 
	MPI_Cart_coords(g_topo_com, g_my_rank, DIM2, coords.data());

	
	auto grid_size = local_grid_size(coords, true);
	size_t NX = grid_size[COORD_X];
	size_t NY = grid_size[COORD_Y];
	FP_TYPE h = 1.0 / (g_resolution - 1);	
#else
	// Serial
	size_t NY = resolution;
	size_t NX = (2.0 * NY) - 1;
	FP_TYPE h = 1.0 / (NY - 1);
	vector<int> coords = {0,0};
#endif	
	
	// -------------------------------------------------------------------------
	// Prepare Stencil and Domain
	// -------------------------------------------------------------------------
	
	
	const auto stencil = Stencil(h);

	// domain cell types
	std::vector<int> domain(NX * NY, Cell::UNKNOWN);
	MatrixView<int> domainView(domain, NX, NY);

	// Fill borders depending on neigfbhours
	auto borders = border_types(coords);
	for (size_t i = 0; i != NX; ++i) 	
	{
		if(borders[TOP] == BORDER_DIR)	// top Border
			domainView.set(i, NY-1) = Cell::DIR;

		if(borders[BOTTOM] == BORDER_DIR)	// Bottom Border
			domainView.set(i, 0) = Cell::DIR;	
	}

	// left/right border
	for (size_t j = 0; j != NY; ++j)
	{		
		if(borders[LEFT] == BORDER_DIR)	// left border
			domainView.set(0, j) = Cell::DIR;

		if(borders[RIGHT] == BORDER_DIR) // righty border
			domainView.set(NX - 1, j) = Cell::DIR;
	}

	// quick sanity check
	if(g_my_rank == MASTER) 
	{	// Master is alway at bottom left
		assert(domainView.get(0,0) == Cell::DIR);
		assert(domainView.get(1,0) == Cell::DIR);	
		assert(domainView.get(0,1) == Cell::DIR);	
		assert(domainView.get(1,1) == Cell::UNKNOWN);		
	}
	
	// right hand side
	std::vector<FP_TYPE> rightHandSide(NX * NY, 0);
	MatrixView<FP_TYPE> rightHandSideView(rightHandSide, NX, NY);
	auto coord_offset = to_global_grid_coords(coords, {(int) 0, (int) 0});
	for (size_t j = 0; j != NY; ++j) 
	{
		for (size_t i = 0; i != NX; ++i) 
		{
			rightHandSideView.set(i, j) =
					ParticularSolution((coord_offset[0] + i) * h, (coord_offset[1] + j) * h) * 4 * M_PI * M_PI;
		}
	}

	// -------------------------------------------------------------------------
	// Jacobi Iteration Cycle 
	// -------------------------------------------------------------------------

	
	auto SolverJacobi = [](std::vector<FP_TYPE> &sol, std::vector<FP_TYPE> &sol2,
												std::vector<FP_TYPE> &rhs, const Stencil &stencil,
												size_t NX, size_t NY, vector<int> & neighbours) {
		MatrixView<FP_TYPE> solView(sol, NX, NY);
		MatrixView<FP_TYPE> sol2View(sol2, NX, NY);
		MatrixView<FP_TYPE> rhsView(rhs, NX, NY);

		// Iterate over entire grid
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
#ifdef USEMPI
	MPI_Request req;
	MPI_Status status;

	// TODO 2D
	// sync borders
	if(neighbours[BOTTOM] != NO_NEIGHBOUR) // send down
		MPI_Isend(&solView.get(1, 1),	NX-2, MPI_FP_TYPE, neighbours[BOTTOM], 0, g_topo_com, &req);	

	if(neighbours[TOP] != NO_NEIGHBOUR) // send up
		MPI_Isend(&solView.get(1, NY-2), NX-2, MPI_FP_TYPE, neighbours[TOP], 0, g_topo_com, &req);	
	
	if(neighbours[BOTTOM] != NO_NEIGHBOUR) // receivce from bot
		MPI_Recv(&solView.get(1, 0), NX-2, MPI_FP_TYPE, neighbours[BOTTOM], 0, g_topo_com, &status);	

	if(neighbours[TOP] 	!= NO_NEIGHBOUR) // receivce from top
		MPI_Recv(&solView.get(1, NY-1), NX-2, MPI_FP_TYPE, neighbours[TOP] , 0, g_topo_com, &status);	
#endif
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

	// -------------------------------------------------------------------------
	// Jacobi Iteration Cycle 
	// -------------------------------------------------------------------------


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


	// -------------------------------------------------------------------------
	// Initial Fill of the grid
	// -------------------------------------------------------------------------


	// solution approximation starting with boundary initialized to dirichlet
	// conditions, else 0
	std::vector<FP_TYPE> solution(NX * NY, 0);
	MatrixView<FP_TYPE> solutionView(solution, NX, NY);
	coord_offset = to_global_grid_coords(coords, {(int) 0, (int) 0});
	for (size_t j = 0; j != NY; ++j) {
		for (size_t i = 0; i != NX; ++i) {
			if (domainView.get(i, j) == Cell::DIR)
				solutionView.set(i, j) = ParticularSolution((coord_offset[0] + i) * h, (coord_offset[1] + j) * h);
		}
	}

	// check who is my neigbhour
	vector<int> neighbours(4);
	for(size_t i = 0; i < neighbours.size(); i++)
		neighbours[i] = get_neighbours(i);


	// ------------------------------------------------------------------------- 
	// Iteration Loop
	// ------------------------------------------------------------------------- 



	std::vector<FP_TYPE> solution2 = solution;
	std::vector<double> iteration_times(g_iterations, -1);	// containts the time for every iteration

	if(g_my_rank == MASTER) std::cout << "solve LSE using stencil jacobi" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();		
	for (size_t iter = 0; iter <= iterations; ++iter) 
	{
		auto iteration_start = std::chrono::high_resolution_clock::now();	

		// computer iteraion
		SolverJacobi(solution, solution2, rightHandSide, stencil, NX, NY, neighbours);

#ifdef USEMPI			
		MPI_Barrier(g_topo_com);
#endif
		auto iteration_end = std::chrono::high_resolution_clock::now();			
		double seconds = std::chrono::duration_cast<std::chrono::duration<FP_TYPE>>(iteration_end - iteration_start).count();
		iteration_times[iter] = seconds;
	}
		
	if(PRINT_LOCAL_SOLUTIONS)
	{
		cout << "writing local solution of rank " << g_my_rank << endl;
		std::string fileName = "out/solution" + std::to_string(g_my_rank) + ".txt";
		print_matrix(solution, "fileName", NX, NY);	
	}			

	// -------------------------------------------------------------------------
	// Collect results on root
	// -------------------------------------------------------------------------

#ifdef USEMPI	
	// check if we need to remove the ghost layer before sending it
	int x_start = borders[LEFT] 	== BORDER_GHOST ? 1 : 0;
	int x_end 	= borders[RIGHT] 	== BORDER_GHOST ? grid_size[COORD_X] - 1 : grid_size[COORD_X];
	int y_start = borders[BOTTOM] 	== BORDER_GHOST ? 1 : 0;
	int y_end 	= borders[TOP] 		== BORDER_GHOST ? grid_size[COORD_Y] - 1 : grid_size[COORD_Y];

	// fill send buffer without ghost layers
	std::vector<FP_TYPE> send_buf;
	for(int y = y_start; y < y_end; y++)
	{
		// copy each line into the send-buffer
		std::vector<FP_TYPE> next_line(&solution[x_start + NX * y], &solution[x_end + NX* y]); // exract line without ghost layers
		send_buf.insert(send_buf.end(), next_line.begin(), next_line.end());	// add to send buffer
	}
	
	// Gather length of all send buffers
	int send_buf_size = send_buf.size();
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

	// Create disaplcement vector on root
	int offset = 0;
	vector<int> displs(g_n_processes);
	for(size_t i = 0; i < recvcounts.size(); i++)	// recvcounts.size() > 0 only on root
	{
		displs[i] = offset;
		offset += recvcounts[i] + 0;
	}

	// Gather all solutions. Gatherv required because grid sizes differ between ranks
	vector<FP_TYPE> rbuf(g_my_rank == MASTER ? offset : 0);	// only allocate memory on master		
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
#ifdef USEMPI
		// ---------------------------------------------------------------------
		//Assemble final solution on root/master
		// ---------------------------------------------------------------------


		// change to global coordinates
		NY = g_resolution;
		NX = (2.0 * NY) - 1;
		h = 1.0 / (NY - 1);

		// Assemble solution
		solution.clear();
		solution.insert(solution.begin(), rbuf.begin(), rbuf.begin() + offset);		
#endif

		MatrixView<FP_TYPE> solution_view(solution, NX, NY);

		// coimputer global right hand side
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

					
		// global referenceSolution
		std::vector<FP_TYPE> global_referenceSolution(NX * NY, 0);
		MatrixView<FP_TYPE> global_referenceSolutionView(global_referenceSolution, NX, NY);

		for (size_t j = 0; j != NY; ++j) 
		{
			for (size_t i = 0; i != NX; ++i) 
			{
				global_referenceSolutionView.set(i, j) = ParticularSolution( i * h, j * h);
			}
		}

		auto stop = std::chrono::high_resolution_clock::now();
		cout << "Calculations finished. collecting data for benchmark" << endl;

		// ---------------------------------------------------------------------
		// Computer Results
		// ---------------------------------------------------------------------

		// errors
		auto residual = ComputeResidual(solution, global_rightHandSide, stencil, NX, NY);
		auto error = ComputeError(solution, global_referenceSolution, NX, NY);
		auto residualNorm = NormL2(residual);
		auto residualMax = NormInf(residual);
		auto errorNorm = NormL2(error);
		auto errorMax = NormInf(error);
		
		// runtimes
		
		auto seconds = std::chrono::duration_cast<std::chrono::duration<FP_TYPE>>(stop - start).count();
		double average_iteration_time = 0;
		for(size_t i = 0; i < iteration_times.size(); i++)
		{
			average_iteration_time += iteration_times[i];
		}
		average_iteration_time = average_iteration_time / iteration_times.size();

		
		// ---------------------------------------------------------------------
		// Write Results
		// ---------------------------------------------------------------------

		// console/log output
#ifdef USEMPI
		std::cout << "used MPI: true" <<  endl;
		log.add("mpi", "true");
#else
		std::cout << "used MPI: false" <<  endl;
		log.add("mpi", "false");
#endif


#ifdef USE_FLOAT
		std::cout << "dtype = FLOAT" << endl;
		log.add("dtype", "float");
#else
		std::cout << "dtype = DOUBLE" << endl;
		log.add("dtype", "double");
#endif
		std::cout << std::scientific << "|residual|=" << residualNorm << std::endl;		
		std::cout << std::scientific << "|residualMax|=" << residualMax << std::endl;			
		std::cout << std::scientific << "|error|2=" << errorNorm << std::endl;		
		std::cout << std::scientific << "|errorMax|inf=" << errorMax << std::endl;
		std::cout << std::scientific << "runtime=" << seconds << std::endl;
		std::cout << std::scientific << "average time/iteration=" << average_iteration_time << std::endl << endl;

		log.add("n_processes", std::to_string(g_n_processes));
		log.add("total runtime", std::to_string(seconds));
		log.add("error", std::to_string(errorNorm));
		log.add("errorMax", std::to_string(errorMax));
		log.add("residual", std::to_string(residualNorm));
		log.add("residualMax", std::to_string(residualMax));
		log.add("average_iteration_time", std::to_string(average_iteration_time));

#ifdef USE_FLOAT
		
#else
		
#endif
		
	}

#ifdef USEMPI
	MPI_Barrier(g_topo_com);
#endif

}




