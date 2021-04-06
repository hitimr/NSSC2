#include <assert.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

#ifdef USEMPI
#include <mpi.h>
#endif

// switch between float and double
#ifdef USE_FLOAT
  #define FP_TYPE float
  #define MPI_FP_TYPE MPI_FLOAT
#else
  #define FP_TYPE double
  #define MPI_FP_TYPE MPI_DOUBLE
#endif


#include "arguments.hpp"
#include "solver.hpp"
#include "common.hpp"


int g_my_rank = 0;
int g_n_processes = -1;
int g_dim = -1;
int g_iterations = -1;
size_t g_resolution = 0;

#ifdef USEMPI
    MPI_Comm g_topo_com;
#endif

 
int main(int argc, char *argv[]) 
{
#ifdef USEMPI
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &g_n_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_my_rank);
#else
  g_my_rank = 0;
  g_n_processes = 1;
#endif

  // redirect cout except for master to the output is less clogged
  streambuf *old = cout.rdbuf();   // save cout
  if(g_my_rank != MASTER)
  {
    stringstream ss;
    cout.rdbuf (ss.rdbuf());  // redirect cout to null
  }
  
  // parse command line arguments
  switch(argc)
  {
    case 3:
      g_dim = 1;
      g_resolution = convertTo<int>(1, 32, argc, argv);
      g_iterations = convertTo<int>(2, 1000, argc, argv);   
      break; 

    case 4:
      g_dim =        convertToDim(argv[1]);
      g_resolution = convertTo<int>(2, 32, argc, argv);
      g_iterations = convertTo<int>(3, 1000, argc, argv);
      break;
      
    default:
      cout << "Error! Invalid number of arguments" << endl;
      exit(-1);
  }
  
  if(g_my_rank == 0)
  {
    std::cout << "resolution=" << g_resolution << std::endl;
    std::cout << "iterations=" << g_iterations << std::endl;
  }

  // sanity checks
  assert(g_resolution > 0);
  assert(g_iterations > 0);
  assert((g_dim == 1) || (g_dim == 2));

  solve(g_resolution, g_iterations);

#ifdef USEMPI
  MPI_Finalize();
#endif

  // restore cout buffer
  if(g_my_rank != MASTER)
  {
    cout.rdbuf (old);
  }

  return 0;
}
