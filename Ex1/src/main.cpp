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
#else
  #define FP_TYPE double
#endif


#include "arguments.hpp"
#include "solver.hpp"

int g_my_rank = 0;
int g_n_processes = -1;
int g_dim = -1;
int g_iterations = -1;
size_t g_resolution = 0;


 
int main(int argc, char *argv[]) 
{
#ifdef USEMPI
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &g_n_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_my_rank);
#endif
  
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

  solve(g_resolution, g_iterations, g_my_rank);

#ifdef USEMPI
  MPI_Finalize();
#endif

  return 0;
}
