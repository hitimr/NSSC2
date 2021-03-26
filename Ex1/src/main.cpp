#include <assert.h>
#include <iomanip>
#include <iostream>
#include <vector>
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


 
int main(int argc, char *argv[]) {

  // MPI stuff
  int rank = 0;
  int numproc = 1;

  // Grid variables
  int resolution = -1;
  int iterations = -1;
  int dim = 1;  // number of dimensions the grid will be decomnposed in (1D or 2D)


#ifdef USEMPI
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  cout << "argc=" << argc << endl;
  // parse command line arguments
  switch(argc)
  {
    case 3:
      resolution =  convertTo<int>(1, 32, argc, argv);
      iterations =  convertTo<int>(2, 1000, argc, argv);   
      break; 

    case 4:
      dim =        convertToDim(argv[1]);
      resolution =  convertTo<int>(2, 32, argc, argv);
      iterations =  convertTo<int>(3, 1000, argc, argv);
      break;
      
    default:
      cout << "Error! Invalid number of arguments" << endl;
      exit(-1);
  }


  
  if(rank == 0)
  {
    std::cout << "numproc=" << numproc << std::endl;
    std::cout << "resolution=" << resolution << std::endl;
    std::cout << "iterations=" << iterations << std::endl;
  }

  // sanity checks
  assert(resolution > 0);
  assert(iterations > 0);
  assert((dim == 1) || (dim == 2));

  solve(resolution, iterations, rank, numproc);

#ifdef USEMPI
  MPI_Finalize();
#endif

  return 0;
}
