#include <iostream>
#include <assert.h>
#include <mpi.h>
#include "splitting.hpp"

#define SUCCESS 0

using namespace std;


int g_my_rank;
int g_n_processes;
int g_dim;
size_t g_resolution;
MPI_Comm g_topo_com;

int test_local_grid_size()
{
    if(g_my_rank == 0) cout << "Testing local_grid_size()" << endl;



    return SUCCESS;
}

int main()
{
    g_my_rank = 0;
    g_n_processes = -1;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &g_n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_my_rank);

    test_local_grid_size();
    
    MPI_Finalize();
}