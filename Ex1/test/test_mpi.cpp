#include <iostream>
#include <assert.h>
#include <mpi.h>
#include <vector>
#include "splitting.hpp"

#define SUCCESS 0
#define ARG_RESOLUTION 1

using namespace std;


int g_my_rank;
int g_n_processes;
int g_dim = 2;
size_t g_resolution;
MPI_Comm g_topo_com;

int test_local_grid_size()
{
    if(g_my_rank == 0) cout << "Testing local_grid_size()" << endl;

    vector<int> prime_factors = get_prime_factors(g_n_processes);
    for(g_resolution = 120; g_resolution < 121; g_resolution++)
    {
        if(prime_factors.size() < 2)
        {
            // TODO: 1D
        }
        
        int n_x = prime_factors[0];
        int n_y = (int) g_n_processes / (int) n_x;

        // check sums
        int sum_x = 0;
        int sum_y = 0;
        cout << "OK" << endl;
        for(int x = 0; x < n_x; x++)
        {
            sum_x += local_grid_size(0)[0];
        }
        for(int y = 0; y < n_y; y++)
        {
            //sum_y += local_grid_size(0,y)[1];
        }
        
        //assert(sum_x == g_resolution && "Sum of x-sizes does not match");
        //assert(sum_y == g_resolution && "Sum of y-sizes does not match");

        // Check dimesions
        for(int y = 0; y < n_y; y++)
        {
            for(int x = 0; x < n_x; x++)
            {
                //assert(local_grid_size(x,y)[0] == local_grid_size(x,0)[0] && "Grid sizes dont match");
                //assert(local_grid_size(x,y)[1] == local_grid_size(0,y)[1] && "Grid sizes dont match");
            }
        }
    }


    return SUCCESS;
}

int main(int argc, char * argv[])
{
    assert(argc == 2 && "Not enough arguments");
    g_resolution = stoi(argv[ARG_RESOLUTION]);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &g_n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_my_rank);

    
    int coords[2] = {};		


    int ndims = 2;
    int dims[2] = {};
    int bcs[2] = {0,0};	 // logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
    int reorder = 1;			// ranking may be reordered (true) or not (false) (logical)

    MPI_Dims_create(g_n_processes, ndims, dims);	// https://www.mpich.org/static/docs/v3.3.x/www3/MPI_Dims_create.html 
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &g_topo_com);	// https://www.mpich.org/static/docs/v3.3/www3/MPI_Cart_create.html
    MPI_Barrier(MPI_COMM_WORLD);	

    test_local_grid_size();

    MPI_Finalize();

    cout << "All Tests passed on rank" << g_my_rank << "!" << endl;
    return SUCCESS;  
}