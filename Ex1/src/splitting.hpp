#include <iostream>
#include <vector>

#define DIM1 1
#define DIM2 2

// globals according to command line
extern int g_n_processes;
extern int g_dim;
extern size_t g_resolution;


std::vector<size_t> grid_size(int rank)
{
    // 1D
    int base_size = (int) g_resolution / (int) g_n_processes;   // integer divinsion required
    int remainder = g_resolution % g_n_processes; // number of grids with size + 1

    size_t x_dim = g_resolution;
    size_t y_dim = base_size;

    // bigger grids are discributed in ascending order
    // i.e. if 2 grids are bigger rank 0 and 1 have increases sizes
    if(rank < remainder)
    {
        base_size++;
    } 
    
    std::vector<size_t> v = {x_dim, y_dim};
    return v;
}