#include <iostream>
#include <vector>

#define DIM1 1
#define DIM2 2
#define COORD_X 0
#define COORD_Y 1

// globals according to command line
extern int g_n_processes;
extern int g_dim;
extern size_t g_resolution;


std::vector<size_t> grid_size(int rank)
{
    std::vector<size_t> size = {0, 0};
    size_t x_dim;
    size_t y_dim;
    int remainder;
    int base_size;

    switch(g_dim)
    {
    case DIM1:
        // 1D
        base_size = (int) g_resolution / (int) g_n_processes;   // integer divinsion required
        remainder = g_resolution % g_n_processes; // number of grids with size + 1

        x_dim = g_resolution;
        y_dim = base_size;

        // bigger grids are discributed in ascending order
        // i.e. if 2 grids are bigger rank 0 and 1 have increases sizes
        if(rank < remainder)
        {
            y_dim++;
        } 
        
        size = {x_dim, y_dim};
        break;

    case DIM2:
        std::cerr <<  "2D is not implemented yet!" << std::endl;
        break;

    default:
        std::cerr <<  "Invalid dimension: " << g_dim << std::endl; 
        break;
    }   
    return size;
}


std::vector<size_t> ghost_layers(int rank)
{
    // 1D
    std::vector<size_t> boundaries = {0, 0};

    return boundaries;
}
