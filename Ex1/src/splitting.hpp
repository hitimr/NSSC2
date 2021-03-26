#include <iostream>
#include <vector>

// Dimensions
#define DIM1 1
#define DIM2 2
#define COORD_X 0
#define COORD_Y 1

// Border Types
// Assignemnt resembles the enumerator from solver.hpp,
// enum Cell { UNKNOWN = 0, DIR = 1, NEU = 2, ROB = 0 };
#define BORDER_UNKNOWN 0
#define BORDER_DIR     1    // Dirichlet
#define BORDER_NEU     2    // Neumann
#define BORDER_ROB     0    // Robin
#define BORDER_GHOST   3    // Ghost Layer

// Directional Assignemtns
#define BOTTOM  0
#define RIGHT   1
#define TOP     2
#define LEFT    3

// globals according to command line
extern int g_n_processes;
extern int g_dim;
extern size_t g_resolution;


std::vector<size_t> local_grid_size(int rank)
{
    std::vector<size_t> size = {0, 0};
    size_t x_dim;
    size_t y_dim;
    int remainder;
    int base_size;

    switch(g_dim)
    {
    case DIM1: // 1D        
        base_size = (int) g_resolution / (int) g_n_processes;   // integer divinsion required
        remainder = g_resolution % g_n_processes; // number of grids with size + 1

        x_dim = g_resolution;
        y_dim = base_size;

        // bigger grids are allocated in ascending order
        // i.e. if 2 grids are bigger rank 0 and 1 have increased sizes
        if(rank < remainder)
        {
            y_dim++;
        } 
        
        size = {x_dim, y_dim};
        break;

    case DIM2:  // 2D
        std::cerr <<  "2D is not implemented yet!" << std::endl;
        break;

    default:
        std::cerr <<  "Invalid dimension: " << g_dim << std::endl; 
        break;
    }   

    return size;
}


std::vector<int> borders_types(int rank)
{
    std::vector<int> boundaries(4, BORDER_UNKNOWN);

    switch(g_dim)
    {
    case DIM1:  // 1D
        boundaries[LEFT] =   BORDER_DIR;
        boundaries[RIGHT] =  BORDER_DIR;
        boundaries[BOTTOM] = BORDER_GHOST;
        boundaries[TOP] =    BORDER_GHOST;

        // special case for bottom grid
        if(rank == 0)
        {
            boundaries[BOTTOM] = BORDER_DIR;
        }

        // special case for top grid
        if(rank == g_n_processes - 1)
        {
            boundaries[TOP] = BORDER_DIR;
        }
        break;

    case DIM2:  // 2D
        std::cerr <<  "2D is not implemented yet!" << std::endl;
        break;

    default:
        std::cerr <<  "Invalid dimension: " << g_dim << std::endl; 
    }

    return boundaries;
}
