#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>

#include "common.hpp"

#ifdef USEMPI
	#include <mpi.h>
#else
    #define MPI_Comm void
#endif

std::vector<int> border_types(int rank);
std::vector<int> get_prime_factors(int n);
size_t split_1D(int global_size, int splits, int pos);


/* calculate the type of borders for a given rank

@param rank: rank of the grid/process who's border has to be calculated

@return: vector of size 4 containing the type of each edge. Edges are assigned by
    [bottom, right, top, left]
    possible values are:
        BORDER_DIR: dirichlet border i.e. fixed
        BORDER_GHOST: ghost layer containing values of the neighbouring grid

*/
std::vector<int> border_types(int rank)
{
    // sanity check
    assert(rank < g_n_processes);   // Error: invalid rank or g_n_processes is not properly set

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
        break;
    }

    return boundaries;
}


/* Calculate the local grid size of a given rank

If the resolution is not exactly divisible by the number of processes, grids of
different sizes are created such that the difference between each grid is minimal

local_grid_size() takes into account that every ghost layer requires an additional
line of grid points to store the data from the neighbouring grid.

@param rank: rank of the grid/process who's size has to be calculated

@return: a vector of size 2 containing the x, and y sizes of the local grid
    including ghost layers
*/
std::vector<size_t> local_grid_size(int rank, int dim=g_dim)
{
    // variables can't be created within a switch case so we need to do it here
    std::vector<size_t> size = {0, 0};
    std::vector<int> borders = border_types(rank);
    std::vector<int> prime_factors;
    size_t x_dim;
    size_t y_dim;
    int remainder;
    int base_size;
    int n_x;
    int n_y;
    int coords[2];	

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

        // we need to add an additional line of of grid points for every ghost
        // layer in the subgrid to store the data from the neighbouring grid
        y_dim += std::count(borders.begin(), borders.end(), BORDER_GHOST);
        
        break;

    case DIM2:  // 2D
        #ifdef USEMPI
            MPI_Cart_coords(g_topo_com, g_my_rank, DIM2, coords);
        #endif

        prime_factors = get_prime_factors(g_n_processes);
        if(prime_factors.size() < 2)
        {
            // number of processes is a prime number. no splits possible
            // use 1D-split instead
            size = local_grid_size(rank, DIM1);
        }

        n_x = prime_factors[0]; // Number of splits in x-Direction
        n_y = g_n_processes / n_x; // Number of splits in y-Direction


        x_dim = split_1D(g_resolution, n_x, coords[COORD_X]);
        y_dim = split_1D(g_resolution, n_x, coords[COORD_Y]);

        break;

    default:
        std::cerr <<  "Invalid dimension: " << g_dim << std::endl; 
    }   

    size = {x_dim, y_dim};

    return size;
}

/* Calculate the prime factors of a give integer n

# original taken from https://gist.githubusercontent.com/rohan-paul/3b0ef7d6ca9bbcfd3625132be1c29cdc/raw/a1937ce0d746ca4002522ec157561e8de4434701/prime-factors-of-number-simple-python.py

@param n: target number. must be > 0

@return: a vector containing all prime factors. factors appear in ascending order
    may contain the same prime multiple times. i.e.: 2*2*3 = 12
    if n is prime then factors = {n}
*/
std::vector<int> get_prime_factors(int n)
{
    assert(n > 0 && "Number mus be >0");

    std::vector<int> factors;

    while((n % 2) == 0)
    {
        factors.push_back(2);
        n = n / 2;
    }

    for(int i = 3; i < (int) sqrt(n) + 1; i += 2)
    {
        while((n % i) == 0)
        {
            factors.push_back(int(i));
            n = n / i;
        }
    }

    if(n > 2)
    {
        factors.push_back(int(n));
    }

    return factors;
}

/* Calculate the size of a given 1D split in an intervall


@param global_size: total length of the interval

@param splits: number of parts that the intervall is split into

@param pos: position of the sub-interval

*/
size_t split_1D(int global_size, int splits, int pos)
{
    assert(global_size > 0);
    assert(splits > 0);
    assert(pos >= 0);
    assert(pos < splits);

    int size = global_size / (int) splits; // integer division required
    int remainder = global_size % splits;

    // bigger grids are allocated in ascending order
    // i.e. if 2 grids are bigger rank 0 and 1 have increased sizes
    if(pos < remainder)
    {
        size++;
    }

    return (size_t) size;
}


