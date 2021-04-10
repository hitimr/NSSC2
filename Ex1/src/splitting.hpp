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
    #define MPI_Comm void   // little hack so no errors appear when compiling without MPI
#endif

std::vector<int> border_types(int rank, int dim);
std::vector<int> get_prime_factors(int n);
size_t split_1D(int global_size, int splits, int pos, bool add_ghost_layers);


/* calculate the type of borders for a given rank

@param coords: cartesian coordinates within the mpi communicator. has only 1 entry in 1D equivalent to its rank

@return: vector of size 4 containing the type of each edge. Edges are assigned by
    [bottom, right, top, left]
    possible values are:
        BORDER_DIR: dirichlet border i.e. fixed
        BORDER_GHOST: ghost layer containing values of the neighbouring grid

*/
std::vector<int> border_types(const std::vector<int> & coords)
{
    // sanity check
    assert((coords.size() == 2) && "Invalid number of coordinates");

    std::vector<int> boundaries(4, BORDER_UNKNOWN);

    switch(g_dim)
    {
    case DIM1:  // 1D
        boundaries[LEFT] =   BORDER_DIR;
        boundaries[RIGHT] =  BORDER_DIR;
        boundaries[BOTTOM] = BORDER_GHOST;
        boundaries[TOP] =    BORDER_GHOST;

        // special case for bottom grid
        if(coords[COORD_Y] == 0)
        {
            boundaries[BOTTOM] = BORDER_DIR;
        }

        // special case for top grid
        if(coords[COORD_Y] == g_n_processes - 1)
        {
            boundaries[TOP] = BORDER_DIR;
        }
        break;

    case DIM2:  // 2D
        std::cerr <<  "2D is not implemented yet!" << std::endl;
        break;

    default:
        std::cerr <<  "Invalid dimension: " << coords.size() << std::endl; 
        break;
    }

    return boundaries;
}

// convenience overload fuction
std::vector<int> border_types(int n)
{
    assert(n >= 0);
    std::vector<int> coords = {0,n};
    return border_types(coords);
    
}


/* Calculate the local grid size of a given rank

If the resolution is not exactly divisible by the number of processes, grids of
different sizes are created such that the difference between each grid is minimal

local_grid_size() takes into account that every ghost layer requires an additional
line of grid points to store the data from the neighbouring grid.

@param coords: cartesian coordinates within the mpi communicator. has only 1 
entry in 1D equivalent to its rank

@return: a vector of size 2 containing the x, and y sizes of the local grid
    including ghost layers
*/
std::vector<size_t> local_grid_size(const std::vector<int> & coords, bool add_ghost_layers)
{
    // sanity check
    assert((coords.size() == 2) && "Invalid number of coordinates");

    // variables can't be created within a switch case so we need to do it here
    std::vector<size_t> size = {0, 0};  // local grid size
    std::vector<int> borders;
    std::vector<int> prime_factors;
    size_t x_dim = -1;
    size_t y_dim = -1;
    int remainder;
    int base_size;
    int n_x;
    int n_y;

    switch(g_dim)
    {
    case DIM1: // 1D  
        borders = border_types(coords);      
        base_size = (int) g_resolution / (int) g_n_processes;   // integer divinsion required
        remainder = g_resolution % g_n_processes; // number of grids with size + 1

        x_dim = g_resolution;
        y_dim = base_size;

        // bigger grids are allocated in ascending order
        // i.e. if 2 grids are bigger rank 0 and 1 have increased sizes
        if(coords[COORD_Y] < remainder)
        {
            y_dim++;
        } 

        // we need to add an additional line of of grid points for every ghost
        // layer in the subgrid to store the data from the neighbouring grid
        if(add_ghost_layers == true)
        {
            y_dim += std::count(borders.begin(), borders.end(), BORDER_GHOST);
        }        
        
        break;

    case DIM2:  // 2D
        prime_factors = get_prime_factors(g_n_processes);
        if(prime_factors.size() < 2)
        {
            // number of processes is a prime number. no splits possible
            // use 1D-split instead
            // HACK: local_grid_size doies not take g_dim as an arguiment so we 
            // retend to be in 1D for a moment to prevent an endless recursion
            g_dim = 1;  
            size = local_grid_size(coords, true);
            g_dim = 2;
            return size;
        }

        n_x = prime_factors[0]; // Number of splits in x-Direction
        n_y = g_n_processes / n_x; // Number of splits in y-Direction

        x_dim = split_1D(g_resolution, n_x, coords[COORD_X], add_ghost_layers);
        y_dim = split_1D(g_resolution, n_y, coords[COORD_Y], add_ghost_layers);

        if(add_ghost_layers == true)
        {
            
        }

        break;

    default:
        std::cerr <<  "Invalid dimension: " << g_dim << std::endl; 
    }   

    size = {2 * x_dim - 1, y_dim};

    return size;
}

// convenience overload function 
std::vector<size_t> local_grid_size(int n, bool add_ghost_layers)
{
    assert(n >= 0);
    std::vector<int> coords = {0, n};
    return local_grid_size(coords, add_ghost_layers);
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
size_t split_1D(int global_size, int splits, int pos, bool add_ghost_layers)
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


std::vector<int> to_global_grid_coords(const std::vector<int> & topo_coords, std::vector<int> local_grid_coords)
{
    int offset_y = 0;
    switch(g_dim)
    {
    case DIM1:
        // add heigths of the grids below
        for(int i = 0; i < topo_coords[1]; i++)
        {
            offset_y += local_grid_size(i, false)[COORD_Y];
        }
        local_grid_coords[COORD_Y] += offset_y;
        return local_grid_coords;

    default:
        assert(false && "Invalid number of dimensions");
    }
}


int get_neighbours(int direction)
{
#ifdef USEMPI
    int topProc, botProc, leftProc, rightProc;
    int returnProc = -1;

	MPI_Cart_shift(g_topo_com, 1, 1, &botProc, &topProc);
	MPI_Cart_shift(g_topo_com, 0, 1, &leftProc, &rightProc);
	if (direction == TOP)
	{
		returnProc=topProc;
	}
	else if (direction == BOTTOM)
	{
		returnProc=botProc;
	}
	else if (direction == LEFT)
    {
        returnProc=leftProc;
    }
	else if (direction == RIGHT)
    {
        returnProc=rightProc;
    }
	return returnProc;

#else
    return NO_NEIGHBOUR;    // serial has never neighbours
#endif
}


std::vector<int> get_topo_shape()
{
    std::vector<int> topo_shape(2);
    std::vector<int> prime_factors;

    switch(g_dim)
    {
    case DIM1:
        topo_shape[COORD_X] = 1;
        topo_shape[COORD_Y] = g_n_processes;
        break;

    case DIM2: 
        // calculate shape of 2D topology by using the prime factors of our number of procs
        prime_factors = get_prime_factors(g_n_processes);

        
        if(prime_factors.size() == 1)
        {
            // g_n_processes is prime -> revert to 1D
            topo_shape[COORD_X] = 1;
            topo_shape[COORD_Y] = g_n_processes;
        } 
        else
        {
            // g_n_processes is not prime. make grid shape such that y is minimal
            topo_shape[COORD_Y] = g_n_processes / prime_factors[COORD_X];
            topo_shape[COORD_X] = g_n_processes / topo_shape[COORD_Y];  
        }
        break;
        // TODO default case
    default:
        assert(false && "Invalid value for g_dim");
        break;
    }

    assert(topo_shape[COORD_X] * topo_shape[COORD_Y] == g_n_processes);

    return topo_shape;
}

