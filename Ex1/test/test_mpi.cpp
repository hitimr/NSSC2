#include <iostream>
#include <assert.h>
#include <mpi.h>
#include <vector>
#include "splitting.hpp"
#include "common.hpp"

#define SUCCESS 0
#define ARG_RESOLUTION 1

#ifndef USEMPI
    #define USEMPI
#endif

using namespace std;


int g_my_rank;
int g_n_processes;
int g_dim = 2;
size_t g_resolution;
MPI_Comm g_topo_com;


int get_neighbours(vector<int> coords, int direction)
{
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
}

int test_getNeighbours_2D()
{
    MPI_Barrier(MPI_COMM_WORLD); 
    if(g_my_rank == MASTER) cout << "Testing getNeighbours for 1D" << endl;

    if(g_my_rank == 5)
    {
    assert(get_neighbours(TOP) == NO_NEIGHBOUR);
    assert(get_neighbours(BOTTOM) == 4);
    }

    if(g_my_rank == 4)
    {   
    assert(get_neighbours(TOP) == 5);
    assert(get_neighbours(BOTTOM) == 3);
    }

    if(g_my_rank == 3)
    {   
    assert(get_neighbours(TOP) == 4);
    assert(get_neighbours(BOTTOM) == 2);
    }

    if(g_my_rank == 2)
    {
    assert(get_neighbours(TOP) == 3);
    assert(get_neighbours(BOTTOM) == 1);
    }

    if(g_my_rank == 1)
    {
    assert(get_neighbours(TOP) == 2);
    assert(get_neighbours(BOTTOM) == 0);
    }
    if(g_my_rank == 0)
    {
    assert(get_neighbours(TOP) == 1);
    assert(get_neighbours(BOTTOM) == NO_NEIGHBOUR);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(g_my_rank == MASTER) cout << "OK" << endl;
    return SUCCESS;
}




int test_getNeighbours_1D()
{
    MPI_Barrier(MPI_COMM_WORLD);
    if(g_my_rank == MASTER) cout << "Testing getNeighbours for 1D" << endl;

	if(g_my_rank == 5)
	{
    assert(get_neighbours(TOP) == NO_NEIGHBOUR);
    assert(get_neighbours(BOTTOM) == 4);
	}

	if(g_my_rank == 4)
    {   
    assert(get_neighbours(TOP) == 5);
    assert(get_neighbours(BOTTOM) == 3);
	}

    if(g_my_rank == 3)
    {   
    assert(get_neighbours(TOP) == 4);
    assert(get_neighbours(BOTTOM) == 2);
	}

	if(g_my_rank == 2)
    {
    assert(get_neighbours(TOP) == 3);
    assert(get_neighbours({0,2}, BOTTOM) == 1);
	}
	
	if(g_my_rank == 1)
    {
    assert(get_neighbours(TOP) == 2);
    assert(get_neighbours(BOTTOM) == 0);
	}

	if(g_my_rank == 0)
    {
    assert(get_neighbours(TOP) == 1);
    assert(get_neighbours(BOTTOM) == NO_NEIGHBOUR);
	}

    for(int i = 0; i < g_n_processes; i++)
    {
        assert(get_neighbours(LEFT) == NO_NEIGHBOUR);
        assert(get_neighbours(RIGHT) == NO_NEIGHBOUR);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(g_my_rank == MASTER) cout << "OK" << endl;
    
    return SUCCESS;
}

int test_1D()
{
    // Tests with 1D Grid
    // rearrange processes to vertically stacked topology
    if(g_my_rank == MASTER) cout << "Rearanging to vertical topology" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    vector<int> dims(2);
    vector<int> coords(2); // holds the rank's coordinates within the topology
    dims[0] = 1;
    dims[1] = g_n_processes;
    vector<int> periods = {false, false};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims.data(), periods.data(), true, &g_topo_com);	// https://www.mpich.org/static/docs/v3.3/www3/MPI_Cart_create.html
    MPI_Cart_coords(g_topo_com, g_my_rank, 2, coords.data());

    cout << "Rank " << g_my_rank << " now has coords (" << coords[0] << "|" << coords[1] << ")" << endl;


    /* Ranks should now be ordered:

    Rank 5 -> Top
    Rank 4
    Rank 3
    Rank 2
    Rank 1
    Rank 0 -> Bot

    */


    // Start tests
    MPI_Barrier(g_topo_com);

    assert(test_getNeighbours_1D() == SUCCESS);

    // TODO:
    //assert(test_bordertypes_1D() == SUCCESS);

    return SUCCESS;
}


int test_2D()
{
    // Tests with 2D Grid
    // Tests with 1D Grid
    // rearrange processes to vertically stacked topology
    if(g_my_rank == MASTER) cout << "Rearanging to vertical topology" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    vector<int> dims(2);
    vector<int> coords(2); // holds the rank's coordinates within the topology
    dims[0] = 2;
    dims[1] = g_n_processes/2;
    vector<int> periods = {false, false};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims.data(), periods.data(), true, &g_topo_com); // https://www.mpich.org/static/docs/v3.3/www3/MPI_Cart_create.html
    MPI_Cart_coords(g_topo_com, g_my_rank, 2, coords.data());

    cout << "Rank " << g_my_rank << " now has coords (" << coords[0] << "|" << coords[1] << ")" << endl;


    /* Ranks should now be ordered:

    Rank 5 -> Top
    Rank 4
    Rank 3
    Rank 2
    Rank 1
    Rank 0 -> Bot

    */


    // Start tests
    MPI_Barrier(g_topo_com);

    assert(test_getNeighbours_1D() == SUCCESS);

    // TODO:
    //assert(test_bordertypes_1D() == SUCCESS);

    // TODO: 2D tests
    return SUCCESS;
}

int main(int argc, char * argv[])
{

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &g_n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_my_rank);

    
    assert(g_n_processes == 6 && "Error: This test requires 6 processes!");

    assert(test_1D() == SUCCESS);
    //assert(test_2D() == SUCCESS);

    MPI_Finalize();

    cout << "All Tests passed on rank" << g_my_rank << "!" << endl;
    return SUCCESS;  
}
