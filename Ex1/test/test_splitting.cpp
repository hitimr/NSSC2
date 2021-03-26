#include <iostream>
#include <assert.h>
#include "splitting.hpp"

using namespace std;

int g_n_processes = -1;
int g_dim = -1;
size_t g_resolution = 0;


int test_local_grid_size()
{
    // check specific case
    g_n_processes = 3;
    g_dim = DIM1;
    g_resolution = 8;

    assert(local_grid_size(0)[COORD_X] == 8);
    assert(local_grid_size(0)[COORD_Y] == 3);

    assert(local_grid_size(1)[COORD_X] == 8);
    assert(local_grid_size(1)[COORD_Y] == 3);

    assert(local_grid_size(2)[COORD_X] == 8);
    assert(local_grid_size(2)[COORD_Y] == 2);


    // test if total is correct
    g_resolution = 500;
    g_n_processes = 31;
    size_t sum = 0;
    for(int rank = 0; rank < g_n_processes; rank++)
    {
        sum += local_grid_size(rank)[COORD_Y];
    }
    assert(sum == g_resolution);

    return 0;
}


int main()
{

    assert(test_local_grid_size() == 0);
    cout << "All Tests passed!" << endl;
    return 0;    
}