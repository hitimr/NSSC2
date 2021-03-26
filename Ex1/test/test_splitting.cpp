#include <iostream>
#include <assert.h>
#include <algorithm>
#include "splitting.hpp"

#define SUCCESS 0

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

    /* Assignemnt should be:
        rank 2: 2 + 1 Ghost Layer   (top)
        rank 1: 3 + 2 Ghost layers  (mid)
        rank 0: 3 + 1 Ghost Layer   (bot)
    */
    assert(local_grid_size(2)[COORD_X] == 8);
    assert(local_grid_size(2)[COORD_Y] == 3);

    assert(local_grid_size(1)[COORD_X] == 8);
    assert(local_grid_size(1)[COORD_Y] == 5);

    assert(local_grid_size(0)[COORD_X] == 8);
    assert(local_grid_size(0)[COORD_Y] == 4);


    // test if total is correct
    g_resolution = 500;
    g_n_processes = 31;
    size_t sum = 0;
    for(int rank = 0; rank < g_n_processes; rank++)
    {
        sum += local_grid_size(rank)[COORD_Y];
    }

    // substract ghost layers and compare to original resolution
    sum -= (g_n_processes * 2);   // substract two ghost layers for evey grid
    sum += 2;    // add 1 ghost layer for top and bottom grid
    assert(sum == g_resolution);

    return SUCCESS;
}


int test_border_types()
{
    // 1D
    // check specific case
    g_n_processes = 4;
    g_dim = DIM1;
    g_resolution = 100; 

    // assignemnt is [bottom, right, top, left]
    vector<int> top = {BORDER_GHOST, BORDER_DIR, BORDER_DIR,   BORDER_DIR};    
    vector<int> mid = {BORDER_GHOST, BORDER_DIR, BORDER_GHOST, BORDER_DIR};
    vector<int> bot = {BORDER_DIR,   BORDER_DIR, BORDER_GHOST, BORDER_DIR};
    
    assert(border_types(3) == top);
    assert(border_types(2) == mid);
    assert(border_types(1) == mid);
    assert(border_types(0) == bot);

    return SUCCESS;
}


int main()
{

    assert(test_local_grid_size() == SUCCESS);
    assert(test_border_types() == SUCCESS);

    cout << "All Tests passed!" << endl;

    return SUCCESS;    
}