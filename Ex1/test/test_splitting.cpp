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

    return SUCCESS;
}


int test_borders_types()
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
    
    assert(borders_types(3) == top);
    assert(borders_types(2) == mid);
    assert(borders_types(1) == mid);
    assert(borders_types(0) == bot);

    return SUCCESS;
}


int main()
{

    assert(test_local_grid_size() == SUCCESS);
    assert(test_borders_types() == SUCCESS);

    cout << "All Tests passed!" << endl;

    return SUCCESS;    
}