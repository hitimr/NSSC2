#include <iostream>
#include <assert.h>
#include "splitting.hpp"

#define SUCCESS 0

using namespace std;

// globals
int g_n_processes = -1;
int g_dim = -1;
size_t g_resolution = 0;


int test_local_grid_size()
{
    cout << "Testing local_grid_size()" << endl;

    // 1D

    // check specific case
    g_n_processes = 3;
    g_dim = DIM1;
    g_resolution = 8;

    /* Assignemnt should be:
        rank 2: 2 + 1 Ghost Layer   (top)
        rank 1: 3 + 2 Ghost layers  (mid)
        rank 0: 3 + 1 Ghost Layer   (bot)
    */
    assert(local_grid_size(2, true)[COORD_X] == 15);
    assert(local_grid_size(2, true)[COORD_Y] == 3);

    assert(local_grid_size(1, true)[COORD_X] == 15);
    assert(local_grid_size(1, true)[COORD_Y] == 5);

    assert(local_grid_size(0, true)[COORD_X] == 15);
    assert(local_grid_size(0, true)[COORD_Y] == 4);


    // check serial case
    g_n_processes = 1;
    g_dim = DIM1;
    g_resolution = 8;
    assert(local_grid_size(0, true)[COORD_X] == 2*g_resolution - 1);
    assert(local_grid_size(0, true)[COORD_Y] == g_resolution);


    // test if total is correct
    g_resolution = 500;
    g_n_processes = 31;
    size_t sum = 0;
    for(int rank = 0; rank < g_n_processes; rank++)
    {
        sum += local_grid_size(rank, true)[COORD_Y];
    }

    // substract ghost layers and compare to original resolution
    sum -= (g_n_processes * 2);   // substract two ghost layers for evey grid
    sum += 2;    // add 1 ghost layer for top and bottom grid
    assert(sum == g_resolution);




    // 2D
    g_n_processes = 6;
    g_dim = DIM2;
    g_resolution = 13;

    // without ghost layers
    assert(local_grid_size({0,2}, false)[COORD_X] == 13);    
    assert(local_grid_size({1,2}, false)[COORD_X] == 12);
    assert(local_grid_size({0,2}, false)[COORD_Y] == 4);    assert(local_grid_size({1,2}, false)[COORD_Y] == 4);

    assert(local_grid_size({0,1}, false)[COORD_X] == 13);    assert(local_grid_size({1,1}, false)[COORD_X] == 12);
    assert(local_grid_size({0,1}, false)[COORD_Y] == 4);    assert(local_grid_size({1,1}, false)[COORD_Y] == 4);

    assert(local_grid_size({0,0}, false)[COORD_X] == 13);    assert(local_grid_size({1,0}, false)[COORD_X] == 12);
    assert(local_grid_size({0,0}, false)[COORD_Y] == 5);    assert(local_grid_size({1,0}, false)[COORD_Y] == 5);


    // with ghost layers
    assert(local_grid_size({0,2}, true)[COORD_X] == 13+1);    assert(local_grid_size({1,2}, true)[COORD_X] == 12+1);
    assert(local_grid_size({0,2}, true)[COORD_Y] == 4+1);    assert(local_grid_size({1,2}, true)[COORD_Y] == 4+1);

    assert(local_grid_size({0,1}, true)[COORD_X] == 13+1);    assert(local_grid_size({1,1}, true)[COORD_X] == 12+1);
    assert(local_grid_size({0,1}, true)[COORD_Y] == 4+2);    assert(local_grid_size({1,1}, true)[COORD_Y] == 4+2);

    assert(local_grid_size({0,0}, true)[COORD_X] == 13+1);    assert(local_grid_size({1,0}, true)[COORD_X] == 12+1);
    assert(local_grid_size({0,0}, true)[COORD_Y] == 5+1);    assert(local_grid_size({1,0}, true)[COORD_Y] == 5+1);
    
 
    

    cout << "OK" << endl;
    return SUCCESS;
}


int test_border_types()
{
    cout << "Testing border_types()..." << endl;
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

    cout << "OK" << endl;
    return SUCCESS;
}


int test_get_prime_factors()
{
    cout << "Testing get_prime_factors()..." << endl;
    int n;
    vector<int> n_factors;

    n = 2*3*5;
    n_factors = {2,3,5};
    assert(n_factors == get_prime_factors(n));

    n = 2*3*3*5*7*7*13;
    n_factors = {2,3,3,5,7,7,13};
    assert(n_factors == get_prime_factors(n));

    n = 13*17;
    n_factors = {13,17};
    assert(n_factors == get_prime_factors(n));

    n = 13;
    n_factors = {13};
    assert(n_factors == get_prime_factors(n));

    cout << "OK" << endl;

    return SUCCESS;
}

int test_split1D()
{
    cout << "Testing split1D()..." << endl;

    assert(split_1D(20, 5, 0, true) == 4);
    assert(split_1D(20, 5, 4, true) == 4);

    assert(split_1D(21, 5, 0, true) == 5);
    assert(split_1D(21, 5, 4, true) == 4);


    assert(split_1D(25, 2, 0, false) == 13);
    assert(split_1D(25, 2, 1, false) == 12);

    assert(split_1D(13, 3, 2, false) == 4);
    assert(split_1D(13, 3, 1, false) == 4);
    assert(split_1D(13, 3, 0, false) == 5);

    cout << "OK" << endl;
    return SUCCESS;
}


int test_to_global_grid_coords()

{
    cout << "to_global_grid_coords() ..." << endl;

    // 1D equal grid sizes
    g_n_processes = 3;
    g_dim = 1;
    g_resolution = 12;
    //                             rank   x  y
    assert(to_global_grid_coords({0, 0}, {1, 1})[1] == 1);
    assert(to_global_grid_coords({0, 1}, {1, 2})[1] == 6);
    assert(to_global_grid_coords({0, 2}, {1, 3})[1] == 11);
    assert(to_global_grid_coords({0, 0}, {1, 1})[0] == 1);
    assert(to_global_grid_coords({0, 1}, {1, 2})[0] == 1);
    assert(to_global_grid_coords({0, 2}, {1, 3})[0] == 1);


    // 1D different grid sizes
    g_n_processes = 3;
    g_dim = 1;
    g_resolution = 14;
    //                            rank   x  y
    assert(to_global_grid_coords({0, 0}, {1, 1})[1] == 1);
    assert(to_global_grid_coords({0, 1}, {1, 2})[1] == 7);
    assert(to_global_grid_coords({0, 2}, {1, 3})[1] == 13);
    assert(to_global_grid_coords({0, 0}, {1, 1})[0] == 1);
    assert(to_global_grid_coords({0, 1}, {1, 2})[0] == 1);
    assert(to_global_grid_coords({0, 2}, {1, 3})[0] == 1);

    cout << "OK" << endl;
    return SUCCESS;
}

int test_get_topo_shape()
{
    cout  << "Testing get_topo_shape() ..." << endl;

    g_n_processes = 6;
    g_dim = 1;
    assert(get_topo_shape()[COORD_X] == 1);    
    assert(get_topo_shape()[COORD_Y] == 6);

    g_n_processes = 5;
    g_dim = 1;
    assert(get_topo_shape()[COORD_X] == 1);    
    assert(get_topo_shape()[COORD_Y] == 5);


    g_n_processes = 6;
    g_dim = 2;

    assert(get_topo_shape()[COORD_X] == 2);
    assert(get_topo_shape()[COORD_Y] == 3);

    g_n_processes = 5;
    assert(get_topo_shape()[COORD_X] == 1);
    assert(get_topo_shape()[COORD_Y] == 5);

    g_n_processes = 12;
    assert(get_topo_shape()[COORD_X] == 2);
    assert(get_topo_shape()[COORD_Y] == 6);

    cout << "OK" << endl;
    return SUCCESS;
}


int main()
{
    cout << endl << "Starting unit tests..." << endl;

    assert(test_local_grid_size() == SUCCESS);
    assert(test_border_types() == SUCCESS);
    assert(test_get_prime_factors() == SUCCESS);
    assert(test_split1D() == SUCCESS);
    assert(test_to_global_grid_coords() == SUCCESS);
    assert(test_get_topo_shape() == SUCCESS);

    cout << "All Tests passed!" << endl;

    return SUCCESS;    
}