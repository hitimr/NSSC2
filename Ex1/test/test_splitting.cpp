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


    // check serial case
    g_n_processes = 1;
    g_dim = DIM1;
    g_resolution = 8;
    assert(local_grid_size(0)[COORD_X] == g_resolution);
    assert(local_grid_size(0)[COORD_Y] == g_resolution);



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

    cout << "OK" << endl;
    return SUCCESS;
}


int test_border_types()
{
    cout << "Testing border_types()" << endl;
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


int main()
{
    cout << endl << "Starting unit tests..." << endl;

    assert(test_local_grid_size() == SUCCESS);
    assert(test_border_types() == SUCCESS);
    assert(test_get_prime_factors() == SUCCESS);

    cout << "All Tests passed!" << endl;

    return SUCCESS;    
}