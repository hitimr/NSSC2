#include <iostream>
#include <assert.h>
#include "splitting.hpp"

using namespace std;

int g_n_processes = -1;
int g_dim = -1;
size_t g_resolution = 0;


int test_grid_size()
{
    g_n_processes = 5;
    g_dim = DIM1;
    g_resolution = 17;

    for(int rank = 0; rank < g_n_processes; rank++)
    {
        cout << grid_size(rank)[1] << endl;
    }
    return 0;
}


int main()
{

    assert(test_grid_size() == 0);
    cout << "All Tests passed!" << endl;
    return 0;    
}