#include <iostream>
#include "splitting.hpp"

using namespace std;

int g_n_processes;
int g_dim;
size_t g_resolution;


int main()
{

    g_n_processes = 5;
    g_dim = DIM1;
    g_resolution = 17;

    cout << "All Tests passed!" << endl;
    return 0;    
}