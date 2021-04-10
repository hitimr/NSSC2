#pragma once

#ifdef USEMPI
	#include <mpi.h>
#endif

// ---- Defines
// Dimensions
#define DIM1 1
#define DIM2 2
#define COORD_X 0
#define COORD_Y 1

// Border Types
// Assignemnt resembles the enumerator from solver.hpp,
// enum Cell { UNKNOWN = 0, DIR = 1, NEU = 2, ROB = 0 };
#define BORDER_UNKNOWN 0
#define BORDER_DIR     1    // Dirichlet
#define BORDER_NEU     2    // Neumann
#define BORDER_ROB     0    // Robin
#define BORDER_GHOST   3    // Ghost Layer

// Directional Assignemtns
#define BOTTOM  0
#define RIGHT   1
#define TOP     2
#define LEFT    3


// Misc
#define MASTER 0
#define DEBUG_RANK -1   // just for debugging
#define NO_NEIGHBOUR -2
#define PRINT_LOCAL_SOLUTIONS false
#define CHECKPOINT std::cout << "Rank " << g_my_rank << " made it here!" << endl;

// ---- Globals
// globals according to command line
extern int g_n_processes;
extern int g_dim;
extern size_t g_resolution;
extern int g_iterations;
extern int g_my_rank;
#ifdef USEMPI
    extern MPI_Comm g_topo_com;
#endif


#include <string>
#include <array>
// taken from https://dev.to/aggsol/calling-shell-commands-from-c-8ej
std::string execCommand(const std::string cmd, int& out_exitStatus)
{
    out_exitStatus = 0;
    auto pPipe = ::popen(cmd.c_str(), "r");
    if(pPipe == nullptr)
    {
        throw std::runtime_error("Cannot open pipe");
    }

    std::array<char, 256> buffer;

    std::string result;

    while(not std::feof(pPipe))
    {
        auto bytes = std::fread(buffer.data(), 1, buffer.size(), pPipe);
        result.append(buffer.data(), bytes);
    }

    auto rc = ::pclose(pPipe);

    if(WIFEXITED(rc))
    {
        out_exitStatus = WEXITSTATUS(rc);
    }

    return result;
}
    