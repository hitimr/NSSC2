#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include "common.hpp"

using namespace std;

extern int g_my_rank;

template<typename T>
T convertTo(const int position, const T init, int argc, char *argv[]) {
  if (argc <= position) {
    if(g_my_rank == MASTER)
    std::cout
        << "Conversion of argument " << position
        << " to 'int' failed, not enough parameters, using default parameter: "
        << init << std::endl;
    return init;
  }
  T arg;
  std::istringstream tmp(argv[position]);
  tmp >> arg;
  return arg;
}

int convertToDim(char *arg)
{
  string dim_string = arg;
  int dim_int;

  dim_string.erase(1,1);  // remove the D
  dim_int = stoi(dim_string);  // transform the remaining number to an integer

  return dim_int;
}