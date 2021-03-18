# requirements on ubuntu
# sudo apt-get build-essentials
# sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# required modules on cluster
# module load mpi/openmpi-x86_64
# module load pmi/pmix-x86_64


# directories
OUT_DIR ?= ./out
BUILD_DIR ?= $(OUT_DIR)/build
LOG_DIR ?= $(OUT_DIR)/logs
SRC_DIRS ?= ./src
HEADERS = $(shell find $(SRC_DIRS) -name *.hpp -or -name *.h)

# create directories if the dont exist
$(shell mkdir -p $(OUT_DIR))
$(shell mkdir -p $(BUILD_DIR))
$(shell mkdir -p $(LOG_DIR))

# compiler settings
CXX=g++
MPICXX?=mpic++
CXXFLAGS := $(CXXFLAGS) -std=c++14 -O0 -g -Wall -pedantic -march=native -ffast-math


all: jacobiSERIAL jacobiMPI jacobiMPI_float run

# Serial Baseline for comparision
jacobiSERIAL: Makefile $(SRC_DIRS)/main.cpp $(HEADERS)
	$(CXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiSERIAL $(CXXFLAGS)

# Task 2
jacobiMPI: Makefile $(SRC_DIRS)/main.cpp $(HEADERS)
	$(MPICXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiMPI -lpthread -DUSEMPI $(CXXFLAGS)

# Task 3
jacobiMPI_float: Makefile $(SRC_DIRS)/main.cpp $(HEADERS)
	$(MPICXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiMPI_float -lpthread -DUSEMPI -DUSE_FLOAT $(CXXFLAGS)

# quick run to check if all programs work properly
run:
	./$(BUILD_DIR)/jacobiSERIAL 100 10
	mpirun -np 2 $(BUILD_DIR)/jacobiMPI 1D 100 10

clean:
	rm -rf $(OUT_DIR)
