# requirements on ubuntu
# sudo apt-get build-essentials
# sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# required modules on cluster
# module load mpi/openmpi-x86_64
# module load pmi/pmix-x86_64


# directories
OUT_DIR = ./out
BUILD_DIR = $(OUT_DIR)/build
LOG_DIR = $(OUT_DIR)/logs
SRC_DIRS = ./src
TEST_DIR = ./test

HEADERS = $(shell find $(SRC_DIRS) -name *.hpp -or -name *.h)


# create directories if the dont exist
$(shell mkdir -p $(OUT_DIR))
$(shell mkdir -p $(BUILD_DIR))
$(shell mkdir -p $(LOG_DIR))
$(shell mkdir -p $(TEST_DIR)/out)

# compiler settings
CXX=g++
MPICXX?=mpic++

# Flags for debugging
CXXFLAGS := -std=c++11 -O0 -g -pedantic -march=native -ffast-math -Wall
XXFLAGS += -ggdb3	# enable valgrind to show code

# Flags for release
#CXXFLAGS := -std=c++11 -O3 -pedantic -march=native -ffast-math -Wall

all:  jacobiMPI jacobiMPI_float jacobiSERIAL run test


# Serial Baseline for comparision
jacobiSERIAL: $(SRC_DIRS)/main.cpp $(HEADERS)
	$(CXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiSERIAL $(CXXFLAGS)

# Task 2
jacobiMPI: $(SRC_DIRS)/main.cpp $(HEADERS)
	$(MPICXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiMPI -lpthread -DUSEMPI $(CXXFLAGS)

# Task 3
jacobiMPI_float: $(SRC_DIRS)/main.cpp $(HEADERS)
	$(MPICXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiMPI_float -lpthread -DUSEMPI -DUSE_FLOAT $(CXXFLAGS)


test_splitting: $(TEST_DIR)/test_splitting.cpp $(HEADERS)
	$(CXX) $(TEST_DIR)/test_splitting.cpp -I$(SRC_DIRS) -o $(TEST_DIR)/out/test_splitting $(CXXFLAGS)
	$(TEST_DIR)/out/test_splitting

test_mpi: $(TEST_DIR)/test_mpi.cpp $(HEADERS)
	$(MPICXX) $(TEST_DIR)/test_mpi.cpp -I$(SRC_DIRS) -o $(TEST_DIR)/out/test_mpi $(CXXFLAGS) -lpthread -DUSEMPI
	mpirun -np 6 --oversubscribe $(TEST_DIR)/out/test_mpi

test: test_splitting test_mpi	# Unit tests

# debugging presets for gdb
GDB_ARGS = -ex "break solver.hpp:388" 
GDB_RUN = -ex "run 4 10"

# quick run to check if all programs work properly
run:
#./$(BUILD_DIR)/jacobiSERIAL 100 10
	mpirun -np 4 $(BUILD_DIR)/jacobiMPI 20 10
#mpirun -np 4 $(BUILD_DIR)/jacobiMPI_float 100 10
#mpirun -np 2 xterm -e gdb $(GDB_ARGS) $(GDB_RUN) $(BUILD_DIR)/jacobiMPI 10 10	# debugging
#python3 ./scripts/plot_solution.py	# execute a small script to display a plot
#mpirun -np 2 $(BUILD_DIR)/jacobiMPI_float 2D 100 10

clean:
	rm -rf $(OUT_DIR)
	rm -rf $(TEST_DIR)/out
