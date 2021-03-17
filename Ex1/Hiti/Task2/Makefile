# requirements on ubuntu
# sudo apt-get build-essentials
# sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# required modules on cluster
# module load mpi/openmpi-x86_64
# module load pmi/pmix-x86_64

CXX=g++
MPICXX?=mpic++
CXXFLAGS := $(CXXFLAGS) -std=c++14 -O3 -Wall -pedantic -march=native -ffast-math

BUILD_DIR ?= ./out
SRC_DIRS ?= ./src
HEADERS = $(shell find $(SRC_DIRS) -name *.hpp -or -name *.h)
$(shell mkdir -p $(BUILD_DIR))

all: jacobiSERIAL jacobiMPI

jacobiSERIAL: Makefile $(SRC_DIRS)/main.cpp $(HEADERS)
	$(CXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiSERIAL $(CXXFLAGS)

jacobiMPI: Makefile $(SRC_DIRS)/main.cpp $(HEADERS)
	$(MPICXX) $(SRC_DIRS)/main.cpp -o $(BUILD_DIR)/jacobiMPI -lpthread -DUSEMPI $(CXXFLAGS)

clean:
	rm $(BUILD_DIR)/*
