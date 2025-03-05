#!/bin/bash

# Exit immediately if a command fails
set -e

# Remove the "build" directory if it exists
rm -rf build

# Create a new "build" directory and navigate into it
mkdir build
cd build

# Run cmake and make
cmake ../
make -j4

# Return to the parent directory
cd ..
