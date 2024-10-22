#!/bin/bash

# This script specifies how to download dependencies for the software needed on Linux system.
# Note: I should test this on VM to see if they can be compatible to a linux system. 

# Dependencies for FastP: 
git clone https://github.com/intel/isa-l.git
cd isa-l
./autogen.sh
./configure --prefix=/usr --libdir=/usr/lib64
make -j
sudo make install
cd ../

git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate
cmake -B build
cmake --build build
cmake --install build
cd ../
