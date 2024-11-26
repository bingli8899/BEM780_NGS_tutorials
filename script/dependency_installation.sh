#!/bin/bash

# This script specifies how to download dependencies for the software needed on Linux system.
# Note: I should test this on VM to see if they can be compatible to a linux system. 

# To run any softwares, several dependencies should be installed and it depends on the operating system of the users. On my MacBook, to download and install dependencies, download homebrew from here: https://brew.sh/ first and then run the below codes on terminals: 

# gcc --version # Check the gcc version. Need to be version 4.8 or later. If not, check gcc website and download the latest version  
# brew install isa-l
# brew install libdeflate
 
# Note: Since MAC could easily install conda (method 1), so I will change the above methods to focus on Linux system without using homebrew. This requires me to change the software_installation.sh later.  

# On Linux system, run the codes below, run the below code first in the root folder that the users want to store the dependencies. The below script should be revised to the directory the user wants to install the dependencies. 

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


