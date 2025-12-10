#!/bin/bash

mkdir build
cd build
cmake ..
make 

cd ..
sudo cp build/imsonpro /usr/local/bin/
sudo cp -r build /usr/local/bin/imsonpro_build