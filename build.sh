#!/bin/bash

mkdir build
cd build
cmake ..
make 

sudo rm -rf /usr/local/bin/imsonpro
sudo rm -rf /usr/local/bin/imsonpro_build

cd ..
sudo cp build/imsonpro /usr/local/bin/
sudo cp -r build /usr/local/bin/imsonpro_build