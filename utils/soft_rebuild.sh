#!/bin/bash

cd build
make 

cd ..
sudo rm -rf /usr/local/bin/imsonpro_build
sudo rm -rf /usr/local/bin/imsonpro
sudo cp build/imsonpro /usr/local/bin/
sudo cp -r build /usr/local/bin/imsonpro_build