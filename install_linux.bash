#!/bin/bash

sudo apt-get update
sudo apt-get install -y mpg123
sudo apt-get install -y imagemagick
sudo apt-get install -y gnuplot
sudo apt-get install -y doxygen graphviz

cd external
git clone https://github.com/dstahlke/gnuplot-iostream.git
