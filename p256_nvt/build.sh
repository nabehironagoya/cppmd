#!/bin/bash
rm -r build
mkdir -p build
cd build
/home/hwatanabe/app/cmake-2.8.12.2/bin/cmake ../ -DOPENMP=OFF
make 

cp ../sample.gro .
