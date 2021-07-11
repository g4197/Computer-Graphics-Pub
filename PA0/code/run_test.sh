#!/usr/bin/env bash

# If project not ready, generate cmake file.
if [[ ! -d build ]]; then
    mkdir -p build
    cd build
    cmake ..
    cd ..
fi

# Build project.
cd build
make -j
cd ..

# Run all testcases.
# You can comment some lines to disable the run of specific examples.
mkdir -p output

for file in `ls testcases/`
do
    bin/PA0 testcases/$file output/$file.bmp
done
