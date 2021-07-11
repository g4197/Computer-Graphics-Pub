# If project not ready, generate cmake file.
if [[ ! -d build_linux ]]; then
    mkdir -p build_linux
    cd build_linux
    cmake ..
    cd ..
else
    rm -r build_linux
    mkdir -p build_linux
    cd build_linux
    cmake ..
    cd ..
fi

# build_linux project.
cd build_linux
make -j
cd ..

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.
mkdir -p output
time bin/FINAL testcases/scene02_office.txt output/scene03_final.bmp 400000
# time bin/FINAL testcases/scene02_depth.txt output/scene03_depth_final.bmp 40000
# bin/PA1 testcases/scene02_cube.txt output/scene02.bmp
# bin/PA1 testcases/scene03_sphere.txt output/scene03.bmp
# bin/PA1 testcases/scene04_axes.txt output/scene04.bmp
# bin/PA1 testcases/scene05_bunny_200.txt output/scene05.bmp
# bin/PA1 testcases/scene06_bunny_1k.txt output/scene06.bmp
# bin/PA1 testcases/scene07_shine.txt output/scene07.bmp

