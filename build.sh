nvcc main.cpp -I/usr/local/include -std=c++20 -x cu -arch=sm_80 -O3  -Xcompiler='-fopenmp' -o route
