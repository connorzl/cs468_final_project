# cs468_final_project
Built using [geometry-central](http://geometry-central.net/) and [Polyscope](http://polyscope.run/).

### How to build

**Unix-like machines**: configure (with cmake) and compile
```
cd cs468_final_project
mkdir build
cd build
cmake ..
make -j6
```

### How to run
```
./bin/cs468_proj /path/to/a/mesh
```

### Code description

`src/main.cpp` performs mesh I/O and visualization, calling the two scripts below. 
`src/vector_heat.cpp` implements scalar diffusion and vector heat parallel transport.
`src/weiszfeld.cpp` implements the Weiszfeld algorithm for computing geometric medians.
