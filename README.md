# Tight-Inclusion Continuous Collision Detection 
A conservative CCD method with support for minimum separation.

[![Build](https://github.com/Continuous-Collision-Detection/Tight-Inclusion/actions/workflows/continuous.yml/badge.svg)](https://github.com/Continuous-Collision-Detection/Tight-Inclusion/actions/workflows/continuous.yml)

You can read more about this work in our ACM Transactions on Graphics paper:

["A Large Scale Benchmark and an Inclusion-Based Algorithm forContinuous Collision Detection"](https://continuous-collision-detection.github.io/)

---

## Compiling Instruction 

To compile the code, first make sure CMake is installed. 

To build the library on Linux or macOS:
```sh
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make -j4
```

---
 
## Usage
Include `#include <tight_inclusion/inclusion_ccd.hpp>`

To check edge-edge ccd, use `bool edgeEdgeCCD_double();`

To check vertex-face ccd, use `bool vertexFaceCCD_double();`

Please read the annotations for the details of the parameter setting

