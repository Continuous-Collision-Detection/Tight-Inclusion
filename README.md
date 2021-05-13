# Tight-Inclusion Continuous Collision Detection 
![](./fig/line-search.jpg)
A conservative Continuous Collision Detection (CCD) method with support for minimum separation.

[![Build](https://github.com/Continuous-Collision-Detection/Tight-Inclusion/actions/workflows/continuous.yml/badge.svg)](https://github.com/Continuous-Collision-Detection/Tight-Inclusion/actions/workflows/continuous.yml)

You can read more about this work in our ACM Transactions on Graphics paper:

["A Large Scale Benchmark and an Inclusion-Based Algorithm forContinuous Collision Detection"](https://continuous-collision-detection.github.io/)
```bash
@article{Wang:2021:Benchmark,
    title   = {A Large Scale Benchmark and an Inclusion-Based Algorithm for Continuous Collision Detection},
    author  = {Bolun Wang and Zachary Ferguson and Teseo Schneider and Xin Jiang and Marco Attene and Daniele Panozzo},
    year    = 2021,
    journal = {ACM Transactions on Graphics}
}
```
## Compiling Instruction 

To compile the code, first make sure CMake is installed. 

To build the library on Linux or macOS:
```sh
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make
```
Then you can run a CCD example:
```bash
./Tight_Inclusion_bin 
```
We also show you an example to run the [Sample Queries](https://github.com/Continuous-Collision-Detection/Sample-Queries) using our CCD method. You may need to install `gmp` before compiling the code. When compiling the code, you need to set the CMake option `TIGHT_INCLUSION_WITH_TESTS` as `ON`:
```sh
cmake ../ -DCMAKE_BUILD_TYPE=Release -DTIGHT_INCLUSION_WITH_TESTS=ON
make
```
Then you can run `./Tight_Inclusion_bin` to test the `handcrafted queries` in the Sample Queries.
## Usage
Include `#include <tight_inclusion/inclusion_ccd.hpp>`

To check edge-edge ccd, use `bool edgeEdgeCCD_double();`

To check vertex-face ccd, use `bool vertexFaceCCD_double();`
💡 If collision is detected, the ccd function will return `true`, otherwise, the ccd function will return `false`. Since our method is conservative, if the returned result is `false`, we guarantee that there is no collision happens. If the result is `true`, it is possible that there is no collision but we falsely report a collision, but we can guarantee that this happens only if the minimal distance between the two premitives in this time step is no larger than `tolerance + ms + err`. We wil explain these parameters below.  
For both edge-edge ccd and vertex-face ccd, the input CCD query is presented by 8 vertices which are in the format of `Eigen::Vector3d`. Please read our code in `tight_inclusion/inclusion_ccd.hpp` for the correct input order of the vertices. 
Beside the input vertices, there are plenty of input and out parameters for users to  

Please read the annotations for the details of the parameter setting

