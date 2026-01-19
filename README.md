# Tight-Inclusion Continuous Collision Detection

[![Build](https://github.com/Continuous-Collision-Detection/Tight-Inclusion/actions/workflows/continuous.yml/badge.svg)](https://github.com/Continuous-Collision-Detection/Tight-Inclusion/actions/workflows/continuous.yml)

![](https://continuous-collision-detection.github.io/assets/tight-inclusion-teaser.png)

A conservative continuous collision detection (CCD) method with support for minimum separation.

To know more about this work, please read our ACM Transactions on Graphics paper:<br>
["A Large Scale Benchmark and an Inclusion-Based Algorithm for Continuous Collision Detection"](https://continuous-collision-detection.github.io/tight_inclusion/) and watch our [SIGGRAPH 2022 presentation](https://www.youtube.com/watch?v=7cRg52cWL8c).

## Installation via CMake

To compile the code, first, make sure CMake is installed.

To build the library on Linux or macOS:

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

Then you can run a CCD example:

```bash
./app/Tight_Inclusion_bin
```
## Installation via Conda
You can also install Tight-Inclusion via 
```sh
conda install conda-forge::tight-inclusion
```
if you have [Anaconda](https://www.anaconda.com/) installed on your machine. 

### Optional

We also provide an example that tests [sample queries](https://github.com/Continuous-Collision-Detection/Sample-Queries) using our CCD method. This requires installing `gmp` on your system before compiling the code. Set the CMake option `TIGHT_INCLUSION_WITH_SAMPLE_QUERIES` to `ON` when compiling:

```sh
cmake .. -DCMAKE_BUILD_TYPE=Release -DTIGHT_INCLUSION_WITH_SAMPLE_QUERIES=ON
make -j4
```

Then you can run `./app/Tight_Inclusion_bin` to test the handcrafted and simulation queries in the Sample Queries.

## Usage

### Overview

* Include: `#include <tight_inclusion/ccd.hpp>`
* Check vertex-face CCD: `bool ticcd::vertexFaceCCD(...)`
* Check edge-edge CCD: `bool ticcd::edgeEdgeCCD(...)`

### Details

:bulb: Each CCD function returns a boolean result corresponding to if a collision is detected. Because our method is *conservative*, we guarantee a result of `false` implies no collision occurs. If the result is `true`, there may not be a collision but we falsely report a collision. However, we can guarantee that this happens only if the minimal distance between the two primitives in this time step is no larger than `tolerance + ms + err` (see below for a description of these parameters).

#### Parameters

For both vertex-face and edge-edge CCD, the input query is given by eight vertices which are in the format of `Eigen::Vector3d`. Please read our code in `tight_inclusion/ccd.hpp` for the correct input order of the vertices.

Besides the input vertices, there are some input and output parameters for users to tune the performance or to get more information from the CCD.

Here is a list of the explanations of the parameters:

##### Input
* `err`: The numerical filters of the $x$, $y$ and $z$ coordinates. It measures the errors introduced by floating-point calculation when solving inclusion functions.
* `ms`: A minimum separation distance (no less than 0). We guarantee a collision will be reported if the distance between the two primitives is less than `ms`.
* `tolerance`: User-specific solving precision. It is the target maximal $x$, $y$, and $z$ length of the inclusion function. We suggest the to use `1e-6`.
* `t_max`: The time range $[0, t_{\max}]$ where we detect collisions. Since the input query implies the motion is in time interval $[0, 1]$, `t_max` should not be larger than 1.
* `max_itr`: The maximum number of iterations our inclusion-based root-finding algorithm can take. This enables early termination of the algorithm. If you set `max_itr < 0`, early termination will be disabled, but this may cause longer running times. We suggest setting `max_itr = 1e6`.
* `no_zero_toi`: For simulators which use non-zero minimum separation distance (`ms > 0`) to make sure intersection-free for each time-step, we have the option `no_zero_toi` to avoid returning a collision time `toi` of 0. The code will continue the refinement in higher precision if the output `toi` is 0 under the given `tolerance`, so the eventual `toi` will not be 0.
* `CCD_TYPE`: Enumeration of possible CCD schemes. The default and recommended type is `BREADTH_FIRST_SEARCH`. If set `DEPTH_FIRST_SEARCH`, the code will switch to a naive conservative CCD algorithm but lacks our advanced features.

##### Output

* `toi`: The time of impact. If multiple collisions happen in this time step, it will return the earliest collision time. If there is no collision, the returned `toi` value will be `std::numeric_limits<double>::infinity()`.
* `output_tolerance`: The resulting solve's precision. If early termination is enabled, the solving precision may not reach the target precision. This parameter will return the resulting solving precision when the code is terminated.

### Tips

:bulb: The input parameter `err` is crucial to guarantee our algorithm is a conservative method not affected by floating-point rounding errors. To run a single query, you can set `err = Eigen::Array3d(-1, -1, -1)` to enable a sub-function to calculate the real numerical filters when solving CCD. If you are integrating our CCD in simulators, you need to:

* Include the headler: `#include <tight_inclusion/interval_root_finder.hpp>`.
* Call
    ```
    std::array<double, 3> err_vf = ticcd::get_numerical_error()
    ```
    and
    ```
    std::array<double, 3> err_ee = ticcd::get_numerical_error()
    ```
* Use the parameter `err_ee` each time you call `bool ticcd::edgeEdgeCCD()` and `err_vf` when you call `bool ticcd::vertexFaceCCD()`.

The parameters for function `ticcd::get_numerical_error()` are:
* `vertices`: Vertices of the axis-aligned bounding box of the simulation scene. Before you run the simulation, you need to conservatively estimate the axis-aligned bounding box in which the meshes will be located during the whole simulation process, and the vertices should be the corners of the AABB.
* `is_vertex_face`: A boolean flag corresponding to if you are checking vertex-face or edge-edge CCD.
* `using_minimum_separation`: A boolean flag corresponding to if you are using minimum-separation CCD (the input parameter `ms > 0`).

To better understand or to get more details of our Tight-Inclusion CCD algorithm, please refer to our paper.

## Citation

If you use this work in your project, please consider citing the original paper:

```bibtex
@article{Wang:2021:Benchmark,
    title        = {A Large Scale Benchmark and an Inclusion-Based Algorithm for Continuous Collision Detection},
    author       = {Bolun Wang and Zachary Ferguson and Teseo Schneider and Xin Jiang and Marco Attene and Daniele Panozzo},
    year         = 2021,
    month        = oct,
    journal      = {ACM Transactions on Graphics},
    volume       = 40,
    number       = 5,
    articleno    = 188,
    numpages     = 16
}
```

