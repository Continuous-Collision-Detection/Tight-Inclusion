# Tight-Inclusion
Tight-Inclusion Continuous Collision Detection 

:warning: CAUTION: The Code Is Tested Using GCC, Probably Also Works On MSVC and Clang


---

## Compiling Instruction 

To compile the code, first you need to install 
* CMake (https://cmake.org/), 
* GMP (https://gmplib.org/) (GMP is required only when the CMake option TIGHT_INCLUSION_WITH_GMP is setted as ON)

in your system. 

To build the executable file, you can use CMake
```sh
cd Tight-Inclusion/
mkdir build
cd build
cmake ../  -DCMAKE_BUILD_TYPE=Release
make
```


---
 
##	Usage:
```sh
<PATH_TO_BINARY>/TI_CCD_bin  
```
**Note that:** Due to the limits of upload file size, we only provide a subset of queries

## Usage
Include `#include <tight_inclusion/inclusion_ccd.hpp>`
To check edge-edge ccd, use `bool edgeEdgeCCD_double();`
To check vertex-face ccd, use `bool vertexFaceCCD_double();`
Please read the annotations for the details of the parameter setting


  Call one of the `is_outside` function with a triangle, point, or segment.
