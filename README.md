# Tight-Inclusion
Tight-Inclusion Continuous Collision Detection 

:warning: CAUTION: The Code Is Tested Using GCC, Probably Also Works On MSVC and Clang

:warning: The folder /data/ under root folder contains a subset of the whole dataset.

:warning: The output information contains: if testing handcrafted/simulation dataset, if testing edge-edge or vertex-face queries, the number of queries, number of false positives, number of false negatives, average running time, percentage of early return, maximum solving tolerance delta and average tolerance.



---

## Compiling Instruction 

To compile the code, first you need to install 
* CMake (https://cmake.org/), 
* GMP (https://gmplib.org/) 
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
 
##	Terminal commands for testing queries:
```sh
<PATH_TO_BINARY>/TI_CCD_bin  
```
**Note that:** Due to the limits of upload file size, we only provide a subset of queries
