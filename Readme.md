# Distance transform computation using the fast sweeping method.

See `distance.h` for a complete documentation.

## Reference ##
 * "A fast sweeping method for Eikonal equations"
   by H. Zhao, Mathematics of computation, 74 (2005),pp. 603–627
   https://www.math.uci.edu/~zhao/homepage/research_files/FSM.pdf
 * "Finding the Skeleton of 2D Shape and Contours: Implementation of Hamilton-Jacobi Skeleton"
   by Yuchen He, Sung Ha Kang, Luis Álvarez (2020)
   https://www.ipol.im/pub/pre/296/preprint.pdf

## Example ##
The `test` directory contains the source code of a small program that generates a distance map from an input image.

<img src="test/data/a.png" width="256px"/> <img src="test/result/a.png" width="256px"/> 

<img src="test/data/hello.png" width="256px"/> <img src="test/result/hello.png" width="256px"/> 

Image reading and writing libraries [stb_image.h, stb_image_write.h](https://github.com/nothings/stb/) by Sean Barrett (public domain).

## Build ##

A CMake configuration file is provided in order to build a static library and
the associated documentation.
A typical usage of CMake may be:
```bash
mkdir build
cd build
cmake ..
cmake --build .
```

## License ##
This project is licensed under the MIT License, see LICENSE for more information.
