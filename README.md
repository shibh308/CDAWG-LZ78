# CDAWG-LZ78

This repository contains the implementation and benchmark script for LZ78 substring compression on a CDAWG.

## Prepare the Text Files
You can use `./prepare_files.sh` to prepare the dataset.
This script will download three text files (`sources`, `dna`, and `english`) from the [Pizza&Chilli corpus](http://pizzachili.dcc.uchile.cl/), create fibonacchi string file (`fib`), and extract the first 128MiB of each text.

## Build

You can use [cmake](https://cmake.org/) to build the code.
This project requires [sdsl](https://github.com/simongog/sdsl-lite/tree/master).
You should install sdsl and set `SDSL_INCLUDE_DIR` and `SDSL_LIBRARY_DIR` in `CMakeLists.txt`
After completing the above steps, you can build the project using the following commands:
```
cmake .
make
```

## Run Benchmarks

To run all benchmarks, use the following script:
```
./run_benchmark_all.sh {path_to_executable}
```

To run only the compression or construction benchmarks, use the corresponding scripts:
```
./run_benchmark_compression.sh {path_to_executable} {filename} 
./run_benchmark_construction.sh {path_to_executable} {filename} 
```
