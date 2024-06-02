# CDAWG-LZ78

This repository contains the implementation and benchmark script for LZ78 substring compression on a CDAWG.

## Download the Text Files
Run `./download_files.sh` to download three text files (`sources`, `dna`, and `english`) from the [Pizza&Chilli corpus](http://pizzachili.dcc.uchile.cl/) and
extract the first 128MiB of each text.

## Build
You can use [cmake](https://cmake.org/) to build the code.
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
