# CDAWG-LZ78

## Download the text files
`./download_files.sh` downloads three text files (`sources`, `dna`, and `english`) from the [Pizza&Chilli corpus](http://pizzachili.dcc.uchile.cl/) and
extract the first 128MiB of each text.

## Compilation

CMake is available (detailed information: TODO)

## Benchmark

To run all benchmarks, use the following script:
```
./run_benchmark_all.sh {path_to_executable}
```

To run only the compression or construction benchmarks, use the corresponding script:
```
./run_benchmark_compression.sh {path_to_executable} {filename} 
./run_benchmark_construction.sh {path_to_executable} {filename} 
```
