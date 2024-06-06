for filename in "sources" "dna" "english" "fib"
do
  ./run_benchmark_compression_cdawg.sh $1 $filename
  ./run_benchmark_compression_suffixtree.sh $1 $filename
  ./run_benchmark_construction.sh $1 $filename
done

