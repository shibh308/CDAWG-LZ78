for filename in "sources" "dna" "english"
do
  ./run_benchmark_compression.sh $1 $filename
  ./run_benchmark_construction.sh $1 $filename
done

