N=8
MAX=$(cat ../data/$1 | wc -c)

# NがMAX以下の間、ループを回す
while [ $N -le $MAX ]
do
    # コマンドの実行
    ./CDAWG_LZ78 $1 $N
    
    # Nを倍にする
    N=$(( N * 2 ))
done

./CDAWG_LZ78 $1 $MAX
