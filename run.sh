N=8
MAX=$(cat ./data/$1.txt | wc -c)

# NがMAX以下の間、ループを回す
while [ $N -le $MAX ]
do
    # コマンドの実行
    ./a.out $1 $N
    
    # Nを倍にする
    N=$(( N * 2 ))
done

./a.out $1 $MAX

