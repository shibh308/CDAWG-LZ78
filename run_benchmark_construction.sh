N=8
MAX=134217728
while [ $N -le $MAX ]
do
    $1 construct $2 $N
    N=$(( N * 2 ))
done
