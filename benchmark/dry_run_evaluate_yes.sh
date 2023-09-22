
for i in $(seq 1 1 100); do
    python3 evaluate.py -o 2 -i $i -r yes -d count
    # python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r no -d relab
    python3 evaluate.py -o 2 -i $i -r yes -d relab

done
