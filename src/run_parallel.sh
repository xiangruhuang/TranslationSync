#!/bin/bash
mkdir -p logs
for i in `seq 1 $1`; do
    j=$((${i}-1))
    ./run_graph.sh ${j} $1 &> logs/log_${j}_$1 &
done
