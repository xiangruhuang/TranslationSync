#!/bin/bash
#graph type 1
for ratio in `seq $1 $2 99`; do
    if [ ! -d "resampling/graph1_uniform_ratio${ratio}" ]; then
        mkdir -p resampling/graph1_uniform_ratio${ratio}
        make synthetic n=2000 edge_density=0.2 bias=1.0 inc=0.0 \
            noise_type=2 noise_ratio=0.${ratio} decay=0.9 \
            output=resampling/graph1_uniform_ratio${ratio}
    fi
done

#graph type 2
for ratio in `seq $1 $2 99`; do
    if [ ! -d "resampling/graph2_uniform_ratio${ratio}" ]; then
        mkdir -p resampling/graph2_uniform_ratio${ratio}
        make synthetic n=2000 edge_density=0.2 bias=0.1 inc=0.3 \
            noise_type=2 noise_ratio=0.${ratio} decay=0.9 \
            output=resampling/graph2_uniform_ratio${ratio}
    fi
done

#graph type 3
for ratio in `seq $1 $2 99`; do
    if [ ! -d "resampling/graph3_uniform_ratio${ratio}" ]; then
        mkdir -p resampling/graph3_uniform_ratio${ratio}
        make synthetic n=20000 edge_density=0.002 bias=1.0 inc=0.0 \
            noise_type=2 noise_ratio=0.${ratio} decay=0.9 \
            output=resampling/graph3_uniform_ratio${ratio}
    fi
done

#graph type 4
for ratio in `seq $1 $2 99`; do
    if [ ! -d "resampling/graph4_uniform_ratio${ratio}" ]; then
        mkdir -p resampling/graph4_uniform_ratio${ratio}
        make synthetic n=20000 edge_density=0.2 bias=0.05 inc=0.15 \
            noise_type=2  noise_ratio=0.${ratio} decay=0.9 \
            output=resampling/graph4_uniform_ratio${ratio}
    fi
done
