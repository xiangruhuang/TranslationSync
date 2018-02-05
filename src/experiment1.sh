#!/bin/bash
#graph type 1
for ratio in `seq 0 9`; do
    mkdir -p results/graph1_uniform_ratio0${ratio}
    for eid in `seq 1 100`; do
        if [ ! -f "results/graph1_uniform_ratio0${ratio}/${eid}" ]; then
            make synthetic n=2000 edge_density=0.2 bias=1.0 inc=0.0 \
                noise_type=2 noise_ratio=0.${ratio} decay=0.9 \
                output=results/graph1_uniform_ratio0${ratio}/${eid}
        fi
    done
done

#graph type 2
for ratio in `seq 0 9`; do
    mkdir -p results/graph2_uniform_ratio0${ratio}
    for eid in `seq 1 100`; do
        if [ ! -f "results/graph2_uniform_ratio0${ratio}/${eid}" ]; then
            make synthetic n=2000 edge_density=0.2 bias=0.05 inc=0.15 \
                noise_type=2 noise_ratio=0.${ratio} decay=0.9 \
                output=results/graph2_uniform_ratio0${ratio}/${eid}
        fi
    done
done

#graph type 3
for ratio in `seq 0 9`; do
    mkdir -p results/graph3_uniform_ratio0${ratio}
    for eid in `seq 1 100`; do
        if [ ! -f "results/graph3_uniform_ratio0${ratio}/${eid}" ]; then
            make synthetic n=20000 edge_density=0.002 bias=1.0 inc=0.0 \
                noise_type=2 noise_ratio=0.${ratio} decay=0.9 \
                output=results/graph3_uniform_ratio0${ratio}/${eid}
        fi
    done
done

#graph type 4
for ratio in `seq 0 9`; do
    mkdir -p results/graph4_uniform_ratio0${ratio}
    for eid in `seq 1 100`; do
        if [ ! -f "results/graph4_uniform_ratio0${ratio}/${eid}" ]; then
            make synthetic n=20000 edge_density=0.2 bias=0.005 inc=0.015 \
                noise_type=2  noise_ratio=0.${ratio} decay=0.9 \
                output=results/graph4_uniform_ratio0${ratio}/${eid}
        fi
    done
done
