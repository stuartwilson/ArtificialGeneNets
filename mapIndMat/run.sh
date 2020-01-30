#!/bin/bash

dst=data/expt_

for ((i = 1; i <= $1; i++ ));
do
    mkdir $dst$i
    cp inputs.h5 $dst$i/inputs.h5
    cp network.h5 $dst$i/network.h5
    ./../build/sim/greig config.json $dst$i $2 $i &
done
