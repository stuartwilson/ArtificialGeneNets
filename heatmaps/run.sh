#!/bin/bash

for ((i = 1; i <= $1; i++ ));
do
    mkdir expt$i
    cp inputs.h5 expt$i/inputs.h5
    cp network.h5 expt$i/network.h5

    ./../build/sim/heatmaps config.json expt$i $2 $i &
done
