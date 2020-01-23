#!/bin/bash

for ((i = 1; i <= $1; i++ ));
do
    mkdir expt$i
    cp ../configs/inputs.h5 expt$i/inputs.h5
    cp ../configs/network.h5 expt$i/network.h5

    ./../build/sim/model ../configs/config.json expt$i $2 $i &
done
