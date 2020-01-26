#!/bin/bash

for ((i = 1; i <= $1; i++ ));
do
    ./../build/sim/heatmaps config.json expt$i 0 $i &
done
