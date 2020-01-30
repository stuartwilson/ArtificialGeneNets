#!/bin/bash
dst=data/expt_
for ((i = 1; i <= $1; i++ ));
do
    ./../build/sim/greig config.json $dst$i 0 $i &
done
