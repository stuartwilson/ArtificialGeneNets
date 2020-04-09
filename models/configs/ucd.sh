#!/bin/bash
FILES=data/ucd/*

for f in $FILES
do
python ucd_maps.py $f
done
