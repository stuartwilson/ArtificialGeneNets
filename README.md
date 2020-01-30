# ArtificialGeneNets
Using supervised artificial neural nets to reverse-engineer gene regulatory networks

Requires morphologica

build in the usual way:
mkdir build
cd build
cmake ..
make
cd ..

Then to run e.g., mapIndMat:
    cd mapIndMap
    python heatmaps.py 
    python netspec.py
    bash run.sh 1 1000000
    bash test.sh 1 0
    python analyse.py 1 0
