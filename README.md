# ArtificialGeneNets
Using supervised artificial neural nets to reverse-engineer gene regulatory networks

Requires morphologica

build in the usual way:
mkdir build
cd build
cmake ..
make
cd ..


To initalize (adds highres.h5 and network.h5 to configs):
cd pyscripts
python highres.py
python netff.py
cd ..

Then to run...

...training
./build/sim/model configs/config.json logs 0 1

...testing
./build/sim/model configs/config.json logs 1 1
