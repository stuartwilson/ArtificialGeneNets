# ArtificialGeneNets

Using supervised artificial neural nets to reverse-engineer gene regulatory networks

## Install

Requires morphologica (https://github.com/ABRG-Models/morphologica)

build in the usual way:
mkdir build
cd build
cmake ..
make
cd ..

## Overview

This project is essentially an implementation of the recurrent backpropagation algorithm for tranining neural networks, in a supervised manner to map a given set of inputs to a given set of outputs. The algorithm has the advantage over the commonly used backpropagation algorithm that it can be used to train networks with any topolgy (backprop is restricted to feed-forward networks only). Networks with recurrent connectivity can display attractor dynamics, hence it is important to allow the dynamics of recurrent networks to settle before connection weights are adjusted, and to detect when the dynamics have failed to settle. The resulting algorithm was introduced by Pineda in the following journal article:

Pineda, FJ. (1987) Generalization of back-propagation to recurrent neural networks. Physical Review Letters, 59, 2229.

The core algorithm, as reported by Pineda (1987) is implemented in RecurrentNetwork.h, which defines a RecurrentNetwork class, the methods of which can be called to iterate the activation and learning dynamics in response to pairs of inputs and target outputs. 

In addition RecurrentNetworkModel.h provides a useful interface to these methods, in order that inputs and target outputs can be imported from external files, and that the activation and learning can be easily visualised. An instance of class RecurrentNetworkModel provides access to:

i)      A vector of Map objects, each of which contains input/target pairs imported from external .h5 files; 
ii)     A Domain object, which contains a 2D hexagonal lattice that can be used to image the network activation;
iii)    The identity of network nodes that will be used as Input or Output nodes, or Context nodes (nodes that take specific input values when a given Map object is used for training/testing).

## Example



 


  
