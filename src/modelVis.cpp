#include "NetsVis.h"

int main (int argc, char **argv){
    if(argc<4){ cout<<"supply logfile trials(or 0 for display) seed"<<endl<<flush; return 0; }
    Net N(argv[1]);
    int T = stoi(argv[2]);
    srand(stoi(argv[3]));
    if(T<1){ // TESTING
        N.loadWeights();
        N.plotMaps();
        N.plotResponses();
    } else { // TRAINING
        N.run(T,1000,100,false);
        N.saveOutputs();
        N.saveWeights();
    }
    return 0;
}
