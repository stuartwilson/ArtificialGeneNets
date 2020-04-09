#include "NetsVis.h"
int main (int argc, char **argv){
    if(argc<2){ cout<<"supply json"<<endl<<flush; return 0; }
    Net N(argv[1]);
    N.plotMaps();
    return 0;
}
