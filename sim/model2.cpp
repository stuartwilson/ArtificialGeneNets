#include "Nets.h"
#include <iostream>
using namespace std;


class Areas: public Net{

public:
    int nMaps, nLocations, nIns;
    vector<vector<double> > Ins;
    vector<vector<int> > Maps;

    int nPatterns, mapIndex;
    vector<vector<double> > Outs;
    vector<vector<vector<int> > > LocationsByField;
    vector<int> knockoutID;
    vector<double> X, Y;

    Areas(string logpath, string nname) : Net(logpath){

        float dt = 1.0;
        float tauW = 2.0;
        float tauX = 1.0;
        float tauY = 1.0;
        float weightNudgeSize = 0.001;
        float divergenceThreshold = 0.000001;
        int maxConvergenceSteps = 400;

        HdfData network(nname,1);
        network.read_contained_vals ("x", X);
        network.read_contained_vals ("y", Y);
        nLocations = X.size();

        vector<int> Ntmp(0);
        network.read_contained_vals ("N", Ntmp);
        int N = Ntmp[0];
        vector<int> inputID, outputID, contextID, pre, post;
        network.read_contained_vals ("inputs", inputID);
        network.read_contained_vals ("outputs", outputID);
        network.read_contained_vals ("knockouts", knockoutID);
        network.read_contained_vals ("pre", pre);
        network.read_contained_vals ("post", post);

        { // Get input values (co-ordinates)
            vector<double> tmp;
            network.read_contained_vals ("inPatterns", tmp);
            nIns = tmp.size()/nLocations;
            Ins.resize(nIns,vector<double>(nLocations));
            int k=0;
            for(int i=0;i<nIns;i++){
                for(int j=0;j<nLocations;j++){
                    Ins[i][j] = tmp[k];
                    k++;
                }
            }
        }

        { // get mappings (each map is a list of indices into the output array)
            vector<int> tmp;
            network.read_contained_vals ("maps", tmp);
            nMaps = tmp.size()/nLocations;
            Maps.resize(nMaps,vector<int>(nLocations));
            int k=0;
            for(int i=0;i<nMaps;i++){
                for(int j=0;j<nLocations;j++){
                    Maps[i][j] = tmp[k];
                    k++;
                }
            }
        }

        { // get output array (target output patterns)
            vector<double> tmp;
            network.read_contained_vals ("outPatterns", tmp);
            nPatterns = tmp.size() / outputID.size();
            Outs.resize(nPatterns,vector<double>(outputID.size()));
            int k=0;
            for(int i=0;i<nPatterns;i++){
                for(int j=0;j<outputID.size();j++){
                    Outs[i][j] = tmp[k];
                    k++;
                }
            }
        }

        { // re-arange target outputs by field identity (for uniform sampling of fields)
            vector<vector<int> > MapsUnique(nMaps);
            for(int i=0;i<nMaps;i++){
                MapsUnique[i] = getUnique(Maps[i]);
            }

            LocationsByField.resize(nMaps);
            for(int i=0;i<nMaps;i++){
                LocationsByField[i].resize(MapsUnique[i].size());
                for(int j=0;j<MapsUnique[i].size();j++){
                    for(int k=0;k<nLocations;k++){
                        if(Maps[i][k]==MapsUnique[i][j]){
                            LocationsByField[i][j].push_back(k);
                        }
                    }
                }
            }
        }

        // initialize network
        P.init (N,inputID,outputID,dt,tauW,tauX,tauY,weightNudgeSize,divergenceThreshold,maxConvergenceSteps);

        for(int i=0;i<pre.size();i++){ P.connect(pre[i],post[i]); }
        P.addBias();
        P.setNet();
        inputs.resize(inputID.size());

    }

    void setRandomPattern(){
        mapIndex = floor(morph::Tools::randDouble()*LocationsByField.size());
        int fieldIndex = floor(morph::Tools::randDouble()*LocationsByField[mapIndex].size());
        int locIndex = floor(morph::Tools::randDouble()*LocationsByField[mapIndex][fieldIndex].size());
        int locationIndex = LocationsByField[mapIndex][fieldIndex][locIndex];
        for(int i=0;i<nIns;i++){ inputs[i] = Ins[i][locationIndex]; }
        P.reset(inputs, Outs[Maps[mapIndex][locationIndex]]);
    }


    void run(int K, int errorSamplePeriod, int errorSampleSize){

        P.randomizeWeights(-1.0, +1.0);
        double errMin = 1e9;
        for(int k=0;k<K;k++){
            if(k%errorSamplePeriod){
                setRandomPattern();
                P.convergeForward(knockoutID[mapIndex-1],true);
                P.convergeBackward(knockoutID[mapIndex-1],true);
                P.weightUpdate();
            } else {
                double err = 0.;
                for(int j=0;j<errorSampleSize;j++){
                    setRandomPattern();
                    P.convergeForward(knockoutID[mapIndex-1],false);
                    err += P.getError();
                }
                err /= (double)errorSampleSize;
                if(err<errMin){
                    errMin = err;
                    P.Wbest = P.W;
                } else {
                    P.W = P.Wbest;
                }
                Error.push_back(errMin);
            }
            if(!(k%10000)){
                logfile<<"steps: "<<(int)(100*(float)k/(float)K)<<"% ("<<k<<")"<<endl;
            }
        }
        P.W = P.Wbest;
    }



    void test(void){
        // TESTING
        logfile<<"Testing..."<<endl;
        response.resize(nMaps*nLocations*P.X.size(),0.);
        int v=0;
        for(int i=0;i<nMaps;i++){
            for(int j=0;j<nLocations;j++){
                for(int k=0;k<nIns;k++){ inputs[k] = Ins[k][j];}
                P.reset(inputs, Outs[Maps[i][j]]);
                P.convergeForward(knockoutID[i-1],false);
                for(int l=0;l<P.X.size();l++){
                    response[v]=P.X[l];
                    v++;
                }
            }
        }
    }

    void disp(void){

        // Field colours
        vector<vector<double> > cols2;
        {const double tmp[] = {0.8,0.8,0.8}; cols2.push_back(makeVector(tmp));} // none
        {const double tmp[] = {1.,0.,0.}; cols2.push_back(makeVector(tmp));} // V1 col
        {const double tmp[] = {0.,0.,1.}; cols2.push_back(makeVector(tmp));} // S1 col
        {const double tmp[] = {0.,1.,0.}; cols2.push_back(makeVector(tmp));} // M1 col
        {const double tmp[] = {1.,0.,1.}; cols2.push_back(makeVector(tmp));} // A1 col
        {const double tmp[] = {0.,0.,0.}; cols2.push_back(makeVector(tmp));} // mixed

        // Displays
        vector<double> fix(3,0.0);
        vector<morph::Gdisplay> displays;
        displays.push_back(morph::Gdisplay(600, 600, 0, 0, "Image", 1.25, 0.0, 0.0));
        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();

        vector<vector<double> > resp;
        {
            stringstream fname; fname << logpath << "/outputs.h5";
            HdfData data(fname.str(),1);
            vector<double> tmp;
            data.read_contained_vals ("responses", tmp);
            int nNodes = tmp.size()/(nMaps*nLocations);
            resp.resize(nMaps*nLocations,vector<double>(nNodes,0.));
            int k=0;
            int l=0;
            for(int i=0;i<nMaps*nLocations;i++){
                for(int j=0;j<nNodes;j++){
                    resp[k][j]=tmp[l];
                    l++;
                }
                k++;
            }
        }

        double pixelwidth = fmax(fabs(X[10]-X[9]),fabs(Y[10]-Y[9])); // HACK: ASSUMES 10 and 9 are adjacent!
        double maxX, maxY = -1e9;
        double minX, minY = 1e9;
        for(int i=0;i<nLocations;i++){
            if(X[i]>maxX){ maxX=X[i];}
            if(Y[i]>maxY){ maxY=Y[i];}
            if(X[i]<minX){ minX=X[i];}
            if(Y[i]<minY){ minY=Y[i];}
        }
        double Xoff = (maxX-minX)*0.5;
        double Yoff = (maxY-minY)*0.5;

        for(int j=0;j<nMaps;j++){
            int ioff = j*nLocations;

            displays[0].resetDisplay(fix,fix,fix);
            for(int i=0;i<nLocations;i++){
                vector<double> q;
                for(int k=0;k<nPatterns;k++){
                    double sum = 0.;
                    for(int l=0;l<P.outputID.size();l++){
                        double diff = resp[ioff+i][P.outputID[l]] - Outs[k][l];
                        sum += diff*diff;
                    }
                    q.push_back(sum);
                }
                int m = getArgmin(q);
                vector<double> col = cols2[m];
                displays[0].drawRect(X[i]-Xoff,Y[i]-Yoff,0.,pixelwidth,pixelwidth,col);
            }
            stringstream ss3; ss3<< logpath << "/MaxAlign_";
            ss3 << j << ".png";
            displays[0].saveImage(ss3.str().c_str());
            displays[0].redrawDisplay();

        }

        displays[0].closeDisplay();
    }
};


int main (int argc, char **argv){
    if(argc<4){
        cout<<"Usage e.g.: ./build/sim/model 100 12 (where 100 is steps (use 0 to run tests), 12 is seed"<<endl<<flush;
        return 0;
    }

    string logpath = argv[1];
    int K=1;
    int mode = stoi(argv[2]);
    if(mode){
        K = mode;
        mode = 1;
    }
    srand(stoi(argv[3]));
    stringstream nname; nname << logpath << "/network.h5";
    Areas M(logpath,nname.str());

    switch(mode){
        case(1):{           // TRAINING

            M.run(K,1000,100);
            M.test();
            M.saveOutputs();
            M.saveWeights();

        } break;

        case(0): {        // PLOTTING
            M.disp();
        } break;

    default: {
            cout<<"Invalid mode selected"<<endl;
        } break;
    }
    return 0;
}
