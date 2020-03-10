#include "Nets.h"
#include <iostream>
using namespace std;


class Context: public Net{

public:

    Context(string logpath, string nname) : Net(logpath){

        // SHOULD BE ABLE TO DEFINE THESE EXTERNALLY
        float dt = 1.0;
        float tauW = 2.0;
        float tauX = 1.0;
        float tauY = 1.0;
        float weightNudgeSize = 0.001;
        float divergenceThreshold = 0.000001;
        int maxConvergenceSteps = 400;

        HdfData network(nname,1);

        M.push_back(Map("/Users/stuartwilson/MODELS/gitprojects/ArtificialGeneNets/model/maps/sighted.h5"));
        M.push_back(Map("/Users/stuartwilson/MODELS/gitprojects/ArtificialGeneNets/model/maps/enucleate.h5"));

        vector<int> Ntmp(0);
        network.read_contained_vals ("N", Ntmp);
        int N = Ntmp[0];
        vector<int> inputID, outputID, knockoutID, contextID, pre, post;
        network.read_contained_vals ("inputs", inputID);
        network.read_contained_vals ("outputs", outputID);
        network.read_contained_vals ("context", contextID);
        network.read_contained_vals ("pre", pre);
        network.read_contained_vals ("post", post);

        // Additional input to represent which map we are using
        inputID.push_back(contextID[0]);

        P.init (N,inputID,outputID,dt,tauW,tauX,tauY,weightNudgeSize,divergenceThreshold,maxConvergenceSteps);

        for(int i=0;i<pre.size();i++){ P.connect(pre[i],post[i]); }

        P.addBias();
        P.setNet();
        inputs.resize(inputID.size());

    }

    void setRandomPattern(void){
        int mapIndex = floor(morph::Tools::randDouble()*M.size());
        int locationIndex = floor(morph::Tools::randDouble()*M[mapIndex].N);
        inputs[0] = M[mapIndex].X[locationIndex];
        inputs[1] = M[mapIndex].Y[locationIndex];
        inputs[2]=(double)mapIndex;
        P.reset(inputs, vector<double>(1,M[mapIndex].F[locationIndex]));
    }

    void run(int K, int errorSamplePeriod, int errorSampleSize, bool resetWeights){

        P.randomizeWeights(-1.0, +1.0);

        double errMin = 1e9;
        for(int k=0;k<K;k++){
            if(k%errorSamplePeriod){
                setRandomPattern();
                P.convergeForward(-1,true);
                P.convergeBackward(-1,false);
                P.weightUpdate();
            } else {
                double err = 0.;
                for(int j=0;j<errorSampleSize;j++){
                    setRandomPattern();
                    P.convergeForward(-1,false);
                    err += P.getError();
                }
                err /= (double)errorSampleSize;
                if(err<errMin){
                    errMin = err;
                    P.Wbest = P.W;
                } else if (resetWeights) {
                    // P.W = P.Wbest; // doesn't seem to help with Val?
                }
                Error.push_back(errMin);
            }
            if(!(k%10000)){
                logfile<<"steps: "<<(int)(100*(float)k/(float)K)<<"% ("<<k<<")"<<endl;
            }
        }
        P.W = P.Wbest;
    }

    vector<vector<double> > test(int i){
        vector<vector<double> > response(P.N,vector<double>(M[i].N,0.));
        for(int j=0;j<M[i].N;j++){
            inputs[0] = M[i].X[j];
            inputs[1] = M[i].Y[j];
            inputs[2]=(double)i;
            P.reset(inputs, vector<double>(1,M[i].F[j]));
            P.convergeForward(-1,false);
            for(int k=0;k<P.N;k++){
                response[k][j] = P.X[k];
            }
        }
        return response;
    }

    void plotResponses(void){

        // Displays
        vector<double> fix(3,0.0);
        vector<Gdisplay> displays;
        displays.push_back(Gdisplay(600, 600, 0, 0, "Image", 1.3, 0.0, 0.0));
        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();

        for(int j=0;j<M.size();j++){

            vector<vector<double> > R = test(j);
            vector<double> F = R[P.outputID[0]];

            double minF = +1e9;
            double maxF = -1e9;
            for(int i=0;i<M[j].N;i++){
                if(F[i]<minF){minF = F[i];};
                if(F[i]>maxF){maxF = F[i];};
            }
            double norm = 1./(maxF-minF);
            for(int i=0;i<M[j].N;i++){
                F[i] = (F[i]-minF)*norm;
            }

            displays[0].resetDisplay(fix,fix,fix);
            displays[0].resetDisplay(fix,fix,fix);

            for(int i=0;i<M[j].N;i++){
                vector<double> rgb = morph::Tools::getJetColor(F[i]);
                displays[0].drawRect(M[j].Xscaled[i]-0.5,M[j].Yscaled[i]-0.5,0.,1./(M[j].ncols-1),1./(M[j].nrows-1),rgb);
            }
            stringstream ss; ss<< logpath << "/fit_";
            ss << j << ".png";
            displays[0].saveImage(ss.str().c_str());
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
    Context Net(logpath,nname.str());

    switch(mode){
        case(1):{           // TRAINING

            Net.run(K,1000,100,false);
            Net.saveOutputs();
            Net.saveWeights();

        } break;

        case(0): {        // PLOTTING
            Net.loadWeights();
            Net.plotMaps();
            Net.plotResponses();
        } break;

    default: {
            cout<<"Invalid mode selected"<<endl;
        } break;
    }
    return 0;
}
