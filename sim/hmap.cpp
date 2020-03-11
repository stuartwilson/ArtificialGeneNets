#include "Nets.h"
#include <iostream>
using namespace std;

#include <morph/Config.h>
using morph::Config;

class Context: public Net{

public:

    vector<int> inputID;

    Context(string logpath, string nname) : Net(logpath){

        Config conf;
        stringstream ss; ss << logpath <<"/config.json";
        conf.init (ss.str());
        float dt = conf.getFloat("dt",1.0);
        float tauW = conf.getFloat("tauW",2.0);
        float tauX = conf.getFloat("tauX",1.0);
        float tauY = conf.getFloat("tauY",1.0);
        float weightNudgeSize = conf.getFloat("weightNudgeSize",0.001);
        float divergenceThreshold = conf.getFloat("divergenceThreshold",0.000001);
        int maxConvergenceSteps = conf.getInt("maxConvergenceSteps",400);

        const Json::Value maps = conf.getArray("maps");
        for(int i=0;i<maps.size();i++){
            string fn = maps[i].get("filename", "unknown map").asString();
            stringstream ss; ss << logpath <<"/"<<fn;
            logfile<<"Map["<<i<<"]:"<<ss.str()<<endl;
            int oID = maps[i].get("outputID",-1).asInt();
            int cID = maps[i].get("contextID",-1).asInt();
            float cVal = maps[i].get("contextVal",-1.0).asFloat();
            M.push_back(Map(ss.str(),oID,cID,cVal));
        }

        const Json::Value inp = conf.getArray("inputID");
        for(int i=0;i<inp.size();i++){
            inputID.push_back(inp[i].asInt());
        }

        vector<int> pre, post;
        HdfData network(nname,1);
        network.read_contained_vals ("pre", pre);
        network.read_contained_vals ("post", post);

        if(pre.size()!=post.size()){ logfile<<"Pre/Post different sizes ("<<pre.size()<<"/"<<post.size()<<")"<<endl; exit(0);}
        if(pre.size()<1){ logfile<<"No connections in network!"<<endl; exit(0);}

        int N = tools::getMax(pre);
        if(tools::getMax(post)>N){
            N=tools::getMax(post);
        }
        N++;

        P.init (N,dt,tauW,tauX,tauY,weightNudgeSize,divergenceThreshold,maxConvergenceSteps);

        for(int i=0;i<pre.size();i++){ P.connect(pre[i],post[i]); }

        P.addBias();
        P.setNet();
        inputs.resize(inputID.size());

    }

    void setInput(void){
        P.reset();
        P.Input[inputID[0]] = M[mapID].X[locID];
        P.Input[inputID[1]] = M[mapID].Y[locID];
        P.Input[M[mapID].contextID] = M[mapID].contextVal;
    }

    void setRandomInput(void){
        setMap(floor(morph::Tools::randDouble()*M.size()));
        sampleMap(floor(morph::Tools::randDouble()*M[mapID].N));
        setInput();
    }

    void run(int K, int errorSamplePeriod, int errorSampleSize, bool resetWeights){

        P.randomizeWeights(-1.0, +1.0);

        double errMin = 1e9;
        for(int k=0;k<K;k++){
            if(k%errorSamplePeriod){
                setRandomInput();
                P.convergeForward(-1,true);
                P.setError(vector<int> (1,M[mapID].outputID), vector<double> (1,M[mapID].F[locID]));
                P.convergeBackward(-1,false);
                P.weightUpdate();
            } else {
                double err = 0.;
                for(int j=0;j<errorSampleSize;j++){
                    setRandomInput();
                    P.convergeForward(-1,false);
                    P.setError(vector<int> (1,M[mapID].outputID), vector<double> (1,M[mapID].F[locID]));
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
        setMap(i);
        vector<vector<double> > response(P.N,vector<double>(M[mapID].N,0.));
        for(int j=0;j<M[mapID].N;j++){
            sampleMap(j);
            setInput();
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
            vector<double> F = R[M[j].outputID];

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
        cout<<"not enough command line arguments"<<endl<<flush;
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
