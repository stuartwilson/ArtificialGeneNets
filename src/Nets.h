/*
 Implementation of recurrent backprop algorithm following Pineda (1987)
 */
#include "morph/HdfData.h"
#include "morph/display.h"
#include "morph/tools.h"
#include <morph/Config.h>
#include "tools.h"
#include "Pineda.h"

using morph::Config;
using namespace morph;
using namespace tools;

class Map{

    /*
      Structure for storing a map (pre-defined X, Y, Z, F, in a HdfData file)
    */

    public:
        int N;
        vector<double> Xscaled, Yscaled, Zscaled, Fscaled;
        vector<double> X, Y, Z, F;
        double minX, maxX, minY, maxY, minZ, maxZ, minF, maxF, xScale, yScale, xSep, ySep;
        int outputID, contextID;
        double contextVal;

    void init(string filename){
        HdfData network(filename,1);
        network.read_contained_vals ("X", X);
        network.read_contained_vals ("Y", Y);
        network.read_contained_vals ("Z", Z);
        network.read_contained_vals ("F", F);

        if((!(X.size()==Y.size())) || (!(X.size()==Z.size())) || (!(X.size()==F.size()))){
            cout<<"X, Y, Z, F not all the same size... all kinds of badness!"<<endl;
        }
        N = X.size();

        minX = tools::getMin(X);
        minY = tools::getMin(Y);
        minZ = tools::getMin(Z);
        minF = tools::getMin(F);
        maxX = tools::getMax(X);
        maxY = tools::getMax(Y);
        maxZ = tools::getMax(Z);
        maxF = tools::getMax(F);
        Xscaled = tools::getRenormedVector(X);
        Yscaled = tools::getRenormedVector(Y);
        Zscaled = tools::getRenormedVector(Z);
        Fscaled = tools::getRenormedVector(F);

        double maxDim = maxX-minX;
        if((maxY-minY)>maxDim){ maxDim = maxY-minY; };
        xScale = (maxX-minX)/maxDim;
        yScale = (maxY-minY)/maxDim;

        // deduce number of rows and columns
        vector<double> uniqueX = tools::getUnique(X);
        vector<double> uniqueY = tools::getUnique(Y);
        xSep = xScale/(uniqueX.size()-1);
        ySep = yScale/(uniqueY.size()-1);
    }

    Map(string filename){
        init(filename);
    }

    Map(string filename,int outputID){
        init(filename);
        this->outputID = outputID;
    }

    Map(string filename,int outputID, int contextID, double contextVal){
        init(filename);
        this->outputID = outputID;
        this->contextID = contextID;
        this->contextVal = contextVal;
    }
};

class Net{

public:
    string logpath;
    ofstream logfile;
    vector<double> inputs, Error, response;
    Pineda P;
    vector<Map> M;
    int mapID, locID;
    vector<int> inputID;

    Net(string logpath){

        // setup log file
        this->logpath = logpath;
        morph::Tools::createDir (logpath);
        { stringstream ss; ss << logpath << "/log.txt"; logfile.open(ss.str());}
        logfile<<"Hello."<<endl;

        // Read in network params
        Config conf;
        { stringstream ss; ss << logpath <<"/config.json"; conf.init (ss.str()); }
        float dt = conf.getFloat("dt",1.0);
        float tauW = conf.getFloat("tauW",2.0);
        float tauX = conf.getFloat("tauX",1.0);
        float tauY = conf.getFloat("tauY",1.0);
        float weightNudgeSize = conf.getFloat("weightNudgeSize",0.001);
        float divergenceThreshold = conf.getFloat("divergenceThreshold",0.000001);
        int maxConvergenceSteps = conf.getInt("maxConvergenceSteps",400);

        // Read in map info
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

        // Setup network connectivity
        vector<int> pre, post;
        stringstream ss; ss << logpath << "/network.h5"; HdfData network(ss.str(),1);
        network.read_contained_vals ("pre", pre);
        network.read_contained_vals ("post", post);

        if(pre.size()!=post.size()){ logfile<<"Pre/Post different sizes ("<<pre.size()<<"/"<<post.size()<<")"<<endl; exit(0);}
        if(pre.size()<1){ logfile<<"No connections in network!"<<endl; exit(0);}

        int N = tools::getMax(pre);
        if(tools::getMax(post)>N){
            N=tools::getMax(post);
        }
        N++;

        // Initiate network
        P.init (N,dt,tauW,tauX,tauY,weightNudgeSize,divergenceThreshold,maxConvergenceSteps);
        for(int i=0;i<pre.size();i++){ P.connect(pre[i],post[i]); }
        P.addBias();
        P.setNet();
        inputs.resize(inputID.size());

    }

    void saveOutputs(void) { // log outputs
        stringstream fname; fname << logpath << "/outputs.h5";
        HdfData outdata(fname.str());
        outdata.add_contained_vals ("error", Error);
        outdata.add_contained_vals ("responses", response);
    }

    void saveWeights(void){ // log weights
        stringstream fname; fname << logpath << "/weights.h5";
        HdfData weightdata(fname.str());
        weightdata.add_contained_vals ("weights", P.W);
        vector<double> flatweightmat = P.getWeightMatrix();
        weightdata.add_contained_vals ("weightmat", flatweightmat);
    }

    void loadWeights(void) { // log outputs
        stringstream fname; fname << logpath << "/weights.h5";
        HdfData loaded(fname.str(),1);
        loaded.read_contained_vals ("weights", P.W);
        P.Wbest = P.W;
    }

    ~Net(void){
        logfile<<"Goodbye."<<endl;
        logfile.close();
    }

    void setMap(int i){ mapID = i; }

    void sampleMap(int j){ locID = j; }

    void plotMaps(void){

        // Displays
        vector<double> fix(3,0.0);
        vector<Gdisplay> displays;
        displays.push_back(Gdisplay(600, 600, 0, 0, "Map", 1.3, 0.0, 0.0));
        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();

        // plot supplied map values
        for(int j=0;j<M.size();j++){

            displays[0].resetDisplay(fix,fix,fix);
            displays[0].resetDisplay(fix,fix,fix);

            for(int i=0;i<M[j].N;i++){
                //vector<double> rgb = morph::Tools::getJetColor(M[j].Fscaled[i]);
                vector<double> rgb = morph::Tools::getGrayScaleColor(M[j].Fscaled[i]);
                displays[0].drawRect(M[j].xScale*(M[j].Xscaled[i]-0.5),M[j].yScale*(M[j].Yscaled[i]-0.5),0.,M[j].xSep,M[j].ySep,rgb);
            }
            stringstream ss; ss<< logpath << "/m_" << j << ".png";
            displays[0].saveImage(ss.str().c_str());
            displays[0].redrawDisplay();
        }
    }

    vector<vector<double> > test(int i){
        // return response of every node at each location in map i

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
        displays.push_back(Gdisplay(600, 600, 0, 0, "Response", 1.3, 0.0, 0.0));
        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();

        for(int j=0;j<M.size();j++){

            vector<vector<double> > R = test(j);
            vector<double> F = R[M[j].outputID];
            F = getRenormedVector(F);

            displays[0].resetDisplay(fix,fix,fix);
            displays[0].resetDisplay(fix,fix,fix);

            for(int i=0;i<M[j].N;i++){
                //vector<double> rgb = morph::Tools::getJetColor(F[i]);
                vector<double> rgb = morph::Tools::getGrayScaleColor(F[i]);
                displays[0].drawRect(M[j].xScale*(M[j].Xscaled[i]-0.5),M[j].yScale*(M[j].Yscaled[i]-0.5),0.,M[j].xSep,M[j].ySep,rgb);
            }
            stringstream ss; ss<< logpath << "/x_" << j << ".png";
            displays[0].saveImage(ss.str().c_str());
            displays[0].redrawDisplay();
        }

        displays[0].closeDisplay();
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
                    P.W = P.Wbest;
                }
                Error.push_back(errMin);
            }
            if(!(k%(K/100))){
                logfile<<"steps: "<<(int)(100*(float)k/(float)K)<<"% ("<<k<<")"<<endl;
            }
        }
        P.W = P.Wbest;
    }

};



