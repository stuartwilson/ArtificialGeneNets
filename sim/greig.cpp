/*
 Implementation of recurrent backprop algorithm by Pineda (1987)
 */

#include "morph/HdfData.h"
#include "morph/display.h"
#include "morph/tools.h"
#include "Pineda.h"
#include <iostream>

using namespace std;
using namespace morph;

template<typename T, size_t N> std::vector<T> makeVector(const T (&data)[N]){return std::vector<T>(data,data+N);}

int getArgmax(vector<double> q){
    double maxV = -1e9;
    int maxI = 0;
    for(int i=0;i<q.size();i++){
        if(q[i]>maxV){
            maxV = q[i];
            maxI = i;
        }
    }
    return maxI;
}

int getArgmin(vector<double> q){
    double minV = 1e9;
    int minI = 0;
    for(int i=0;i<q.size();i++){
        if(q[i]<minV){
            minV = q[i];
            minI = i;
        }
    }
    return minI;
}


int main (int argc, char **argv){

    if(argc<5){
        cout<<"Usage e.g.: ./build/sim/model configs/config.json logs 100 12 (where 0 is steps (use 0 to run tests), 12 is seed"<<endl<<flush;
        return 0;
    }

    string paramsfile (argv[1]);
    string logpath = argv[2];
    morph::Tools::createDir (logpath);
    int K=1;
    int mode = stoi(argv[3]);
    if(mode){
        K = mode;
        mode = 1;
    }
    srand(stoi(argv[4]));


    ofstream logfile;
    {
        stringstream ss; ss << logpath << "/log.txt";
        logfile.open(ss.str());
    }
    logfile << "Hello! "<<endl<<"Running for "<<K<<" iterations."<<endl;

    // JSON stuff
    ifstream jsonfile_test;
    int srtn = system ("pwd");
    if (srtn) { cerr << "system call returned " << srtn << endl;}
    jsonfile_test.open (paramsfile, ios::in);
    if (jsonfile_test.is_open()) { jsonfile_test.close();}
    else { cerr << "json config file " << paramsfile << " not found." << endl; return 1;}
    ifstream jsonfile (paramsfile, ifstream::binary);
    Json::Value root;
    string errs;
    Json::CharReaderBuilder rbuilder;
    rbuilder["collectComments"] = false;
    bool parsingSuccessful = Json::parseFromStream (rbuilder, jsonfile, &root, &errs);
    if (!parsingSuccessful) { cerr << "Failed to parse JSON: " << errs; return 1;}

    // Get Params

    const float dt = root.get("dt",1.0).asFloat();
    const float tauW = root.get("tauW",2.0).asFloat();
    const float tauX = root.get("tauX",1.0).asFloat();
    const float tauY = root.get("tauY",1.0).asFloat();
    const int errorSamplePeriod = root.get("errorSamplePeriod",1000).asInt();
    const int errorSampleSize = root.get("errorSampleSize",100).asInt();

    const float weightNudgeSize = root.get("weightNudgeSize",0.001).asFloat();
    const float divergenceThreshold = root.get("divergenceThreshold",0.000001).asFloat();
    const int maxConvergenceSteps = root.get("maxConvergenceSteps",400).asInt();

    stringstream iname; iname << logpath << "/inputs.h5";
    HdfData input(iname.str(),1);
    vector<double> X, Y;
    input.read_contained_vals ("x", X);
    input.read_contained_vals ("y", Y);

    int nLocations = X.size();

    int nIns;
    vector<vector<double> > Ins;
    {
        vector<double> tmp;
        input.read_contained_vals ("inPatterns", tmp);
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


    int nMaps;
    vector<vector<int> > Maps;
    {
        vector<int> tmp;
        input.read_contained_vals ("maps", tmp);
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


    stringstream nname; nname << logpath << "/network.h5";
    HdfData network(nname.str(),1);
    vector<int> Ntmp(0);
    network.read_contained_vals ("N", Ntmp);
    int N = Ntmp[0];
    vector<int> inputID, outputID, knockoutID, pre, post;
    network.read_contained_vals ("inputs", inputID);
    network.read_contained_vals ("outputs", outputID);
    network.read_contained_vals ("knockouts", knockoutID);
    network.read_contained_vals ("pre", pre);
    network.read_contained_vals ("post", post);


    vector<vector<double> > Outs;

    int nPatterns;
    {
        vector<double> tmp;
        input.read_contained_vals ("outPatterns", tmp);
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


        Pineda P (N,inputID,outputID,dt,tauW,tauX,tauY,weightNudgeSize, divergenceThreshold,maxConvergenceSteps);



    for(int i=0;i<pre.size();i++){ P.connect(pre[i],post[i]); }

    P.addBias();
    P.setNet();
    vector<double> inputs(inputID.size());


    switch(mode){

    case(1):{           // TRAINING

        P.randomizeWeights(-1.0, +1.0);

        vector<double> Error;
        double errMin = 1e9;

        for(int k=0;k<K;k++){

            if(k%errorSamplePeriod){

                int mapIndex = floor(morph::Tools::randDouble()*nMaps);
                int locationIndex = floor(morph::Tools::randDouble()*nLocations);
                for(int i=0;i<nIns;i++){ inputs[i] = Ins[i][locationIndex]; }
                P.reset(inputs, Outs[Maps[mapIndex][locationIndex]]);
                P.convergeForward(knockoutID[mapIndex-1],true);
                P.convergeBackward(knockoutID[mapIndex-1],true);
                P.weightUpdate();

            } else {
                double err = 0.;
                for(int j=0;j<errorSampleSize;j++){
                    int mapIndex = floor(morph::Tools::randDouble()*nMaps);
                    int locationIndex = floor(morph::Tools::randDouble()*nLocations);
                    for(int i=0;i<nIns;i++){ inputs[i] = Ins[i][locationIndex]; }
                    P.reset(inputs, Outs[Maps[mapIndex][locationIndex]]);
                    P.convergeForward(knockoutID[mapIndex-1],false);
                    err += P.getError();
                }
                err /= (double)errorSampleSize;
                if(err<errMin){
                    errMin = err;
                    //P.Wbest = P.W;
                }
                Error.push_back(err);
            }

            if(!(k%10000)){
                logfile<<"steps: "<<(int)(100*(float)k/(float)K)<<"% ("<<k<<")"<<endl;
            }

        }

        //P.W = P.Wbest;

        // TESTING
        logfile<<"Testing..."<<endl;
        vector<double> response;
        for(int i=0;i<nMaps;i++){
            for(int j=0;j<nLocations;j++){
                for(int k=0;k<nIns;k++){ inputs[k] = Ins[k][j];}
                P.reset(inputs, Outs[Maps[i][j]]);
                P.convergeForward(knockoutID[i-1],false);
                for(int l=0;l<P.X.size();l++){
                    response.push_back(P.X[l]);
                }
            }
        }

        { // log outputs
            stringstream fname; fname << logpath << "/outputs.h5";
            HdfData outdata(fname.str());
            outdata.add_contained_vals ("error", Error);
            outdata.add_contained_vals ("responses", response);
        }

        { // log weights
            stringstream fname; fname << logpath << "/weights.h5";
            HdfData weightdata(fname.str());
            weightdata.add_contained_vals ("weights", P.W);
            vector<double> flatweightmat = P.getWeightMatrix();
            weightdata.add_contained_vals ("weightmat", flatweightmat);
        }

    } break;





    case(0): {        // PLOTTING

        // Displays
        vector<double> fix(3,0.0);
        vector<morph::Gdisplay> displays;
        displays.push_back(morph::Gdisplay(600, 600, 0, 0, "Image", 1.7, 0.0, 0.0));
        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();

        {
            // Loading
            stringstream fname; fname << logpath << "/weights.h5";
            HdfData data(fname.str(),1);
            data.read_contained_vals ("weights", P.W);
        }

        vector<vector<double> > response;
        {
            stringstream fname; fname << logpath << "/outputs.h5";
            HdfData data(fname.str(),1);
            vector<double> tmp;
            data.read_contained_vals ("responses", tmp);
            int nNodes = tmp.size()/(nMaps*nLocations);
            response.resize(nMaps*nLocations,vector<double>(nNodes,0.));
            int k=0;
            int l=0;
            for(int i=0;i<nMaps*nLocations;i++){
                for(int j=0;j<nNodes;j++){
                    response[k][j]=tmp[l];
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

            // max field id
            vector<vector<double> > cols2;
            {const double tmp[] = {0.8,0.8,0.8}; cols2.push_back(makeVector(tmp));} // none
            {const double tmp[] = {1.,0.,0.}; cols2.push_back(makeVector(tmp));} // V1 col
            {const double tmp[] = {0.,0.,1.}; cols2.push_back(makeVector(tmp));} // S1 col
            {const double tmp[] = {0.,1.,0.}; cols2.push_back(makeVector(tmp));} // M1 col
            {const double tmp[] = {1.,0.,1.}; cols2.push_back(makeVector(tmp));} // A1 col
            {const double tmp[] = {0.,0.,0.}; cols2.push_back(makeVector(tmp));} // mixed
            displays[0].resetDisplay(fix,fix,fix);
            for(int i=0;i<nLocations;i++){
                vector<double> q;
                for(int j=0;j<nPatterns;j++){
                    double sum = 0.;
                    for(int k=0;k<outputID.size();k++){
                        double diff = response[ioff+i][outputID[k]] - Outs[j][k];
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

    } break;

    default: {
            cout<<"Invalid mode selected"<<endl;
        } break;
    }

    logfile<<"Goodbye."<<endl;
    logfile.close();
    return 0;
}
