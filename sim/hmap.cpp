/*
 Implementation of recurrent backprop algorithm following Pineda (1987)
 */
#include "morph/HdfData.h"
#include "morph/display.h"
#include "morph/tools.h"
#include "Pineda.h"
#include <iostream>

using namespace std;
using namespace morph;

vector<double> getUnique(vector<double> x){
    vector<double> unique;
    for(int i=0;i<x.size();i++){
        bool uni = true;
        for(int k=0;k<unique.size();k++){
            if(x[i]==unique[k]){ uni = false; break; }
        } if(uni){ unique.push_back(x[i]);}
    }
    return unique;
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

    vector<int> dims;
    input.read_contained_vals ("dims", dims);

    int nMaps = 2; // map is a condition, e.g., sighted/enucleated
    int nLocations = dims[1];
    int nIns = nMaps*2; // (x and y for each map)

    vector<vector<double> > Ins;
    {
        vector<double> tmp;
        input.read_contained_vals ("inPatterns", tmp);
        Ins.resize(nIns,vector<double>(nLocations));
        int k=0;
        for(int i=0;i<nIns;i++){
            for(int j=0;j<nLocations;j++){
                Ins[i][j] = tmp[k];
                k++;
            }
        }
    }

    vector<vector<double> > Maps;
    {
        vector<double> tmp;
        input.read_contained_vals ("maps", tmp);
        Maps.resize(nMaps,vector<double>(nLocations));
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
    vector<int> inputID, outputID, knockoutID, contextID, pre, post;
    network.read_contained_vals ("inputs", inputID);
    network.read_contained_vals ("outputs", outputID);
    network.read_contained_vals ("context", contextID);
    network.read_contained_vals ("pre", pre);
    network.read_contained_vals ("post", post);

    // Additional input to represent which map we are using
    inputID.push_back(contextID[0]);

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
                for(int i=0;i<2;i++){
                    inputs[i] = Ins[mapIndex*2+i][locationIndex];
                }
                inputs[2]=(double)mapIndex;

                P.reset(inputs, vector<double>(1,Maps[mapIndex][locationIndex]));
                P.convergeForward(-1,true);
                P.convergeBackward(-1,false);
                P.weightUpdate();

            } else {

                double err = 0.;
                for(int j=0;j<errorSampleSize;j++){
                    int mapIndex = floor(morph::Tools::randDouble()*nMaps);
                    int locationIndex = floor(morph::Tools::randDouble()*nLocations);
                    for(int i=0;i<2;i++){
                        inputs[i] = Ins[mapIndex*2+i][locationIndex];
                    }
                    inputs[2]=(double)mapIndex;

                    P.reset(inputs, vector<double>(1,Maps[mapIndex][locationIndex]));
                    P.convergeForward(-1,false);
                    err += P.getError();
                }
                err /= (double)errorSampleSize;
                if(err<errMin){
                    errMin = err;
                    P.Wbest = P.W;
                }
                Error.push_back(err);
            }

            if(!(k%10000)){
                logfile<<"steps: "<<(int)(100*(float)k/(float)K)<<"% ("<<k<<")"<<endl;
            }

        }

        P.W = P.Wbest;

        // TESTING
        logfile<<"Testing..."<<endl;
        vector<double> response;
        for(int i=0;i<nMaps;i++){
            for(int j=0;j<nLocations;j++){
                for(int k=0;k<2;k++){
                    inputs[k] = Ins[i*2+k][j];
                }
                inputs[2]=(double)i;
                P.reset(inputs, vector<double>(1,Maps[i][j]));
                P.convergeForward(-1,false);
                response.push_back(P.X[outputID[0]]);
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
        vector<Gdisplay> displays;
        displays.push_back(Gdisplay(600, 600, 0, 0, "Image", 1.3, 0.0, 0.0));
        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();

        // deduce number of rows and columns
        vector<vector<int> > nunique(nMaps,vector<int>(2,0));
        for (int j=0;j<nMaps;j++){
            vector<double> uniqueX = getUnique(Ins[j*2]);
            vector<double> uniqueY = getUnique(Ins[j*2+1]);
            nunique[j][0] = uniqueX.size();
            nunique[j][1] = uniqueY.size();
        }

        // joint normalize x, y extrema for plotting
        double maxX, maxY = -1e9;
        double minX, minY = 1e9;
        {
            double x,y;
            for (int j=0;j<nMaps;j++){
                for(int i=0;i<nLocations;i++){
                    x=Ins[j*2][i];
                    y=Ins[j*2+1][i];
                    if(x>maxX){ maxX=x;}
                    if(x<minX){ minX=x;}
                    if(y>maxY){ maxY=y;}
                    if(y<minY){ minY=y;}
                }
            }
        }

        // individually normalize z values (color) for each map
        vector<vector<double> > normedMaps = Maps;
        for(int i=0;i<Maps.size();i++){
            double maxM = -1e9;
            double minM =  1e9;
            for(int j=0;j<Maps[i].size();j++){
                if(maxM<Maps[i][j]){ maxM=Maps[i][j];}
                if(minM>Maps[i][j]){ minM=Maps[i][j];}
            }
            double normM = 1./(maxM-minM);
            for(int j=0;j<Maps[i].size();j++){
                normedMaps[i][j] = (Maps[i][j]-minM)*normM;
            }
        }

        // load pre-computed responses
        vector<double > response;
        {
            stringstream fname; fname << logpath << "/outputs.h5";
            HdfData data(fname.str(),1);
            data.read_contained_vals ("responses", response);
        }

        // individually normalize z values (color) for each fit
        vector<vector<double> > normedFits = Maps;
        for(int j=0;j<Maps.size();j++){
            int ioff = j*nLocations;
            double maxZ=-1e9;
            double minZ=+1e9;
            for(int i=0;i<nLocations;i++){
                if(response[ioff+i]>maxZ){maxZ=response[ioff+i];}
                if(response[ioff+i]<minZ){minZ=response[ioff+i];}
            }
            double normZ = 1./(maxZ-minZ);
            if(maxZ==minZ){ normZ = 0.; }
            for(int i=0;i<nLocations;i++){
                normedFits[j][i] = (response[ioff+i]-minZ)*normZ;
            }
        }

        // plot supplied map values
        for(int j=0;j<nMaps;j++){

            displays[0].resetDisplay(fix,fix,fix);
            displays[0].resetDisplay(fix,fix,fix);

            for(int i=0;i<nLocations;i++){
                double x = Ins[j*2][i]-0.5;
                double y = Ins[j*2+1][i]-0.5;
                vector<double> rgb = morph::Tools::getJetColor(normedMaps[j][i]);
                displays[0].drawRect(x,y,0.,1./(double)nunique[j][0],1./(double)nunique[j][1],rgb);
            }
            stringstream ss2; ss2<< logpath << "/map_";
            ss2 << j << ".png";
            displays[0].saveImage(ss2.str().c_str());
            displays[0].redrawDisplay();
        }


        // plot responses
        for(int j=0;j<nMaps;j++){

            displays[0].resetDisplay(fix,fix,fix);
            displays[0].resetDisplay(fix,fix,fix);

            for(int i=0;i<nLocations;i++){
                double x = Ins[j*2][i]-0.5;
                double y = Ins[j*2+1][i]-0.5;
                vector<double> rgb = morph::Tools::getJetColor(normedFits[j][i]);
                displays[0].drawRect(x,y,0.,1./nunique[j][0],1./nunique[j][1],rgb);
            }
            stringstream ss2; ss2<< logpath << "/fit_";
            ss2 << j << ".png";
            displays[0].saveImage(ss2.str().c_str());
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
