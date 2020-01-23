/*
 Implementation of (continous) recurrent backprop algorithm by Pineda (1987).
 */

#include "morph/display.h"
#include "morph/tools.h"
#include "morph/HdfData.h"
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

    if(argc<3){
        cout<<"Usage e.g.: ./build/sim/model configs/config.json logs 0 12 (where 0 is mode and 12 is optional seed)"<<endl<<flush;
        return 0;
    }

    string paramsfile (argv[1]);
    string logpath = argv[2];
    morph::Tools::createDir (logpath);
    int mode = stoi(argv[3]);
    if(argc>4){ srand(stoi(argv[4])); }

    stringstream ss; ss << logpath << "/log.txt";
    ofstream logfile;
    logfile.open(ss.str());
    logfile << "Hello! "<<endl;

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

    //const int mode = root.get ("mode", -1).asInt();
    const unsigned int T = root.get ("T", 400).asUInt();
    const unsigned int K = root.get ("K", 50000).asUInt();
    const float taux = root.get ("taux", 1).asFloat();
    const float tauz = root.get ("tauz", 1).asFloat();
    const float tauw = root.get ("tauw", 32).asFloat();
    const float dt = root.get ("dt", 0.02).asFloat();

    string inputfilepath = root.get("inputfilepath", "no input file selected").asString();
    logfile<<inputfilepath<<endl;

    string networkfilepath = root.get("networkfilepath", "no network file selected").asString();
    logfile<<networkfilepath<<endl;

    string weightfilepath = root.get("weightfilepath", "no weight file selected").asString();
    logfile<<weightfilepath<<endl;

    stringstream iname; iname << inputfilepath;
    HdfData input(iname.str(),1);
    vector<double> X, Y;
    input.read_contained_vals ("x", X);
    input.read_contained_vals ("y", Y);

    int nLocations = X.size();

    int nIns;
    vector<vector<int> > Ins;
    {
        vector<int> tmp;
        input.read_contained_vals ("inPatterns", tmp);
        nIns = tmp.size()/nLocations;
        Ins.resize(nIns,vector<int>(nLocations));
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

    stringstream nname; nname << networkfilepath;
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
    vector<double> weightbounds;
    network.read_contained_vals ("weightbounds", weightbounds);

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

    Pineda P (N,inputID,outputID, taux, tauz, tauw, dt);

    for(int i=0;i<pre.size();i++){ P.connect(pre[i],post[i]); }

    P.addBias();

    P.setNet();

    vector<double> inputs(inputID.size());

    switch(mode){

    case(0):{           // TRAINING

        P.randomizeWeights(weightbounds[0], weightbounds[1]);

        vector<double> Error;

        for(int k=0;k<K;k++){

            if(!(k%1000)){logfile<<"steps: "<<k<<endl;}

            int mapIndex = floor(morph::Tools::randDouble()*nMaps);
            int locationIndex = floor(morph::Tools::randDouble()*nLocations);

            inputs[0] = X[locationIndex];
            inputs[1] = Y[locationIndex];
            P.reset(inputs, Outs[Maps[mapIndex][locationIndex]]);

            for(int t=0;t<T;t++){
                P.forward();
                if(mapIndex>0){
                    P.X[knockoutID[mapIndex-1]] = 0.;
                }
                P.backward();
                if(mapIndex>0){
                    P.Z[knockoutID[mapIndex-1]] = 0.;
                }
                P.weightUpdate();
            }
            Error.push_back(P.getError());
        }

        {   // log outputs
            stringstream fname; fname << logpath << "/out.h5";
            HdfData data(fname.str());
            stringstream ss; ss<<"error";
            data.add_contained_vals (ss.str().c_str(), Error);;
        }

        { // log weights
            stringstream fname; fname << logpath << "/net.h5";
            HdfData data(fname.str());
            stringstream ss; ss<<"weights";
            data.add_contained_vals (ss.str().c_str(), P.W);
        }
    break;
    }





    case(1): {        // TESTING

        // Displays
        vector<double> fix(3, 0.0);
        vector<morph::Gdisplay> displays;
        displays.push_back(morph::Gdisplay(600, 600, 0, 0, "Image", 1.7, 0.0, 0.0));

        stringstream nname; nname << weightfilepath;
        HdfData network(nname.str(),1);
        vector<double> tmp;
        network.read_contained_vals ("weights", tmp);
        for(int i=0;i<tmp.size();i++){
            P.W[i] = tmp[i];
        }

        // TESTING
        vector<vector<double> > response;
        for(int i=0;i<nMaps;i++){
            for(int j=0;j<nLocations;j++){
                inputs[0] = X[j];
                inputs[1] = Y[j];
                P.reset(inputs, Outs[Maps[i][j]]);
                for(int t=0;t<T;t++){
                    P.forward();
                    if(i>0){
                        P.X[knockoutID[i-1]] = 0.;
                    }
                }
                response.push_back(P.X);
            }
        }


        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();
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

        vector<vector<double> > cols;
        {const double tmp[] = {1.,0.,0.}; cols.push_back(makeVector(tmp));} // V1 col
        {const double tmp[] = {0.,0.,1.}; cols.push_back(makeVector(tmp));} // S1 col
        {const double tmp[] = {0.,1.,0.}; cols.push_back(makeVector(tmp));} // M1 col
        {const double tmp[] = {1.,0.,1.}; cols.push_back(makeVector(tmp));} // A1 col
        {const double tmp[] = {1.,1.,0.}; cols.push_back(makeVector(tmp));}


        for(int j=0;j<nMaps;j++){
            int ioff = j*nLocations;

            // V1, M1, S1 as RGB
            displays[0].resetDisplay(fix,fix,fix);
            for(int i=0;i<nLocations;i++){
                displays[0].drawHex(X[i]-Xoff,Y[i]-Yoff,0.,0.008,response[ioff+i][outputID[0]],response[ioff+i][outputID[2]],response[ioff+i][outputID[1]]);
            }
            stringstream ss1; ss1<< logpath << "/VMS_RGB_";
            ss1 << j << ".png";
            displays[0].saveImage(ss1.str().c_str());

            // max field id
            displays[0].resetDisplay(fix,fix,fix);
            for(int i=0;i<nLocations;i++){

                vector<double> q;
                for(int k=0;k<outputID.size();k++){
                    q.push_back(response[ioff+i][outputID[k]]);
                }
                int m = getArgmax(q);
                vector<double> col = cols[m];
                displays[0].drawHex(X[i]-Xoff,Y[i]-Yoff,0.,0.008,col[0],col[1],col[2]);
            }
            stringstream ss2; ss2<< logpath << "/MaxOut_";
            ss2 << j << ".png";
            displays[0].saveImage(ss2.str().c_str());


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
                displays[0].drawHex(X[i]-Xoff,Y[i]-Yoff,0.,0.008,col[0],col[1],col[2]);
            }
            stringstream ss3; ss3<< logpath << "/MaxAlign_";
            ss3 << j << ".png";
            displays[0].saveImage(ss3.str().c_str());

            vector<vector<double> > weightmat(P.N+1,vector<double>(P.N+1));
            for(int i=0;i<P.W.size();i++){
                weightmat[P.Pre[i]][P.Post[i]] = P.W[i];
            }

            vector<double> flatweightmat;
            for(int i=0;i<P.N+1;i++){
                for(int j=0;j<P.N+1;j++){
                    flatweightmat.push_back(weightmat[i][j]);
                }
            }

            // log weights
            stringstream fname; fname << logpath << "/weightmat.h5";
            HdfData data(fname.str());
            stringstream ss; ss<<"weightmat";
            data.add_contained_vals (ss.str().c_str(), flatweightmat);

        }
        break;
    }

    default: {
            cout<<"Invalid mode selected"<<endl;
            break;
        }
    }

    logfile<<"Goodbye."<<endl;
    logfile.close();
    return 0;
}
