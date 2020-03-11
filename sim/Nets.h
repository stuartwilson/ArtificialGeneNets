/*
 Implementation of recurrent backprop algorithm following Pineda (1987)
 */
#include "morph/HdfData.h"
#include "morph/display.h"
#include "morph/tools.h"
#include "tools.h"
#include "Pineda.h"

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
        double minX, maxX, minY, maxY, minZ, maxZ, minF, maxF;
        int ncols, nrows, outputID, contextID;
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
        Xscaled.resize(N);
        Yscaled.resize(N);
        Zscaled.resize(N);
        Fscaled.resize(N);

        minX = minY = minZ = minF = +1e9;
        maxX = maxY = maxZ = maxF = -1e9;
        for(int i=0;i<N;i++){
            if(X[i]<minX){ minX = X[i]; }
            if(Y[i]<minY){ minY = Y[i]; }
            if(Z[i]<minZ){ minZ = Z[i]; }
            if(F[i]<minF){ minF = F[i]; }
            if(X[i]>maxX){ maxX = X[i]; }
            if(Y[i]>maxY){ maxY = Y[i]; }
            if(Z[i]>maxZ){ maxZ = Z[i]; }
            if(F[i]>maxF){ maxF = F[i]; }
        }
        for(int i=0;i<N;i++){
            Xscaled[i] = (X[i]-minX)/(maxX-minX);
            Yscaled[i] = (Y[i]-minY)/(maxY-minY);
            Zscaled[i] = (Z[i]-minZ)/(maxZ-minZ);
            Fscaled[i] = (F[i]-minF)/(maxF-minF);
        }
        // deduce number of rows and columns
        vector<double> uniqueX = getUnique(X);
        vector<double> uniqueY = getUnique(Y);
        ncols = uniqueX.size();
        nrows = uniqueY.size();
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

    Net(string logpath){
        this->logpath = logpath;
        morph::Tools::createDir (logpath);
        stringstream ss; ss << logpath << "/log.txt";
        logfile.open(ss.str());
        logfile<<"Hello."<<endl;
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
        displays.push_back(Gdisplay(600, 600, 0, 0, "Image", 1.3, 0.0, 0.0));
        displays[0].resetDisplay(fix,fix,fix);
        displays[0].redrawDisplay();

        // plot supplied map values
        for(int j=0;j<M.size();j++){

            displays[0].resetDisplay(fix,fix,fix);
            displays[0].resetDisplay(fix,fix,fix);

            for(int i=0;i<M[j].N;i++){
                vector<double> rgb = morph::Tools::getJetColor(M[j].Fscaled[i]);
                displays[0].drawRect(M[j].Xscaled[i]-0.5,M[j].Yscaled[i]-0.5,0.,1./(M[j].ncols-1),1./(M[j].nrows-1),rgb);
            }
            stringstream ss; ss<< logpath << "/map_";
            ss << j << ".png";
            displays[0].saveImage(ss.str().c_str());
            displays[0].redrawDisplay();
        }
    }

};

