/*
 Implementation of recurrent backprop algorithm following Pineda (1987)
 */

// With the correct OpenGL definitions (-DGL3_PROTOTYPES etc) you probably don't need this for Apple
#ifdef __OSX__
# include "OpenGL/gl3.h"
#endif

#include "morph/HdfData.h"
#include "morph/Visual.h"
#include "morph/QuadsVisual.h"
#include "morph/HexGridVisual.h"
#include "morph/ColourMap.h"
#include "morph/tools.h"
#include <morph/Config.h>
#include "morph/Scale.h"
#include "tools.h"
#include "Pineda.h"

#include "morph/ReadCurves.h"
#include "morph/RD_Base.h"

using morph::Config;
using morph::Visual;
using morph::ColourMapType;
using morph::Tools;
using morph::QuadsVisual;
using morph::HexGridVisual;
using morph::Scale;
using morph::HdfData;

using namespace tools;


template <class Flt>
class Domain : public morph::RD_Base<Flt> {
public:
    double ellipseA=1.0;
    double ellipseB=1.0;
    virtual void init (void) {
        this->stepCount = 0;
    }
    void setEllipse(double ellipseA, double ellipseB, double hextohex_d){
        this->ellipseA = ellipseA;
        this->ellipseB = ellipseB;
        this->hextohex_d = hextohex_d;
    }
    virtual void allocate (void) {

        this->hg = new HexGrid (this->hextohex_d, this->hexspan, 0, morph::HexDomainShape::Boundary);
            DBG ("Initial hexagonal HexGrid has " << this->hg->num() << " hexes");
            this->hg->setEllipticalBoundary (ellipseA, ellipseB);
            // Compute the distances from the boundary
            this->hg->computeDistanceToBoundary();
            // Vector size comes from number of Hexes in the HexGrid
            this->nhex = this->hg->num();
            DBG ("After setting boundary, HexGrid has " << this->nhex << " hexes");
            // Spatial d comes from the HexGrid, too.
            this->set_d(this->hg->getd());
            DBG ("HexGrid says d = " << this->d);
            this->set_v(this->hg->getv());
            DBG ("HexGrid says v = " << this->v);
    }
    virtual void step (void) {
        this->stepCount++;
    }
};


class Map{

    /*
      Structure for storing a map (pre-defined X, Y, Z, F, in a HdfData file)
    */

    public:
        int N;
        vector<double> Xscaled, Yscaled, Zscaled, Fscaled;
        vector<double> X, Y, Z, F;
        double minX, maxX, minY, maxY, minZ, maxZ, minF, maxF, xScale, yScale, xSep, xOff, yOff;// ySep;
        vector<double> ySep;
        int outputID, contextID;
        double contextVal;
        vector<array<float, 12>> quads;

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
        Xscaled = tools::normalize(X);
        Yscaled = tools::normalize(Y);
        Zscaled = tools::normalize(Z);
        Fscaled = tools::normalize(F);

        double maxDim = maxX-minX;
        if((maxY-minY)>maxDim){ maxDim = maxY-minY; };
        xScale = (maxX-minX)/maxDim;
        yScale = (maxY-minY)/maxDim;

        xOff = -0.5*(maxX-minX);
        yOff = -0.5*(maxY-minY);

        vector<double> uniqueX = tools::getUnique(X);
        int cols = uniqueX.size();
        vector<int> colID(N,0);
        vector<vector<double> > yByCol(cols);
        vector<int> count(cols,0);
        for(int i=0;i<N;i++){
            for(int j=0;j<cols;j++){
                if(X[i]==uniqueX[j]){
                    colID[i]=j;
                    yByCol[j].push_back(Y[i]);
                    count[j]++;
                }
            }
        }

        vector<double> yRange(cols);
        for(int i=0;i<cols;i++){
            yRange[i] = tools::getMax(yByCol[i])-tools::getMin(yByCol[i]);
        }

        xSep = 0.5*(maxX-minX)/((double)cols-1);

        for(int i=0;i<N;i++){
            ySep.push_back(0.5*yRange[colID[i]]/((double)count[colID[i]]-1));
        }

        array<float, 12> sbox;
        for (int i=0; i<N; i++) {
            // corner 1 x,y,z
            sbox[0] = xScale*(xOff+X[i]-xSep);
            sbox[1] = yScale*(yOff+Y[i]-ySep[i]);
            sbox[2] = 0.0;
            // corner 2 x,y,z
            sbox[3] = xScale*(xOff+X[i]-xSep);
            sbox[4] = yScale*(yOff+Y[i]+ySep[i]);
            sbox[5] = 0.0;
            // corner 3 x,y,z
            sbox[6] = xScale*(xOff+X[i]+xSep);
            sbox[7] = yScale*(yOff+Y[i]+ySep[i]);
            sbox[8] = 0.0;
            // corner 4 x,y,z
            sbox[9] = xScale*(xOff+X[i]+xSep);
            sbox[10]= yScale*(yOff+Y[i]-ySep[i]);
            sbox[11]= 0.0;
            quads.push_back(sbox);
        }

        outputID = -1; // flag for not set
        contextID = -1; // flag for not set
        contextVal = 0.0; // default

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
    Domain<double> domain;
    double xScale, yScale, xOffset, yOffset;
    int mapID, locID;
    vector<int> inputID;
    vector<int> outputID;
    vector<int> contextIDs;
    vector<double> contextVals;
    int nContext;
    morph::ColourMapType colourMap;


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

        float dx = conf.getFloat("dx",0.02);
        float yAspect = conf.getFloat("yAspect",0.75);
        float scaleDomain = conf.getFloat("scaleDomain",1.5);


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

        for(int i=0;i<M.size();i++){
            if(M[i].outputID!=-1){
                outputID.push_back(M[i].outputID);
            }
        }
        outputID = tools::getUnique(outputID);

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

        // Define the domain over which it can be evaluated
        //{ stringstream ss; ss << logpath <<"/domain.svg"; domain.svgpath = ss.str(); }
        domain.init();
        domain.setEllipse(1.0,yAspect,dx);
        domain.allocate();

        double maxX=-1e9;
        double maxY=-1e9;
        double minX= 1e9;
        double minY= 1e9;
        for(int i=0;i<M.size();i++){
            if(M[i].maxX>maxX){maxX=M[i].maxX;};
            if(M[i].maxY>maxY){maxY=M[i].maxY;};
            if(M[i].minX<minX){minX=M[i].minX;};
            if(M[i].minY<minY){minY=M[i].minY;};
        }
        xOffset = minX+(maxX-minX)*0.5;
        yOffset = minY+(maxY-minY)*0.5;
        xScale = (maxX-minX)/(2.0*domain.ellipseA);
        yScale = xScale;

        xScale *= scaleDomain;
        yScale *= scaleDomain;

        // identify the unique context conditions
        vector<int> allContextIDs;
        vector<double> allContextVals;
        for(int i=0;i<M.size();i++){
            if(M[i].contextID!=-1){
                allContextIDs.push_back(M[i].contextID);
                allContextVals.push_back(M[i].contextVal);
            }
        }

        for(int i=0;i<allContextIDs.size();i++){
            bool uni = true;
            for(int k=0;k<contextIDs.size();k++){
                if((allContextIDs[i]==contextIDs[k]) && (allContextVals[i]==contextVals[k])){
                    uni = false; break;
                }
            } if(uni){
                contextIDs.push_back(allContextIDs[i]);
                contextVals.push_back(allContextVals[i]);
            }
        }
        nContext = contextIDs.size();

        setColourMap(morph::ColourMapType::Viridis);

	setMap(0);
	sampleMap(0);

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

    void setInput(void){
        P.reset();
        P.Input[inputID[0]] = M[mapID].X[locID]; // Training on the supplied x-values
        P.Input[inputID[1]] = M[mapID].Y[locID];
	if(M[mapID].contextID != -1){
       	   P.Input[M[mapID].contextID] = M[mapID].contextVal;
	}
 }

    void setRandomInput(void){
        setMap(floor(morph::Tools::randDouble()*M.size()));
        sampleMap(floor(morph::Tools::randDouble()*M[mapID].N));
        setInput();
    }

    vector<vector<double> > testMap(int i){
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

    vector<vector<double> > testDomainContext(int i){

        vector<vector<double> > R(P.N,vector<double>(domain.nhex,0.));

        for (unsigned int j=0; j<domain.nhex; ++j) {
            P.reset();
            P.Input[inputID[0]] = domain.hg->vhexen[j]->x*xScale+xOffset;
            P.Input[inputID[1]] = domain.hg->vhexen[j]->y*yScale+yOffset;
            P.Input[contextIDs[i]] = contextVals[i];
            P.convergeForward(-1,false);
            for(int k=0;k<P.N;k++){
                R[k][j] = P.X[k];
            }
        }
        return R;
    }

    vector<vector<vector<double> > > testDomains(void){

        vector<vector<vector<double> > > R;
        for(int i=0;i<nContext;i++){
            R.push_back(testDomainContext(i));
        }
        return R;
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
	
            if(fmod(k,K/100)==0){
                logfile<<"steps: "<<(int)(100*(float)k/(float)K)<<"% ("<<k<<")"<<endl;
            }
        }
        P.W = P.Wbest;

    }


    /*
        PLOTTING
    */

    void setColourMap(morph::ColourMapType cmap){
        colourMap = cmap;
    }

    void plotMapValues(vector<double> F, string fname, int mapIndex){
        if(M[mapIndex].N != F.size()){ cout<<"Field supplied not correct size (Map)."<<endl;}
        Visual v (500, 500, "Map");
        v.backgroundWhite();
        v.zNear = 0.001;
        v.zFar = 20;
        v.fov = 45;
        v.sceneLocked = false;
        v.setZDefault(-3.7f);
        v.setSceneTransXY (0.0f,0.0f);
        array<float, 3> offset = { 0., 0., 0.0 };
        Scale<float> scale;
        scale.do_autoscale = true;
        vector<float> fFlt;
        for (unsigned int i=0; i<M[mapIndex].N; i++){ fFlt.push_back (static_cast<float>(F[i])); }
        v.addVisualModel (new QuadsVisual<float> (v.shaderprog, &M[mapIndex].quads, offset, &fFlt, scale, colourMap));
        v.render();
        v.render();
        v.saveImage(fname);
    }


    void plotDomainValues(vector<double> F, string fname){
        if(domain.nhex != F.size()){ cout<<"Field supplied not correct size (domain)"<<endl;}
        Visual v (500, 500, "Response");
        v.backgroundWhite();
        v.zNear = 0.001;
        v.zFar = 20;
        v.fov = 45;
        v.sceneLocked = false;
        v.setZDefault(-2.7f);
        v.setSceneTransXY (0.0f,0.0f);
        array<float, 3> offset = { 0.0, 0.0, 0.0 };
        Scale<float> scale;
        scale.do_autoscale = true;
        Scale<float> zscale; zscale.setParams (0.0f, 0.0f);
        vector<float> fFlt;
        for (unsigned int k=0; k<domain.nhex; k++){ fFlt.push_back (static_cast<float>(F[k])); }
        v.addVisualModel (new HexGridVisual<float> (v.shaderprog, domain.hg, offset, &fFlt, zscale, scale, colourMap));
        v.render();
        v.render();
        v.saveImage(fname);
    }


    // **************************************************** //

    void plotMapTargets(void){
        for(int i=0;i<M.size();i++){
            stringstream ss; ss<< logpath << "/targ_map_" << i << ".png";
            plotMapTarget(i, ss.str().c_str());
        }
    }

    void plotMapTarget(int i, string fname){
        plotMapValues(M[i].Fscaled, fname, i);
    }

    void plotMapResponsesAllMaps(void){
        for(int i=0;i<M.size();i++){
            plotMapResponses(i);
        }
    }

    void plotMapResponses(int i){
        vector<vector<double> > R = testMap(i);
        for(int j=0;j<R.size();j++){
            vector<double> F = normalize(R[j]);
            stringstream ss; ss<< logpath << "/resp_map_" << i << "_node_" << j << ".png";
            plotMapValues(F, ss.str().c_str(), i);
        }
    }

    void plotDomainContext(int i){
        vector<vector<double> > R = testDomainContext(i);
        R = tools::normalize(R);
        for(int j=0;j<R.size();j++){
            stringstream ss; ss<< logpath << "/context_" << contextIDs[i] << "_val_" << contextVals[i] << "_Node_" << j << "_(indivNorm).png";
            plotDomainValues(R[j],ss.str().c_str());
        }
    }

    void plotDomainsAllContexts(void){
        vector<vector<vector<double> > > R = testDomains();
        R = tools::normalize(R);
        for(int i=0;i<R.size();i++){
            for(int j=0;j<R[i].size();j++){
                stringstream ss; ss<< logpath << "/context_" << contextIDs[i] << "_val_" << contextVals[i] << "_Node_" << j << "_(jointNorm).png";
                plotDomainValues(R[i][j],ss.str().c_str());
            }
        }
    }

    void plotDomainNodeDiff(int contextIndex, int nodeA, int nodeB){

        if(contextIndex>=nContext){cout<<"Invalid context ID "<<contextIndex<<". Only "<<nContext<<"contexts."<<endl;}
        if(nodeA>=P.N){cout<<"Invalid node ID (A) "<<nodeA<<". Only "<<P.N<<"nodes."<<endl;}
        if(nodeB>=P.N){cout<<"Invalid node ID (B) "<<nodeB<<". Only "<<P.N<<"nodes."<<endl;}

        vector<vector<double> >  A = testDomainContext(contextIndex);
        vector<double> diff = A[nodeA];
        for(int i=0;i<domain.nhex;i++){
            diff[i] -= A[nodeB][i];
        }
        diff = tools::normalize(diff);

        stringstream ss; ss<< logpath << "/DIFF_node_"<<nodeA<<"_minus_node_"<<nodeB<<"_context_" << contextIDs[contextIndex] << "_val_" << contextVals[contextIndex]<<".png";
        plotDomainValues(diff,ss.str().c_str());

    }

    void plotDomainContextDiff(int nodeIndex, int contextA, int contextB){

        if(contextA>=nContext){cout<<"Invalid context ID (A) "<<contextA<<". Only "<<nContext<<"contexts."<<endl;}
        if(contextB>=nContext){cout<<"Invalid context ID (B) "<<contextB<<". Only "<<nContext<<"contexts."<<endl;}
        if(nodeIndex>=P.N){cout<<"Invalid node ID "<<nodeIndex<<". Only "<<P.N<<"nodes."<<endl;}

        vector<vector<double> >  A = testDomainContext(contextA);
        vector<vector<double> >  B = testDomainContext(contextB);
        vector<double> diff = A[nodeIndex];
        for(int i=0;i<domain.nhex;i++){
            diff[i] -= B[nodeIndex][i];
        }
        diff = tools::normalize(diff);

        stringstream ss; ss<< logpath << "/DIFF_context_("<<contextIDs[contextA]<<","<<contextVals[contextA]<<")_minus_context_("<<contextIDs[contextB]<<","<<contextVals[contextB]<<")_node"<<nodeIndex<<".png";
        plotDomainValues(diff,ss.str().c_str());

    }

    void plotDomainContextDiffOutputNodes(int contextA, int contextB){
        for(int i=0;i<outputID.size();i++){
            plotDomainContextDiff(outputID[i], contextA, contextB);
        }

    }

};
