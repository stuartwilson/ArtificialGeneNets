#include <vector>
#include <math.h>

using namespace std;
class Pineda{
public:
    int N, Nweight, Nplus1, Nouts, Nins, maxConvergenceSteps;
    vector<double> W, X, Input, Target, U, Wbest;
    vector<int> Pre, Post, inputID, outputID;
    double zero, etax2, weightNudgeSize, divergenceThreshold;
    vector<double*> Wptr;

    Pineda(int N, vector<int> inputID, vector<int> outputID, double eta, double weightNudgeSize, double divergenceThreshold, int maxConvergenceSteps){
        this->N=N;
        X.resize(N,0.);
        U.resize(N,0.);
        this->inputID = inputID;
        this->outputID = outputID;
        etax2 = 2.*eta;
        zero = 0.0;
        Nplus1 = N; // overwrite if bias
        this->weightNudgeSize=weightNudgeSize;
        this->divergenceThreshold=divergenceThreshold;
        this->maxConvergenceSteps=maxConvergenceSteps;
    }

    void addBias(void){
        for(int i=0;i<N;i++){
            W.push_back(0.);
            Pre.push_back(N);
            Post.push_back(i);
        }
        X.push_back(1.0);
        Nplus1 = N+1;
    }

    void connect(int pre, int post){
        W.push_back(0.);
        Pre.push_back(pre);
        Post.push_back(post);
    }

    void randomizeWeights(double weightMin, double weightMax){
        double weightRange = weightMax-weightMin;
        for(int i=0;i<W.size();i++){
            W[i] = morph::Tools::randDouble()*weightRange+weightMin;
        }
    }

    void setNet(void){
        Nweight = W.size();
        Nouts = outputID.size();
        Nins = inputID.size();
        Wbest = W;
        Wptr.resize(Nplus1*Nplus1,&zero);
        for(int i=0;i<Nweight;i++){
            Wptr[Pre[i]*Nplus1+Post[i]] = &W[i];
        }
    }

    void reset(vector<double> input, vector<double> target){

        std::fill(X.begin(),X.end(),0.); // NOTE: DW randomized the initial state
        Input = input;
        Target = target;
    }


    void step(void){

        std::fill(U.begin(),U.end(),0.);

        for(int k=0;k<Nweight;k++){
            U[Post[k]] += X[Pre[k]] * W[k];
        }
        for(int i=0;i<N;i++){
            X[i] = 1./(1.+exp(-U[i]));
        }
        for(int i=0;i<Nins;i++){
            X[inputID[i]] += Input[i];
        }
    }

    void step(int ko){
        X[ko] = 0.;
        step();
    }

    double getError(void){
        double error = 0.;
        for(int i=0;i<Nouts;i++){
            error += (Target[i]-X[outputID[i]])*(Target[i]-X[outputID[i]]);
        }
        return error * 0.5;
    }

    vector<double> getOutput(void){
        vector<double> outp(Nouts);
        for(int i=0;i<Nouts;i++){
            outp[i] = X[outputID[i]];
        }
        return outp;
    }

    vector<double> getWeightMatrix(void){
        vector<double> flatweightmat(Wptr.size());
        for(int i=0;i<Wptr.size();i++){
            flatweightmat[i] = *Wptr[i];
        }
        return flatweightmat;
    }


    void converge(int ko){
        vector<double> Xpre(Nweight,0.);
        double total = 1.;
        int count = 0;
        bool knockout = false;
        if(ko>=0){ // -1 is no ko flag
            knockout = true;
        }
        while(total>divergenceThreshold*N){
            count++;
            if(count>maxConvergenceSteps){
                W = Wbest;
                for(int k=0;k<Nweight;k++){
                    W[k] += (morph::Tools::randDouble()*2-1)*weightNudgeSize;
                }
                count = 0;
                break;
            }
            Xpre=X;
            if(knockout){
                step(ko);
            } else {
                step();
            }
            total=0;
            for(int i=0;i<N;i++){
                total +=(X[i]-Xpre[i])*(X[i]-Xpre[i]);
            }
        }
        X[ko]=0.;
    }

    void weightUpdate(void){

        vector<bool> fixed(N,false);
        for(int i=0;i<Nouts;i++){
            fixed[outputID[i]] = true;
        }
        vector<double> deltas(N,0.);
        for(int i=0;i<Nouts;i++){
            int j = outputID[i];
            deltas[j] = etax2*(Target[i]-X[j])*X[j]*(1.-X[j]);
        }
        vector<bool> visits(N,false);
        int visitcount=0;
        int test=0;
        while(true){
            test++;
            visitcount=0;
            for(int j=0;j<N;j++){
                if(!visits[j]){
                    visitcount++;
                    if(fixed[j]){
                        for(int i=0;i<N;i++){
                            if(!fixed[i]){
                                deltas[i] += *Wptr[i*(N+1)+j]*deltas[j];
                            }
                        }
                        visits[j]=true;
                    }
                }
            }
            for(int k=0;k<N;k++){
                if(!fixed[k] and deltas[k]!=0){
                    deltas[k] *= X[k]*(1.-X[k]);
                    fixed[k]=true;
                }
            }
            if(visitcount==0){ break; }
            if(test>10){ cout<<"stuck"<<endl; return;}
        }
        for(int k=0;k<Nweight;k++){
            W[k] += X[Pre[k]] * deltas[Post[k]];
        }
    }
};
