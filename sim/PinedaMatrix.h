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
        Nplus1 = N; // overwrite if bias
        this->weightNudgeSize= weightNudgeSize;
        this->divergenceThreshold= divergenceThreshold * N;
        this->maxConvergenceSteps= maxConvergenceSteps;
        zero = 0.0;
    }

    ~Pineda(void){

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

    void randomizeState(void){
        for(int i=0;i<N;i++){
            X[i] = morph::Tools::randDouble()*2.0-1.0;
        }
    }

    void reset(vector<double> input, vector<double> target){

        std::fill(X.begin(),X.end(),0.); // NOTE: DW randomized the initial state
        //randomizeState();
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


    /* ORIGINAL - BASED ON DAN
    void converge(int ko){

        vector<double> Xpre(Nweight,0.);
        double total = 1.;
        int count = 0;
        bool knockout = false;
        if(ko>=0){ // -1 is no ko flag
            knockout = true;
        }

        while(total>divergenceThreshold){
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
        if(knockout){
            X[ko]=0.;
        }
    }

*/

/*
    void converge(int ko){

        vector<double> Xpre(Nweight,0.);
        double total = 1.;
        int count = 0;
        bool knockout = false;
        if(ko>=0){ // -1 is no ko flag
            knockout = true;
        }

        while(total>divergenceThreshold){
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
        if(knockout){
            X[ko]=0.;
        }
    }


        void convergeNoNudge(int ko){

        vector<double> Xpre(Nweight,0.);
        double total = 1.;
        int count = 0;
        bool knockout = false;
        if(ko>=0){ // -1 is no ko flag
            knockout = true;
        }

        while(total>divergenceThreshold){
            count++;
            if(count>maxConvergenceSteps){
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
        if(knockout){
            X[ko]=0.;
        }
    }
*/



    void converge(int ko, bool nudge){
        bool knockout = (ko>=0);
        vector<double> Xpre(N,0.);
        double total = N;
        for(int t=0;t<maxConvergenceSteps;t++){
            if(total>divergenceThreshold){
                Xpre=X;
                if(knockout){ step(ko); } else { step(); }
                total = 0.0;
                for(int i=0;i<N;i++){
                    total +=(X[i]-Xpre[i])*(X[i]-Xpre[i]);
                }
            } else {
                if(nudge){
                    W = Wbest;
                    for(int k=0;k<Nweight;k++){
                        W[k] += (morph::Tools::randDouble()*2-1)*weightNudgeSize;
                    }
                }
                break;
            }
        }
    if(knockout){ X[ko]=0.; }
    }


/*
    void converge(int ko){

        bool knockout = (ko>=0);

        vector<double> Xpre(N,0.);

        for(int t=0;t<maxConvergenceSteps;t++){

            Xpre=X;
            if(knockout){ step(ko); } else { step(); }
            double total = 0.0;
            for(int i=0;i<N;i++){
                total +=(X[i]-Xpre[i])*(X[i]-Xpre[i]);
            }

            if(total<divergenceThreshold){
                W = Wbest;
                for(int k=0;k<Nweight;k++){
                    W[k] += (morph::Tools::randDouble()*2-1)*weightNudgeSize;
                }
                break;
            }
        }

        if(knockout){
            X[ko]=0.;
        }
    }
*/

    void weightUpdate(void){

        // declare output nodes as fixed
        vector<bool> fixed(N,false);
        for(int i=0;i<Nouts;i++){
            fixed[outputID[i]] = true;
        }

        // give output nodes additional delta term
        vector<double> deltas(Nplus1,0.);
        for(int i=0;i<Nouts;i++){
            int j = outputID[i];
            deltas[j] = etax2*(Target[i]-X[j])*X[j]*(1.-X[j]);
        }

        int test=0;
        bool someNodesFixed = false;
        while(!someNodesFixed){

            for(int k=0;k<Nweight;k++){
                if(!fixed[Pre[k]] && fixed[Post[k]]){
                    deltas[Pre[k]] += W[k]*deltas[Post[k]];
                }
            }

            for(int k=0;k<N;k++){
                if(!fixed[k] && deltas[k]!=0){
                    deltas[k] *= X[k]*(1.-X[k]);
                    fixed[k]=true;
                }
            }

            someNodesFixed = false;
            for(int i=0;i<N;i++){
                if(fixed[i]){
                    someNodesFixed = true;
                    break;
                }
            }

            if(test>10){
                cout<<"stuck"<<endl;
                return;
            }
            test++;
        }

        for(int k=0;k<Nweight;k++){
            W[k] += X[Pre[k]] * deltas[Post[k]];
        }
    }



/* //EQUIVALENT TO DAN'S
void weightUpdate(void){

        vector<bool> fixed(N,false);
        for(int i=0;i<Nouts;i++){
            fixed[outputID[i]] = true;
        }
        vector<double> deltas(Nplus1,0.);
        for(int i=0;i<Nouts;i++){
            int j = outputID[i];
            deltas[j] = etax2*(Target[i]-X[j])*X[j]*(1.-X[j]);
        }
        vector<bool> visits(N,false);

        int test=0;

        bool someNodesUnfixed = false;
        while(!someNodesUnfixed){


            for(int k=0;k<Nweight;k++){
                if(!fixed[Pre[k]] && fixed[Post[k]] && !visits[Post[k]]){
                    deltas[Pre[k]] += W[k]*deltas[Post[k]];
                }
            }

            for(int j=0;j<N;j++){
                if(!visits[j] && fixed[j]){
                    visits[j]=true;
                }
            }

            for(int k=0;k<N;k++){
                if(!fixed[k] && deltas[k]!=0){
                    deltas[k] *= X[k]*(1.-X[k]);
                    fixed[k]=true;
                }
            }

            someNodesUnfixed = false;
            for(int i=0;i<N;i++){
                if(visits[i]){
                    someNodesUnfixed = true;
                    break;
                }
            }

            if(test>10){
                cout<<"stuck"<<endl;
                return;
            }
            test++;
        }
        for(int k=0;k<Nweight;k++){
            W[k] += X[Pre[k]] * deltas[Post[k]];
        }
    }
*/

/*
// DAN ORIGINAL
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
*/


};
