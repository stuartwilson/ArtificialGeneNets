#include <vector>
#include <math.h>

using namespace std;
class Pineda{
public:
    int N, Nweight, Nplus1, Nouts, Nins;
    vector<double> W, X, Z, Input, Target, F, Fprime, U, V;
    vector<int> Pre, Post, inputID, outputID;
    double dtOverTauX, dtOverTauZ, dtOverTauW, dt;

    Pineda(int N, vector<int> inputID, vector<int> outputID, double taux, double tauz, double tauw, double dt){
        this->N=N;
        Nplus1 = N+1;
        X.resize(N,0.);
        Z.resize(N,0.);
        F.resize(N,0.);
        U.resize(N,0.);
        V.resize(Nplus1,0.);
        Fprime.resize(N,0.);
        Input.resize(N,0.);
        this->inputID = inputID;
        this->outputID = outputID;
        this->dt = dt;
        dtOverTauX = dt/taux;
        dtOverTauZ = dt/tauz;
        dtOverTauW = dt/tauw;

    }

    void addBias(void){
        for(int i=0;i<N;i++){
            W.push_back(0.);
            Pre.push_back(N);
            Post.push_back(i);
        }
        X.push_back(1.0);
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
    }

    void reset(vector<double> input, vector<double> target){

        std::fill(X.begin(),X.end(),0.);
        std::fill(Z.begin(),Z.end(),0.);

        for(int i=0;i<Nins;i++){
            Input[inputID[i]] = input[i];
        }
        Target = target;

    }

    void forward(void){


        std::fill(U.begin(),U.end(),0.);

        // Don't OMP this loop - buggy!
        for(int k=0;k<Nweight;k++){
            U[Post[k]] += X[Pre[k]] * W[k];
        }
        //#pragma omp parallel for
        for(int i=0;i<N;i++){
            F[i] = 1./(1.+exp(-U[i]));
        }
        //#pragma omp parallel for
        for(int i=0;i<N;i++){
            X[i] +=dtOverTauX* ( -X[i] + F[i] + Input[i] );
        }
    }


    void backward(void){

        //#pragma omp parallel for
        for(int i=0;i<N;i++){
            Fprime[i] = F[i]*(1.0-F[i]);
        }

        std::fill(V.begin(), V.end(),0.);

        //#pragma omp parallel for
        for(int k=0;k<Nweight;k++){
            V[Pre[k]] += Fprime[Post[k]] * W[k] * Z[Post[k]];
        }
        //#pragma omp parallel for
        for(int i=0;i<N;i++){
            Z[i] +=dtOverTauZ * (V[i] - Z[i]);
        }
        //#pragma omp parallel for
        for(int i=0;i<Nouts;i++){
            Z[outputID[i]] +=dtOverTauZ* (Target[i]-X[outputID[i]]);
        }
    }

    void weightUpdate(void){

        //#pragma omp parallel for
        for(int k=0;k<Nweight;k++){
            W[k] +=dtOverTauW* (X[Pre[k]] * Z[Post[k]] * Fprime[Post[k]]);
        }
    }

    double getError(void){
        double error = 0.;
        double diff;
        for(int i=0;i<Nouts;i++){
            diff = Target[i] - X[outputID[i]];
            error += diff*diff;
        }
        return error * 0.5;
    }

    vector<double> getOutput(void){
        vector<double> out(Nouts);
        for(int i=0;i<Nouts;i++){
            out[i] = X[outputID[i]];
        }
        return out;
    }

    vector<double> getWeightMatrix(void){
        vector<vector<double> > weightmat(N+1,vector<double>(N+1));
        for(int i=0;i<Nweight;i++){
            weightmat[Pre[i]][Post[i]] = W[i];
        }

        vector<double> flatweightmat;
        for(int i=0;i<Nplus1;i++){
            for(int j=0;j<Nplus1;j++){
                flatweightmat.push_back(weightmat[i][j]);
            }
        }
        return flatweightmat;
    }

};

