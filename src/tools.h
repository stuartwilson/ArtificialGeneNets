#include <vector>
using std::vector;

namespace tools{

    template <class T>
    vector<T> getUnique(vector<T> x){
        vector<T> unique;
        for(int i=0;i<x.size();i++){
            bool uni = true;
            for(int k=0;k<unique.size();k++){
                if(x[i]==unique[k]){ uni = false; break; }
            } if(uni){ unique.push_back(x[i]);}
        }
        return unique;
    }

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

    template <class T>
    T getMin(vector<T> x){
        T min = 1e9;
        for(int i=0;i<x.size();i++){
            if(x[i]<min){
                min = x[i];
            }
        }
        return min;
    }

    template <class T>
    T getMax(vector<T> x){
        T max = -1e9;
        for(int i=0;i<x.size();i++){
            if(x[i]>max){
                max = x[i];
            }
        }
        return max;
    }

    template <class T>
    vector<T> getRenormedVector(vector<T> X){
        T minX = getMin(X);
        T maxX = getMax(X);
        T norm = 1./(maxX-minX);
        for(int i=0;i<X.size();i++){
            X[i] = (X[i]-minX)*norm;
        }
        return X;
    }

};
