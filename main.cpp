#include <vector>
#include <iostream>
#include "cosine_similarity.h"
#include "euclidean_distance.h"

using namespace std;

int main(){

    vector<double> a;
    a.push_back(1.0);
    a.push_back(2.0);
    a.push_back(3.0);
    vector<double> b;
    b.push_back(2.0);
    b.push_back(3.0);
    b.push_back(4.0);

    CosineSimilarity cos;
    double cosine_sim = cos.cosine_similarity(a, b);

    cout << "Cosine similarity: " << cosine_sim << endl;

    EuclideanDistance eucd;
    double distance = eucd.euclidean_distance(a, b);
    
    cout << "The Euclidean distance between the vectors is " << distance << endl;

    cout << "Original vector: ";
    for (double d : a) {
        cout << d << " ";
    }
    cout << endl;

    cout<<"Normalized Vector : "<<"\n";
    vector<double> normalize = eucd.normalize(a);
    for (double d : normalize){
        cout<< d << " ";
    }
    cout<< endl;

    return 0;
}