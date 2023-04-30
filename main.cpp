#include <vector>
#include <iostream>
#include "Cosine_similarity.h"

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

    return 0;
}