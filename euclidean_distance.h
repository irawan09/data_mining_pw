#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class EuclideanDistance{
    public :

        vector<double> normalize(vector<double> v) {
            double length = 0.0;
            for (double d : v) {
                length += d * d;
            }
            length = sqrt(length);
            for (double& d : v) {
                d /= length;
            }
            return v;
        }

        double euclidean_distance(const vector<double>& v1, const vector<double>& v2) {
            if (v1.size() != v2.size()) {
                throw runtime_error("Vectors must have the same dimensionality");
            }
            
            double distance = 0.0;
            
            for (size_t i = 0; i < v1.size(); ++i) {
                distance += pow(v1[i] - v2[i], 2);
            }
            
            return sqrt(distance);
        }

};

