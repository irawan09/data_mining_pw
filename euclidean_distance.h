#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

class EuclideanDistance{
    public :

        vector<double> normalize(vector<double>& vec) {
            vector<double> x_normalize;
            double min_val = *std::min_element(vec.begin(), vec.end());
            double max_val = *std::max_element(vec.begin(), vec.end());
            double range = max_val - min_val;
            if (range == 0) {
                // Avoid division by zero
                cout<<"the substraction of max and min value will result zero\n";
            }
            for (auto& x : vec) {
                x = (x - min_val) / range;
                x_normalize.push_back(x);
            }
            return x_normalize;
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

