#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class CosineSimilarity{

    public:

    // Function to calculate dot product of two vectors
    double dot_product(const vector<double>& a, const vector<double>& b) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }

    // Function to calculate magnitude of a vector
    double magnitude(const vector<double>& a) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * a[i];
        }
        return sqrt(result);
    }

    // Function to calculate cosine similarity between two vectors
    double cosine_similarity(const vector<double>& a, const vector<double>& b) {
        double dot_product_result = dot_product(a, b);
        double magnitude_a = magnitude(a);
        double magnitude_b = magnitude(b);
        return dot_product_result / (magnitude_a * magnitude_b);
    }

};