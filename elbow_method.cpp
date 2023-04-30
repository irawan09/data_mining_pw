#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

double calculate_distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

int main() {
    vector<double> x, y;
    string line;

    // open CSV file
    fstream file ("coba.csv", ios::in);
    if (file.is_open()) {
        // read data from CSV file
        while (getline(file, line)) {
            stringstream ss(line);
            string value;
            getline(ss, value, ',');
            try{
                double number_x = stod(value);
                x.push_back(number_x);
            }catch (const invalid_argument& e) {

            }
            getline(ss, value, ',');
            try{
                double number_y = stod(value);
                y.push_back(number_y);
            }catch (const invalid_argument& e) {

            }
        }
    } else {
        cout << "Failed to open file." << endl;
        return 1;
    }

    int num_clusters = 10;
    vector<double> wcss(num_clusters);

    // calculate WCSS for different number of clusters
    for (int k = 1; k <= num_clusters; k++) {
        double wcss_k = 0;
        vector<double> centroids(k, 0);

        // initialize centroids
        for (int i = 0; i < k; i++) {
            centroids[i] = y[i];
            
        }

        // calculate WCSS for current number of clusters
        for (int i = 0; i < x.size(); i++) {
            double min_distance = INFINITY;

            // find nearest centroid
            for (int j = 0; j < centroids.size(); j++) {
                double distance = calculate_distance(x[i], y[i], x[j], centroids[j]);
                if (distance < min_distance) {
                    min_distance = distance;
                }
            }

            // add squared distance to WCSS
            wcss_k += pow(min_distance, 2);
        }

        // store WCSS for current number of clusters
        wcss[k-1] = wcss_k;
    }

    // find optimal number of clusters using elbow method
    int num_optimal_clusters = 1;
    double elbow = wcss[0];
    for (int k = 1; k < num_clusters; k++) {
        double diff = elbow - wcss[k];
        if (diff < elbow * 0.1) {
            break;
        }
        elbow = wcss[k];
        num_optimal_clusters++;
    }

    cout << "Number of optimal clusters: " << num_optimal_clusters << endl;

    return 0;
}
