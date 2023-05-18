#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>

using namespace std;

class NearestBasedClustering{
    public :

    // Function to normalize data using min-max normalization
    vector<DataPoint> normalizeData(const vector<DataPoint>& data) {
        vector<DataPoint> normalizedData;
        vector<double> x_vect, y_vect; 
        double minValX, maxValX, minValY, maxValY;
        DataPoint normalizedVec;

        // extracting the vector data
        for (const auto& vec : data) {
            double x = vec[0].first;
            x_vect.push_back(x);
            double y = vec[0].second;
            y_vect.push_back(y);
        }

        // Find the minimum and maximum values within each vector
        minValX = *min_element(x_vect.begin(), x_vect.end());
        maxValX = *max_element(x_vect.begin(), x_vect.end());
        minValY = *min_element(y_vect.begin(), y_vect.end());
        maxValY = *max_element(y_vect.begin(), y_vect.end());
        
        // Normalize each data point within the vector
        for (const auto& vec : data) {        
            for (const auto& pair : vec) {
                double normalizedValX = (pair.first - minValX) / (maxValX - minValX);
                double normalizedValY = (pair.second - minValY) / (maxValY - minValY);
                // cout<<"Normalized Value : "<<normalizedValY<<endl;
                normalizedVec.push_back(make_pair(normalizedValX, normalizedValY));
                // cout<<normalizedVec.size()<<endl;
            }
        }
        normalizedData.push_back(normalizedVec);
        // cout<<normalizedData.size()<<endl;
    
    return normalizedData;
}

// Returns the Euclidean distance between two data points
double distance(const DataPoint& p1, const DataPoint& p2) {
    double d = 0;
    for (int i = 0; i < p1.size(); i++) {
        double x = p2[i].first - p1[i].first;
        double y = p2[i].second - p1[i].second;
        d += pow(x, 2) + pow(y, 2);
        // cout<<"Euclidean Distance between ("<<p1[i].first<<","<<p1[i].second<<") and ("<<p2[i].first<<","<<p2[i].second<<"): "<<sqrt(d)<< endl;
    }
    return sqrt(d);
}

// Builds a graph based on the NBC algorithm
vector< vector<int> > buildGraph(const vector<DataPoint>& data, double r, vector< vector<int> >& graph) {
    vector< vector<int> > temp_graph;
    vector<DataPoint> norm_data = normalizeData(data);

    int n = data.size();
    temp_graph.resize(n);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = distance(data[i], data[j]);
            if (d <= r) {
                temp_graph[i].push_back(j);
                temp_graph[j].push_back(i);
            } 
            // else {
            //     cout<<"The distance is greater than the radius"<<endl;
            // }
        }
    }

    return temp_graph;
}

// Assigns data points to clusters based on the graph
vector< set<int> > assignClusters(const vector< vector<int> >& graph, vector< set<int> >& clusters) {
    // Visualize the graph
    // for(int l=0; l< graph.size(); l++){
    //     cout<<"Node number "<<l<<" : ";
    //     for(int m=0; m<graph[l].size(); m++){
    //         cout<<graph[l][m]<<", " ;
    //     }
    //     cout<< endl;
    // }

    int n = graph.size();
    vector<bool> visited(n, false);
    for(int i=0; i< graph.size();i++){
        int size_knn = graph[i].size();
        
        int value_to_count = i;
        int count_rknn = 0;
        for (const auto& inner_vec : graph) {
            for (const auto& value : inner_vec) {
                if (value == value_to_count) {
                    count_rknn++;
                }
            }
        }
        // cout << "Value " << value_to_count << " appears " << count_rknn << " times in nested vector" << endl;

        double ndf;
        ndf = count_rknn/size_knn;
        if( ndf >= 1){
            // cout<<"Point "<<i<<" is core point"<< endl;
            if(!visited[i]){
                visited[i] = true;
                Cluster cluster;
                cluster.insert(i);
                visited[i] = true;
                for (double j : graph[i]) {
                    if (!visited[j]) {
                        visited[j] = true;
                        cluster.insert(j);
                    }
                }

                // for (auto j : cluster){
                //     cout<<"cluster no "<<i<<" : "<<j<< endl; 
                // }

                clusters.push_back(cluster);
            }    
        } 
        // else{
        //     cout<<"Point "<<i<<" is not core point"<< endl;
        // }
    }
}

// Prints the clusters
void printClusters(const vector<DataPoint>& data, const vector<Cluster>& clusters) {
    int i = 1;
    for (const Cluster& cluster : clusters) {
        cout << "Cluster " << i++ << ":" << endl;
        for (int j : cluster) {
            cout << "(";
            for (pair<double, double> x : data[j]) {
                cout << x.first;
            }
            cout << ");" << endl;
        }
    }
}

// Function to calculate the Rand index
double calculateRandIndex(const vector<int>& clusterLabels, const vector<int>& trueLabels) {
    int n = clusterLabels.size();
    int a = 0, b = 0, c = 0, d = 0;

    // Make sure the input vectors have the same size
    if (clusterLabels.size() != trueLabels.size()) {
        cerr << "Error: true and predicted labels have different sizes!" << endl;
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (clusterLabels[i] == clusterLabels[j] && trueLabels[i] == trueLabels[j]) {
                a++;
            } else if (clusterLabels[i] != clusterLabels[j] && trueLabels[i] != trueLabels[j]) {
                b++;
            } else if (clusterLabels[i] == clusterLabels[j] && trueLabels[i] != trueLabels[j]) {
                c++;
            } else if (clusterLabels[i] != clusterLabels[j] && trueLabels[i] == trueLabels[j]) {
                d++;
            }
        }
    }

    double randIndex = static_cast<double>(a + b) / (a + b + c + d);
    return randIndex;
}


}; 



