#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include <time.h>

using namespace std;

// A data point is represented as a vector of double values
typedef vector< pair<double, double> > DataPoint;

// A cluster is represented as a set of indices of data points
typedef set<int> Cluster;

// A graph is represented as an adjacency list
typedef vector< vector<int> > Graph;

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
void buildGraph(const vector<DataPoint>& data, double r, Graph& graph) {
    vector<DataPoint> norm_data = normalizeData(data);

    int n = data.size();
    graph.resize(n);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = distance(data[i], data[j]);
            if (d <= r) {
                graph[i].push_back(j);
                graph[j].push_back(i);
            } 
            // else {
            //     cout<<"The distance is greater than the radius"<<endl;
            // }
        }
    }
}

// Assigns data points to clusters based on the graph
void assignClusters(const Graph& graph, vector<Cluster>& clusters) {
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

int main() {
    // Sample data
    vector<DataPoint> data;
    DataPoint data_point1, data_point2, data_point3, data_point4, data_point5, data_point6, data_point7, data_point8, data_point9, data_point10, data_point11, data_point12, data_point13, data_point14, data_point15;
    data_point1.push_back(pair<double, double>(1.1, 2.3));
    data.push_back(data_point1);
    data_point2.push_back(pair<double, double>(2.2, 1.9));
    data.push_back(data_point2);
    data_point3.push_back(pair<double, double>(2.3, 3.5));
    data.push_back(data_point3);
    data_point4.push_back(pair<double, double>(3.7, 3.8));
    data.push_back(data_point4);
    data_point5.push_back(pair<double, double>(4.7, 4.1));
    data.push_back(data_point5);
    data_point6.push_back(pair<double, double>(5.6, 5.1));
    data.push_back(data_point6);
    data_point7.push_back(pair<double, double>(5.9, 6.5));
    data.push_back(data_point7);
    data_point8.push_back(pair<double, double>(5.0, 5.0));
    data.push_back(data_point8);
    data_point9.push_back(pair<double, double>(7.0, 7.0));
    data.push_back(data_point9);
    data_point10.push_back(pair<double, double>(8.0, 8.0));
    data.push_back(data_point10);
    data_point11.push_back(pair<double, double>(8.6, 9.0));
    data.push_back(data_point11);
    data_point12.push_back(pair<double, double>(9.0, 8.0));
    data.push_back(data_point12);
    data_point13.push_back(pair<double, double>(15.0, 14.0));
    data.push_back(data_point13);
    data_point14.push_back(pair<double, double>(15.1, 14.4));
    data.push_back(data_point14);
    data_point15.push_back(pair<double, double>(15.9, 14.2));
    data.push_back(data_point15);

    // parameter for measure the computation time
    time_t start = clock();

    // Parameters
    double r = 5.2;

    // Clustering
    Graph graph;
    buildGraph(data, r, graph);
    vector<Cluster> clusters;
    assignClusters(graph, clusters);

    // Output
    // printClusters(data, clusters);

    time_t end = clock();
    double elapsed = double(end - start)/ CLOCKS_PER_SEC;
    
    vector<int> pred_clust(data.size(), 0);
    vector<int> true_clust;

    int i = 1;
    for (const Cluster& cluster : clusters) {
        int cluster_size = cluster.size();
        cout<<"Size cluster "<<i<<" : "<<cluster_size<<endl;
        if(cluster_size>0){
            for(int k : cluster){
                pred_clust[k] = i;
            }
        }
        i++;
        for(int j : cluster){
            cout<<"Index Data : "<<j<<endl;
        }
    }

    cout<<"Predicted Cluster : "<<endl;
    for (int i = 0; i < pred_clust.size(); ++i) {
        cout<<"Index Data "<<i<<" is cluster "<< pred_clust[i] << endl;
    }

    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(2);
    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(2);
    true_clust.push_back(2);
    true_clust.push_back(1);
    true_clust.push_back(2);
    true_clust.push_back(3);
    true_clust.push_back(3);
    true_clust.push_back(3);

    //Using RAND index for the clustering evaluation
    double rand_index = calculateRandIndex(pred_clust, true_clust);
    cout<<"RAND index value : "<<rand_index<<endl;
    
    cout<<"Computation Time : "<<elapsed<<" seconds"<<endl;

    return 0;
}
