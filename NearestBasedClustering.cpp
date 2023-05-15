#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include "euclidean_distance.h"

using namespace std;

// A data point is represented as a vector of double values
typedef vector< pair<double, double> > DataPoint;

// A cluster is represented as a set of indices of data points
typedef set<int> Cluster;

// A graph is represented as an adjacency list
typedef vector< vector<int> > Graph;

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

// Builds a graph based on the RNC algorithm
void buildGraph(const vector<DataPoint>& data, double r, Graph& graph) {
    int n = data.size();
    graph.resize(n);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = distance(data[i], data[j]);
            if (d <= r) {
                graph[i].push_back(j);
                graph[j].push_back(i);
            }
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
            // cout<<"Point "<<i<<" is not core point"<< endl;
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

    // Parameters
    double r = 5.5;

    // Clustering
    Graph graph;
    buildGraph(data, r, graph);
    vector<Cluster> clusters;
    assignClusters(graph, clusters);

    // Output
    // printClusters(data, clusters);

    int i = 1;
    for (const Cluster& cluster : clusters) {
        cout<<"Size cluster "<<i<<" : "<<cluster.size()<<endl;
        i++;
        for(int j : cluster){
            cout<<"Index Data : "<<j<<endl;
        }
    }

    return 0;
}
