#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <tuple>
#include <map>
#include <unordered_map>
#include <cmath>
#include <time.h>
#include <queue>

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
    DataPoint temp_normalizedVec;
    pair<double, double> normalizedPair;

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
            normalizedVec.clear();
            normalizedVec.push_back(make_pair(normalizedValX, normalizedValY));
            normalizedData.push_back(normalizedVec);
        }
    }
    
    return normalizedData;
}

bool hasDuplicate(const vector<double>& distances) {
    set<double> uniqueValues;

    for (double value : distances) {
        if (uniqueValues.find(value) != uniqueValues.end()) {
            // Duplicate value found
            return true;
        }
        uniqueValues.insert(value);
    }

    // No duplicate values found
    return false;
}

unordered_map<double, int> countDuplicates(const vector<double>& distances) {
    unordered_map<double, int> countMap;

    for (double value : distances) {
        countMap[value]++;
    }

    return countMap;
}

double findValueAtIndex2(const vector<tuple<int, int, double>>& data, double target, int sourcePoint) {
    for (const auto& tuple : data) {
        double valueAtIndex2 = get<2>(tuple);
        if (valueAtIndex2 == target) {
            int destinationPoint = get<1>(tuple);
            if(sourcePoint == destinationPoint){
                // cout<<"Index 0"<<endl;
                //
                return std::get<0>(tuple);
            } 
            // cout<<"Index 1"<<endl;
            return  destinationPoint; // Return the value at index 1
        }
    }
    return -1.0;  // Return a default value (-1.0) if target value is not found
}

vector<int> findValueDuplicates(const vector<tuple<int, int, double>>& database, vector<double> targets) {
    vector<int> results;

    for (size_t i = 0; i < database.size(); ++i) {
        double value = get<2>(database[i]);

        for (const auto& target : targets){
            if (value == target) {
                int index0 = get<0>(database[i]);
                results.push_back(index0);
            }
        }
    }

    // Use std::unique to remove adjacent duplicates
    auto last = unique(results.begin(), results.end());

    // Erase the redundant elements from the vector
    results.erase(last, results.end());

    return results;
}

// Remove duplicates the value from the vector
void removeDuplicates(std::vector<int>& vec) {
    std::sort(vec.begin(), vec.end());

    auto uniqueEnd = std::unique(vec.begin(), vec.end());
    vec.erase(uniqueEnd, vec.end());
}


// Returns the Euclidean distance between two data points from vector nested with Pair data
double distanceVectorPair(const DataPoint& p1, const DataPoint& p2) {
    double d = 0;
    for (int i = 0; i < p1.size(); i++) {
        double x = p2[i].first - p1[i].first;
        double y = p2[i].second - p1[i].second;
        d += pow(x, 2) + pow(y, 2);
        // cout<<"Euclidean Distance between ("<<p1[i].first<<","<<p1[i].second<<") and ("<<p2[i].first<<","<<p2[i].second<<"): "<<sqrt(d)<< endl;
    }
    return sqrt(d);
}

// Returns the Euclidean distance between two data points from Pair data
double calculateDistancePair(const pair<double, double>& p1, const pair<double, double>& p2) {
    double dx = p2.first - p1.first;
    double dy = p2.second - p1.second;
    return sqrt(dx * dx + dy * dy);
}

// To find the neareset distance between the target point and all the point in database
vector<double> findNearestDistances(const vector<vector<pair<double, double>>>& data, const pair<double, double>& targetPoint, int k) {
    vector<double> distances;

    for (const auto& innerVector : data) {
        for (const auto& point : innerVector) {
            double distance = calculateDistancePair(targetPoint, point);
            if (distance != 0){
                distances.push_back(distance);
            }
        }
    }

    sort(distances.begin(), distances.end());

    bool hasDuplicates = hasDuplicate(distances);

    if (hasDuplicates) {
        // cout << "There are duplicate values." << endl;
        int number_of_duplicates;

        unordered_map<double, int> duplicateCounts = countDuplicates(distances);

        for (const auto& entry : duplicateCounts) {
            if (entry.second > 1) {
                number_of_duplicates = entry.second - 1;
                 distances.resize(k+number_of_duplicates);
            }
        }
    } else {
        // cout << "There are no duplicate values." << endl;
        if (distances.size() > k) {
            distances.resize(k);
        }
    }

    return distances;
}

vector<bool> calculatingNDF(const Graph& graph){

    int n = graph.size();
    vector<bool> corePoint(n, false);

    for(int i=0; i< graph.size();i++) {
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
            corePoint[i] = true;
         }
    }

    return corePoint;
}

// Builds a graph based on the NBC algorithm
void buildGraph(const vector<DataPoint>& data, int k, Graph& graph) {
    vector<DataPoint> norm_data = normalizeData(data);

    int n = norm_data.size();
    // Graph represents the KNN value connected into the source points
    graph.resize(n);


    // // Normal Euclidean Distance
    // vector<tuple<int, int, double>> database;

    // //create a distance database with source and destination point
    // for (int i = 0; i < n; i++) {
    //     for (int j = i + 1; j < n; j++) {
    //         double d = distanceVectorPair(norm_data[i], norm_data[j]);
    //         database.push_back(make_tuple(i, j, d));
    //     }
    // }

    // vector<double> nearestDistances;
    // vector<int> duplicateIndex;

    // for (int i = 0; i < n; i++) {
    //         pair<double, double> targetPoint;
    //         DataPoint target_data = norm_data[i];

    //         for (const auto& pair : target_data) {
    //             targetPoint = make_pair(pair.first, pair.second);
    //         }

    //         nearestDistances = findNearestDistances(norm_data, targetPoint, k);

    //         // Print the nearest distance
    //         // cout<<"---------------"<<i<<"----------------"<<endl;
    //         // for (const auto& pair : nearestDistances) {
    //         //     cout<<"Nearest Distance : "<<pair<<endl;
    //         // }
    //         // cout<<"--------------------------------"<<endl;


    //         for (const auto& targetValue : nearestDistances) {
    //             double result = findValueAtIndex2(database, targetValue, i);
    //             if (result != -1.0) {

    //                 if (hasDuplicate(nearestDistances)){
    //                     duplicateIndex = findValueDuplicates(database, nearestDistances);

    //                     for(const auto& duplicate : duplicateIndex ){
    //                         if (graph[i].size() < nearestDistances.size()){
    //                             graph[i].push_back(duplicate);
    //                         }
    //                     }             
    //                 } else {
    //                      graph[i].push_back(result);
    //                 }
    //             } else {
    //                 cout << "Target value " << targetValue << " not found." << endl;
    //             }
    //         }
    // }
  

    // Triangle Inequality Property
    double epsilon = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d, dist1, dist2, dist3, distance4;

            d = distanceVectorPair(norm_data[i], norm_data[j]);
            // cout<<"distance data "<<i<<" and data "<<j<<" : "<<d<<endl;

            dist1 = distanceVectorPair(norm_data[i], norm_data[0]);
            // cout<<"Distance point "<<i<<" into reference point : "<<dist1<<endl;

            dist2 = distanceVectorPair(norm_data[1], norm_data[i]);
            // cout<<"Distance point "<<i<<" into reference point : "<<dist2<<endl;

            if(epsilon == 0){
                epsilon = max(dist1, dist2);
                // cout<<"Epsilon value : "<<epsilon<<endl;
            } else {
                dist3 = distanceVectorPair(norm_data[j], norm_data[0]);
                distance4 = dist1 - dist3;
                if (distance4 < epsilon){
                    epsilon = d;
                    // cout<<"Epsilon value : "<<epsilon<<endl;
                    graph[i].push_back(j);
                    graph[j].push_back(i);
                }
            }
        }
    }
}

// Assigns data points to clusters based on the graph
void assignClusters(const Graph& graph, vector<Cluster>& clusters) {

    // // Visualize the graph
    // for(int l=0; l<= graph.size(); l++){
    //     cout<<"Node number "<<l<<" : ";
    //     for(int m=0; m<graph[l].size(); m++){
    //         cout<<graph[l][m]<<", " ;
    //     }
    //     cout<< endl;
    // }

    // Clustering process
    int n = graph.size();

    // Marks whether the point is a core point or not
    vector<bool> corePoint(n, false); 
    // Marks whether the point has been visited
    vector<bool> visited(n, false);
    // Seed for the clustering process
    queue<int> seeds; 

    int clusterID = 1;
    clusters.resize(n);

    corePoint = calculatingNDF(graph);
    
    for (int i = 0; i < graph.size(); i++) {
        if (!visited[i]) {
            visited[i] = true;

            // Process 1: Iterate from index 0 data until you find the core point
            seeds.push(i);
            while (!seeds.empty()) {
            
                int currentNode = seeds.front();
                seeds.pop();

                // Process 2: If you find a core point, give the data point the same cluster as the node being processed
                if (corePoint[currentNode]) {
                    clusters[currentNode].insert(clusterID);
                }

                // Process 3: Fill in the seed with the KNN value of that point
                vector<int> knn = graph[currentNode];

                // Process 4: Take the first seed point that is being processed
                for (int nextNode : knn) {
                    // Process 5: Repeat the second process and so on. until there are no more points in the seed
                    if (!visited[nextNode]) {
                        visited[nextNode] = true;
                        seeds.push(nextNode);
                        clusters[nextNode].insert(clusterID);
                    }
                    clusters[nextNode].insert(clusterID);
                }
            }

            // Process 6: Repeat the first process and so on. until all clusters are filled
            clusterID++;
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
    DataPoint data_point1, data_point2, data_point3, data_point4, data_point5, 
    data_point6, data_point7, data_point8, data_point9, data_point10, 
    data_point11, data_point12;

    data_point1.push_back(pair<double, double>(4.2, 4.0));
    data.push_back(data_point1);
    data_point2.push_back(pair<double, double>(5.9, 3.9));
    data.push_back(data_point2);
    data_point3.push_back(pair<double, double>(2.8, 3.5));
    data.push_back(data_point3);
    data_point4.push_back(pair<double, double>(12.0, 1.3));
    data.push_back(data_point4);
    data_point5.push_back(pair<double, double>(10.0, 1.3));
    data.push_back(data_point5);
    data_point6.push_back(pair<double, double>(1.1, 3.0));
    data.push_back(data_point6);
    data_point7.push_back(pair<double, double>(0.0, 2.4));
    data.push_back(data_point7);
    data_point8.push_back(pair<double, double>(2.4, 2.0));
    data.push_back(data_point8);
    data_point9.push_back(pair<double, double>(11.5, 1.8));
    data.push_back(data_point9);
    data_point10.push_back(pair<double, double>(11.0, 1.0));
    data.push_back(data_point10);
    data_point11.push_back(pair<double, double>(0.9, 0.0));
    data.push_back(data_point11);
    data_point12.push_back(pair<double, double>(1.0, 1.5));
    data.push_back(data_point12);

    // parameter for measure the computation time
    time_t start = clock();

    // Parameters k point
    int k = 3;

    // Clustering
    Graph graph;
    buildGraph(data, k, graph);
    vector<Cluster> clusters;
    assignClusters(graph, clusters);

    time_t end = clock();
    double elapsed = double(end - start)/ CLOCKS_PER_SEC;
    
    vector<int> pred_clust(data.size(), 0);
    vector<int> true_clust;

    for(int i=0; i< clusters.size();i++) {
        set<int> :: iterator iter;
        for (iter = clusters[i].begin(); iter != clusters[i].end(); iter++) {
            cout<< "Node "<<i<<" : cluster "<<*iter << endl;
            pred_clust[i] = *iter;
        }
    }

    true_clust.push_back(1); 
    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(2);
    true_clust.push_back(2);
    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(1);
    true_clust.push_back(2);
    true_clust.push_back(2);
    true_clust.push_back(1);
    true_clust.push_back(1);

    //Using RAND index for the clustering evaluation
    double rand_index = calculateRandIndex(pred_clust, true_clust);
    cout<<"RAND index value : "<<rand_index<<endl;
    cout<<"Computation Time : "<<elapsed<<" seconds"<<endl;

    return 0;
}