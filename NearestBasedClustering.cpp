#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <time.h>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace std;

typedef vector<vector<double>> DataFromCSV;

// A graph is represented as an adjacency list
typedef vector<vector<int>> Graph;

// Normal Euclidean Distance database
vector<tuple<int, int, double>> database;

vector<bool> corePointStatus;

int count_rknn;

vector<vector<double>> readCSVFile(const string &filename)
{
  vector<vector<double>> columns;
  ifstream file(filename);

  if (!file.is_open())
  {
    cout << "Failed to open the CSV file" << endl;
    return columns;
  }

  string line;
  getline(file, line); // Read the first row (labels) and discard

  while (getline(file, line))
  {
    stringstream ss(line);
    string cell;
    size_t columnIndex = 0;

    while (getline(ss, cell, ','))
    {
      // If it's the first row, add a new column vector
      if (columns.size() <= columnIndex)
      {
        columns.push_back(vector<double>());
      }

      // Convert the cell value to a double and add it to the respective column
      // vector
      double value = stod(cell);
      columns[columnIndex].push_back(value);
      columnIndex++;
    }
  }

  return columns;
}

vector<vector<double>> normalizeData(const vector<vector<double>> &data)
{
  size_t numColumns = data.size();
  size_t numRows = data[0].size();

  vector<vector<double>> normalizedColumns(numColumns,
                                           vector<double>(numRows, 0.0));

  for (size_t i = 0; i < numColumns; i++)
  {
    double magnitude = 0.0;

    // Calculate the magnitude of the vector
    for (size_t j = 0; j < numRows; j++)
    {
      magnitude += pow(data[i][j], 2);
    }

    magnitude = sqrt(magnitude);

    // Normalize each value in the column and convert to double
    for (size_t j = 0; j < numRows; j++)
    {
      normalizedColumns[i][j] = static_cast<double>(data[i][j]) / magnitude;
    }
  }

  return normalizedColumns;
}

double calculateEuclideanDistance(const vector<double> &point1,
                                  const vector<double> &point2)
{

  if (point1.size() != point2.size())
  {
    throw invalid_argument("Vectors must have the same size");
  }

  double distance = 0.0;
  size_t numFeatures = point1.size();

  for (size_t i = 0; i < numFeatures; i++)
  {
    double diff = pow(point1[i], 2) - pow(point2[i], 2);
    distance += diff * diff;
  }

  return sqrt(distance);
}

bool hasDuplicate(const vector<double> &values)
{
  unordered_map<double, int> countMap;
  for (const auto &value : values)
  {
    countMap[value]++;
    if (countMap[value] > 1)
    {
      return true;
    }
  }
  return false;
}

unordered_map<double, int> countDuplicates(const vector<double> &values)
{
  unordered_map<double, int> countMap;
  for (const auto &value : values)
  {
    countMap[value]++;
  }
  return countMap;
}

vector<double> extractColumn(const vector<vector<double>> &data,
                             size_t column_index)
{
  vector<double> column;

  if (!data.empty())
  {
    if (column_index < data[0].size())
    {
      for (const auto &row : data)
      {
        column.push_back(row[column_index]);
      }
    }
    else
    {
      throw invalid_argument("Invalid column index");
    }
  }
  else
  {
    throw invalid_argument("Empty dataset");
  }

  return column;
}

vector<double> findNearestDistances(const vector<vector<double>> &data,
                                    const vector<double> &targetPoint,
                                    size_t k)
{
  vector<double> distances;

  try
  {
    for (int i = 0; i < data[0].size(); i++)
    {
      vector<double> innerVector = extractColumn(data, i);
      double distance = calculateEuclideanDistance(targetPoint, innerVector);
      if (distance != 0)
      {
        distances.push_back(distance);
      }
    }
  }
  catch (const exception &e)
  {
    cout << "Error: " << e.what() << endl;
  }

  sort(distances.begin(), distances.end());

  bool hasDuplicates = hasDuplicate(distances);

  if (hasDuplicates)
  {
    int number_of_duplicates;

    unordered_map<double, int> duplicateCounts = countDuplicates(distances);

    for (const auto &entry : duplicateCounts)
    {
      if (entry.second > 1)
      {
        number_of_duplicates = entry.second - 1;
        distances.resize(k + number_of_duplicates);
      }
    }
  }
  else
  {
    if (distances.size() > k)
    {
      distances.resize(k);
    }
  }

  return distances;
}

double findValueAtIndex2(const vector<tuple<int, int, double>> &data,
                         double target, int sourcePoint)
{
  for (const auto &tuple : data)
  {
    double valueAtIndex2 = get<2>(tuple);
    if (valueAtIndex2 == target)
    {
      int destinationPoint = get<1>(tuple);
      if (sourcePoint == destinationPoint)
      {
        return get<0>(tuple);
      }
      
      return destinationPoint;
    }
  }
  
  return -1.0;
}

vector<int> findValueDuplicates(const vector<tuple<int, int, double>> &database,
                                vector<double> targets)
{
  vector<int> results;

  for (size_t i = 0; i < database.size(); ++i)
  {
    double value = get<2>(database[i]);

    for (const auto &target : targets)
    {
      if (value == target)
      {
        int index0 = get<0>(database[i]);
        results.push_back(index0);
      }
    }
  }

  auto last = unique(results.begin(), results.end());

  results.erase(last, results.end());

  return results;
}

void buildGraph(const vector<vector<double>> &data, int k, Graph &graph)
{
  vector<vector<double>> norm_data = normalizeData(data);

  // for(int i=0;i<norm_data.size();i++){
  //     cout<<"Normalize data index "<<i<<" ";
  //     for(int j=0;j<norm_data[i].size();j++){
  //         cout<<norm_data[i][j]<<" ";
  //     }
  //     cout<<endl;
  // }

  size_t n = norm_data[0].size();

  // Graph represents the KNN value connected into the source points
  graph.resize(n);

  // create a distance database with source and destination point
  try
  {
    for (size_t i = 0; i < n; i++)
    {
      for (size_t j = i + 1; j < n; j++)
      {
        vector<double> row1 = extractColumn(norm_data, i);
        vector<double> row2 = extractColumn(norm_data, j);
        double d = calculateEuclideanDistance(row1, row2);
        // cout<<"Distance of "<<i<<" towards "<<j<<" : "<<d<<endl;
        database.push_back(make_tuple(i, j, d));
      }
    }
  }
  catch (const exception &e)
  {
    cout << "Error: " << e.what() << endl;
  }

  vector<double> nearestDistances;
  vector<int> duplicateIndex;

  for (int i = 0; i < n; i++)
  {
    vector<double> targetData = extractColumn(norm_data, i);

    nearestDistances = findNearestDistances(norm_data, targetData, k);

    // cout<<"Node "<<i<<" : ";
    // for(int j=0;j< nearestDistances.size();j++) {
    //     cout<<nearestDistances[j]<<" ";
    // }
    // cout<<endl;

    for (const auto &targetValue : nearestDistances)
    {
      double result = findValueAtIndex2(database, targetValue, i);

      if (result != -1.0)
      {
        if (hasDuplicate(nearestDistances))
        {
          duplicateIndex = findValueDuplicates(database, nearestDistances);

          for (const auto &duplicate : duplicateIndex)
          {
            if (graph[i].size() < nearestDistances.size())
            {
              graph[i].push_back(duplicate);
            }
          }
        }
        else
        {
          graph[i].push_back(result);
        }
      }
      else
      {
        cout << "Target value " << targetValue << " not found." << endl;
      }
    }
  }
}

vector<bool> calculatingNDF(const Graph &graph)
{

  int n = graph.size();
  vector<bool> corePoint(n, false);

  for (int i = 0; i < graph.size(); i++)
  {
    int size_knn = graph[i].size();

    int value_to_count = i;
    count_rknn = 0;
    for (const auto &inner_vec : graph)
    {
      for (const auto &value : inner_vec)
      {
        if (value == value_to_count)
        {
          count_rknn++;
        }
      }
    }

    double ndf;
    ndf = count_rknn / size_knn;
    if (ndf >= 1)
    {
      corePoint[i] = true;
    }
  }

  return corePoint;
}

// Assigns data points to clusters based on the graph
void assignClusters(const Graph &graph, vector<vector<int>> &clusters)
{

  // for(int i=0; i<graph.size();i++){
  //     cout<<"Node "<<i<<" : ";
  //     for(int j=0;j<graph[0].size();j++){
  //         cout<<" "<<graph[i][j]<<", ";
  //     }
  //     cout<<endl;
  // }

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

  for (int i = 0; i < graph.size(); i++)
  {
    if (!visited[i])
    {
      visited[i] = true;

      // Process 1: Iterate from index 0 data until you find the core point
      seeds.push(i);
      while (!seeds.empty())
      {

        int currentNode = seeds.front();

        seeds.pop();

        // Process 2: If you find a core point, give the data point the same
        // cluster as the node being processed
        if (corePoint[currentNode])
        {
          if (clusters[currentNode].size() < 1)
          {
            clusters[currentNode].push_back(clusterID);
          }
        }

        // cout<<"Node "<<currentNode<<" : Cluster "<<clusterID<<endl;

        // Process 3: Fill in the seed with the KNN value of that point
        vector<int> knn = graph[currentNode];

        // Process 4: Take the first seed point that is being processed
        for (int nextNode : knn)
        {
          // Process 5: Repeat the second process and so on. until there are no
          // more points in the seed
          if (!visited[nextNode] || clusters[nextNode].empty())
          {
            visited[nextNode] = true;
            seeds.push(nextNode);
            if (clusters[nextNode].size() < 1)
            {
              clusters[nextNode].push_back(clusterID);
            }
          }
        }
      }

      // Process 6: Repeat the first process and so on. until all clusters are
      // filled
      clusterID++;
    }
  }
  corePointStatus.assign(corePoint.begin(), corePoint.end());
}

// Function to calculate the Rand index
double calculateRandIndex(const vector<int> &clusterLabels,
                          const vector<int> &trueLabels)
{
  int n = clusterLabels.size();
  int a = 0, b = 0, c = 0, d = 0;

  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      if (clusterLabels[i] == clusterLabels[j] &&
          trueLabels[i] == trueLabels[j])
      {
        a++;
      }
      else if (clusterLabels[i] != clusterLabels[j] &&
               trueLabels[i] != trueLabels[j])
      {
        b++;
      }
      else if (clusterLabels[i] == clusterLabels[j] &&
               trueLabels[i] != trueLabels[j])
      {
        c++;
      }
      else if (clusterLabels[i] != clusterLabels[j] &&
               trueLabels[i] == trueLabels[j])
      {
        d++;
      }
    }
  }

  double randIndex = static_cast<double>(a + b) / (a + b + c + d);
  return randIndex;
}

int main()
{

  string filename =
      "/home/irawan/Documents/WUT/Data Mining/update_project/small_dataset.csv";
  DataFromCSV dataFromCSV = readCSVFile(filename);

  // parameter for measure the computation time
  time_t start = clock();

  int k = 3;
  Graph graph;
  buildGraph(dataFromCSV, k, graph);
  vector<vector<int>> clusters;
  assignClusters(graph, clusters);

  time_t end = clock();
  double elapsed = double(end - start) / CLOCKS_PER_SEC;

  vector<int> pred_clust(dataFromCSV.size(), 0);
  vector<int> true_clust;

  ofstream results("results.txt");

  if (results.is_open())
  {
    results << "cluster" << endl;
    for (const auto &cluster : clusters)
    {
      if (cluster.size() == 0)
      {
        results << -1 << endl;
      }
      for (const auto &data : cluster)
      {
        results << data << endl;
      }
    }

    results.close();
    cout << "File Results.txt created successfully." << endl;
  }
  else
  {
    cout << "Failed to create file OUT.txt" << endl;
  }

  string clusterfile = "/home/irawan/Documents/WUT/Data\ Mining/update_project/"
                       "output/results.txt";
  DataFromCSV dataClusters = readCSVFile(clusterfile);

  for (int i = 0; i < dataClusters.size(); i++)
  {
    for (int j = 0; j < dataClusters[0].size(); j++)
    {
      true_clust.push_back(dataClusters[i][j]);
    }
  }

  // Using RAND index for the clustering evaluation
  double rand_index = calculateRandIndex(pred_clust, true_clust);
  cout << "RAND index value : " << rand_index << endl;
  cout << "Computation Time : " << elapsed << " seconds" << endl;

  ofstream outputFile1("OUT_NBC-EuclideanDistance_fname_D2_R10000_k3.txt");
  ofstream outputFile2("STAT_NBC-EuclideanDistance_fname_D2_R10000_k3.txt");
  ofstream outputFile3("k+NN_NBC-EuclideanDistance_fname_D2_R10000_k3.txt");

  if (outputFile1.is_open())
  {
    outputFile1 << "Data dimension : " << dataFromCSV.size() << " x "
                << dataFromCSV[0].size() << ", ";
    for (int i = 0; i < dataFromCSV.size(); i++)
    {
      for (int j = 0; j < dataFromCSV[0].size(); j++)
      {
        vector<double> row = extractColumn(dataFromCSV, j);
        outputFile1 << "Data ID : " << j << ", Feature " << i << " : " << row[j]
                    << " ,";
      }
      outputFile1 << endl;
    }

    for (int j = 0; j < database.size(); j++)
    {
      int source = get<0>(database[j]);
      int dest = get<1>(database[j]);
      double distance = get<2>(database[j]);
      outputFile1 << "Data ID source : " << source
                  << ", Data ID destination : " << dest
                  << ", distance : " << distance << endl;
    }

    for (int i = 0; i < dataFromCSV[0].size(); i++)
    {
      outputFile1 << "Data ID :" << i
                  << ", Status Core point : " << corePointStatus[i] << endl;
    }

    int dataId = 0;
    for (const auto &cluster : clusters)
    {
      if (cluster.size() == 0)
      {
        dataId++;
        outputFile1 << "Data ID " << dataId << " : cluster " << -1 << endl;
      }
      for (const auto &data : cluster)
      {
        dataId++;
        outputFile1 << "Data ID " << dataId << " : cluster " << data << endl;
      }
    }

    outputFile1.close();
    cout << "File OUT.txt created successfully." << endl;
  }
  else
  {
    cout << "Failed to create file OUT.txt" << endl;
  }

  if (outputFile2.is_open())
  {
    outputFile1 << "Data dimension : " << dataFromCSV.size() << " x "
                << dataFromCSV[0].size() << ", ";
    outputFile2 << "k parameter : " << k << endl;
    outputFile2 << "Total runtime : " << elapsed << " seconds" << endl;
    outputFile2 << "RAND index value : " << rand_index << endl;

    outputFile2.close();
    cout << "File STAT.txt created successfully." << endl;
  }
  else
  {
    cout << "Failed to create file STAT.txt" << endl;
  }

  if (outputFile3.is_open())
  {
    for (int i = 0; i < graph.size(); i++)
    {
      int size_knn = graph[i].size();

      outputFile3 << "Data ID " << i << " KNN : ";
      for (int j = 0; j < graph[i].size(); j++)
      {
        outputFile3 << graph[i][j] << ", ";
      }
      outputFile3 << endl;

      count_rknn = 0;
      for (const auto &inner_vec : graph)
      {
        for (const auto &value : inner_vec)
        {
          if (value == i)
          {
            count_rknn++;
          }
        }
      }
      outputFile3 << "RkNN of Data ID  " << i << " : " << count_rknn << endl;

      double ndf;
      ndf = count_rknn / size_knn;
      outputFile3 << "NDF of Data ID " << i << " : " << ndf << endl;
    }

    outputFile3.close();
    cout << "File k+NN.txt created successfully." << endl;
  }
  else
  {
    cout << "Failed to create file  k+NN.txt" << endl;
  }

  return 0;
}