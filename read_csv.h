#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

class ReadCSV{
    vector< vector<double> > content;
    vector<double> x, y;
    string line;

    public :
        vector< vector<double> > read_file(const string& fname){
            fstream file (fname, ios::in);
            if(file.is_open())
            {
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
                content.push_back(x);
                content.push_back(y);
            } else {
                cout<<"Could not open the file \n";
            }

            return content;
        }

};