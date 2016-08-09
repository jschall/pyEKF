#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>

#include "PLEKF.h"

class Copter {
public:
    float pos[3] {0.0f, 2.0f, -5.0f};
    float vel[3] {0.0f, 0.0f, 0.0f};
};

using namespace std;
int main(int argc, char* argv[])
{
//     if (argc < 2) {
//         cerr << argv[0] << " <csvfile>" << endl;
//         return 1;
//     }
    
//     char* csvfn = argv[1];
//     ifstream file(csvfn);
// 
//     vector< vector<string> > lines;
//     while (file.good()) {
//         string s;
//         getline(file,s);
//         stringstream ss(s);
//         lines.push_back(vector<string>());
//         string item;
//         while (getline(ss, item, ',')) {
//             lines.back().push_back(item);
//         }
//     }
//     
//     int i = 0;
//     for(vector< vector<string> >::iterator line = lines.begin(); line != lines.end(); ++line) {
//         cout << i << endl;
//         for (vector<string>::iterator item = line->begin(); item != line->end(); ++item) {
//             cout << "|-> " << *item << endl;
//         }
//         i++;
//     }
    
    float copterPos[3];
    float copterVel[3];
    float dt = 1.0f/60.0f;
    
    Copter copter;
    PLEKF pl_ekf;
    
    pl_ekf.initialize();
}