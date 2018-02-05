#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

int main(){
    ifstream finx("TL2.x");
    int n = 6001;
    double* x = new double[n];
    for (int i = 0; i < n; i++){
        finx >> x[i];
    }
    finx.close();
    
    ifstream finy("TL2.y");
    double* y = new double[n];
    for (int i = 0; i < n; i++){
        finy >> y[i];
    }
    finy.close();
   
    int m = 140381;
    ifstream fin("sample.txt");

    fin >> n >> m;
    for (int i = 0; i < n; i++){
        double t;
        fin >> t;
        fin >> t;
    }

    vector<double> dists;
    for (int e = 0; e < m; e++){
        int i, j;
        double t_ij_x;
        double t_ij_y;
        fin >> i >> j >> t_ij_x >> t_ij_y;
        i--; j--;
        dists.push_back(fabs(t_ij_x-(x[i]-x[j])));
        dists.push_back(fabs(t_ij_y-(y[i]-y[j])));
    }
    fin.close();
    
    sort(dists.begin(), dists.end());
    ofstream fout("TL2.dist");
    for (int e = 0; e < m*2; e++){
        fout << dists[e] << endl;
    }
    fout.close();

    delete x;
}
