#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cassert>
#include "graph.h"
#include <algorithm>
#include <time.h>
#include <omp.h>
#include "Truncated_L2.h"
#include "CoordinateDescent.h"

using namespace std;

struct point{
    double x, y;
};

int main(){
    ifstream fin("sample.txt");

    Graph graph;
    int n;
    int m;
    fin >> n >> m;
    graph.n = n;
    graph.m = m;
    graph.x = new double[n];
    graph.x0 = new double[n];
    graph.adj = new vector<pair<double, int>>[n];

    cerr << "0" << endl;
    point* p = new point[n];
    for (int i = 0; i < n; i++){
        fin >> p[i].x >> p[i].y;
        graph.x[i] = p[i].y;
        graph.x0[i] = 0.0;
    }
    cerr << "1" << endl;
    point* edge = new point[m];
    int* d = new int[n];
    for (int e = 0; e < m; e++){
        int i, j;
        fin >> i >> j >> edge[e].x >> edge[e].y;
        i--;
        j--;
        d[i]++;
        d[j]++;
        graph.adj[i].push_back(make_pair(edge[e].y,j));
        graph.adj[j].push_back(make_pair(-edge[e].y,i));
    }
    cerr << "2" << endl;
    for (int i = 0; i < n; i++){
        if (d[i] == 0){
            cerr << "i=" << i << endl;
        }
    }
    cerr << "hey" << endl;
    graph.normalize(n, d, graph.x);
    CoordinateDescent CD(10000);
    double ans_CD = CD.solve(graph, "dataI.CDY");
    cout << "CD=" << ans_CD << endl;
    Truncated_L2 TL2(10000, 0.8);
    double ans_TL2 = TL2.solve(graph, "dataI.TL2Y");
    cout << "TL2=" << ans_TL2 << endl;

}
