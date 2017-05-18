#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <stdlib.h>
#include <string>
#include <algorithm>

using namespace std;

class Params{
    public:
    //Graph specific parameters:
    int n; //#nodes
    double p; //Probability of each edge
    double noise_ratio; // probability that a sample has noise
    int noise_type; 
    // possible noise types:
    // 0: gaussian with mean 0 and std_dev=0.001
    // 1: uncentered gaussian with mean 1 and std_dev=1
    // 2: uncentered uniform [-1, 3]
    bool load_graph = false; // if true, load graph from file
    char* graph_file_name;

    double bias;
    double inc;
    double* s;
    //weight of each node, s[i] = bias + inc*i/(n-1)
    //edge (i,j) exists w.p. p*s[i]*s[j]
    
    //sample specific parameters:
    double a, b;
    
    default_random_engine generator;

    //solver specific parameters:
    int max_iter = 1000;
    double decay = 0.9;

    Params(){

    }

    Params(int argc, char** argv){
        n = atoi(argv[1]);
        p = atof(argv[2]);
        bias = atof(argv[3]);
        inc = atof(argv[4]);
        noise_type = atoi(argv[5]);
        noise_ratio = atof(argv[6]);
        s = new double[n];
        for (int i = 0; i < n; i++){
            s[i] = bias + inc*i/(n-1);
        }

        max_iter = atoi(argv[7]);
        decay = atof(argv[8]);
        
        a = atof(argv[9]);
        b = atof(argv[10]);

        if (argc > 12){
            load_graph = true;
            graph_file_name = argv[12];
        } else {
            load_graph = false;
        }
        
        

        cerr << "Initiating Parameters: ";
        cerr << "n=" << n;
        cerr << ", edge_density=" << p;
        cerr << ", bias=" << bias;
        cerr << ", incremental=" << inc;
        cerr << ", noise_type=" << noise_type;
        cerr << ", noise_ratio=" << noise_ratio;
        cerr << ", noise_parameters=(" << a << "," << b << ")";
        cerr << endl;
    }

    ~Params(){
        delete s;
    }

    double noise(){
        normal_distribution<double> small_noise(0.0, 0.01);
        normal_distribution<double> gaussian(a, b);
        uniform_real_distribution<double> uniform(a, b);
        uniform_real_distribution<double> dice(0.0, 1.0);
        double noise = small_noise(generator);
        if (dice(generator) <= noise_ratio){
            if (noise_type == 1){
                //gaussian
                noise = gaussian(generator);
            }
            if (noise_type == 2){
                //uniform
                noise = uniform(generator);
            }
        }
        return noise;    
    }

    double ground_truth(){
        normal_distribution<double> ground_truth(0.0, 1.0); // should not matter, can be just zero
        return ground_truth(generator); 
    }
};

class Graph{
    public:
    int n; // #nodes
    int m; // #edges
    double* x;
    double* x0;
    vector<pair<double, int>>* adj; // #edges in terms of adjacency matrix
    vector<pair<int, int>> edges;
    Params* params;
    bool* visited; // used to check connectivity
    int num_visited;
    Graph(){
    }

    Graph(Params* _params){
        params = _params;
        n = params->n;
        //generate ground truth x
        x = new double[n];
        x0 = new double[n];
        for (int i = 0; i < n; i++){
            x[i] = params->ground_truth();
        }
        adj = new vector<pair<double, int>>[n];

        //generate samples, possibly with noise
        if (params->load_graph){
            edges.clear();
            ifstream fin(params->graph_file_name);
            fin >> n >> m;
            for (int e = 0; e < m; e++){
                int i, j;
                fin >> i >> j;
                assert(i < j);
                edges.push_back(make_pair(i, j));
                adj[i].push_back(make_pair(0.0, j));
                adj[j].push_back(make_pair(0.0, i));
            }
            fin.close();
            cerr << "Done Loading Graph From " << params->graph_file_name;
            cerr << "#nodes=" << n;
            cerr << ", #edges=" << m;
            cerr << endl;
        } else {
            m = 0;
            double p = params->p;
            double* s = params->s;
            edges.clear();
            for (int i = 0; i < n; i++){
                double s_i = s[i];
                for (int j = i+1; j < n; j++){
                    double s_j = s[j];
                    double dice = rand()*1.0 / RAND_MAX;
                    if (dice <= p * s_i * s_j){
                        edges.push_back(make_pair(i, j));
                        double t_ij = x[i] - x[j] + params->noise();
                        adj[i].push_back(make_pair(t_ij, j));
                        adj[j].push_back(make_pair(-t_ij, i));
                        m++;
                    }
                }
            }
            visited = new bool[n];
            num_visited = 0;
            memset(visited, false, sizeof(bool)*n);
            dfs(0);
            int num_extra_edges = 0;
            for (int i = 1; i < n; i++){
                if (!visited[i]){
                    int dice = rand() % num_visited;
                    int count = 0;
                    for (int j = 0; j < n; j++){
                        if (visited[j]){
                            if (count == dice){
                                //add edge to (i, j) and then dfs(j)
                                double t_ij = x[i] - x[j] + params->noise();
                                adj[i].push_back(make_pair(t_ij, j));
                                adj[j].push_back(make_pair(-t_ij, i));
                                edges.push_back(make_pair(min(i, j), max(i, j)));
                                dfs(i);
                                m++;
                                num_extra_edges++;
                                break;
                            } else {
                                count++;
                            }
                        }
                    }
                }
            }
            delete visited;
            cerr << "Done Generating Graph: ";
            cerr << "#nodes=" << n;
            cerr << ", #edges=" << m << " (#extras=" << num_extra_edges << ")";
            cerr << endl;
        }
        //ground truth normalization/shift
        int* d = new int[n];
        memset(d, 0, sizeof(int)*n);
        for (auto e = edges.begin(); e != edges.end(); e++){
            int i = e->first, j = e->second;
            d[i]++;
            d[j]++;
        }
        normalize(n, d, x);
        normalize(n, d, x0);
    }

    ~Graph(){
        for (int i = 0; i < n; i++){
            adj[i].clear();
        }
        delete x0;
        delete x;
    }

    inline void normalize(int n, int* d, double* z){
        double up = 0.0, down = 0.0;
        for (int i = 0; i < n; i++){
            up += sqrt(d[i])*z[i];
            down += sqrt(d[i]);
        }
        double shift = up/down;
        for (int i = 0; i < n; i++){
            z[i] -= shift;
        }
    }

    inline Graph copy(){
        Graph g;
        g.n = this->n;
        g.m = this->m;
        g.adj = new vector<pair<double, int>>[g.n];
        g.edges = this->edges;
        g.x = new double[g.n];
        g.x0 = new double[g.n];
        for (int i = 0; i < g.n; i++){
            g.x[i] = x[i];
            g.x0[i] = x0[i];
        }
        g.params = this->params;
        return g;
    }

    void resample(){
        for (int i = 0; i < n; i++){
            adj[i].clear();
            x0[i] = (rand()*1.0 / RAND_MAX) * 100.0;
        }
        for (auto e = edges.begin(); e != edges.end(); e++){
            int i = e->first, j = e->second;
            double t_ij = x[i] - x[j] + params->noise();
            adj[i].push_back(make_pair(t_ij, j));
            adj[j].push_back(make_pair(-t_ij, i));
        }
    }

    void dump(char* filename){
        ofstream fout(filename);
        fout << n << " " << m << endl;
        int count = 0;
        for (int e = 0; e < edges.size(); e++){
            int i = edges[e].first, j = edges[e].second;
            assert(i < j);
            fout << i << " " << j << endl;
        }
        fout.close();
    }

    private:
    void dfs(int i){
        num_visited++;
        visited[i] = true;
        for (auto it = adj[i].begin(); it != adj[i].end(); it++){
            int j = it->second;
            if (!visited[j]){
                dfs(j);
            }
        }
    }
};

class Solver{
    virtual inline double solve(Graph& graph, string output_name){
    };
};

#endif
