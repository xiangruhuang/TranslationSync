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

int main(int argc, char** argv){
    if (argc <= 12){
        cerr << "./TranSync [n] [edge_density] [bias] [inc] [noise_type] [noise_ratio] [max_iter] [decay] [a/mean] [b/std_dev] [sigma] [output_name] (graph_file)" << endl;
        return 0;
    }
    auto seed = time(NULL);
    srand(seed);
    cerr << "random seed=" << seed << endl;
    Params* params = new Params(argc, argv);
    cerr << params->n << endl;
    Graph graph(params);
    //graph.dump("graph4.meta");
    //return 0;
    Truncated_L2 TL(params);

    CoordinateDescent CD(params);
  
    string output_name(argv[12]);
    vector<double> result_TL2;
    vector<double> result_CD;
    vector<double> time_TL2;
    vector<double> time_CD;

    int T = 100;
    int P = 100;
    for (int t = 0; t < T; t++){
        result_TL2.push_back(t);
        result_CD.push_back(t);
        time_TL2.push_back(t);
        time_CD.push_back(t);
    }
    int nr = round((params->noise_ratio)*P);
    string nr_str = to_string(nr);
    #pragma omp parallel for
    for (int t = 0; t < T; t++){
        string str = to_string(t);
        Graph graph_t = graph.copy();
        graph_t.resample();

        double time2 = -omp_get_wtime();
        double r2 = CD.solve(graph_t, output_name+"/ratio"+nr_str+"_"+str+".CD");
        time2 += omp_get_wtime();
        
        double time1 = -omp_get_wtime();
        double r1 = TL.solve(graph_t, output_name+"/ratio"+nr_str+"_"+str+".TL2");
        time1 += omp_get_wtime();

        result_TL2[t] = r1;
        time_TL2[t] = time1;
        result_CD[t] = r2;
        time_CD[t] = time2;
    }
    ofstream fout(output_name+"/ratio"+nr_str+"_summary");
    for (int t = 0; t < T; t++){
        fout << result_TL2[t] << " " << time_TL2[t] << " " << result_CD[t] << " " << time_CD[t] << endl;
    }
    //double zp_TL2 = zero_prob(result_TL2);
    //double zp_CD = zero_prob(result_CD);
    cerr << "noise_ratio=" << params->noise_ratio << ", TL2=(" << min(result_TL2) << "," << median(result_TL2) << "," << max(result_TL2) << ")";
    cerr << ", CD=(" << min(result_CD) << "," << median(result_CD) << "," << max(result_CD) << ")";
    cerr << ", time mean: TL2=" << mean(time_TL2) << ", CD=" << mean(time_CD) << endl;
    fout.close();
    
    return 0;
}
