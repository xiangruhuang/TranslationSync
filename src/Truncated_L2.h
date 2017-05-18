#ifndef TRUNCATED_L2_H
#define TRUNCATED_L2_H

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
#include "util.h"
#include <string>
#include <random>

using namespace std;

class Truncated_L2 : public Solver{
    public:
    int max_iter;
    double decay;

    Truncated_L2(Params* params){
        this->max_iter = params->max_iter;
        this->decay = params->decay;
    }
    Truncated_L2(int max_iter, double decay){
        this->max_iter = max_iter;
        this->decay = decay;
    }

    inline double solve(Graph& graph, string filename){
        ofstream fout(filename);
        double start_time = omp_get_wtime();
        int n = graph.n;
        int m = graph.m;
        vector<pair<double, int>>* v = graph.adj;
        double* x = new double[n];
        for (int i = 0; i < n; i++){
            x[i] = graph.x0[i];
        }
        double* new_x = new double[n];
        double* old_x = new double[n];
        int iter = 0;
        double min_loss = 1e100;
        int* index = new int[n];
        for (int i = 0; i < n; i++){
            index[i] = i;
        }
        bool disconnected = false;
        double threshold = 1e100;

        int* d = new int[n];
        memset(d, 0, sizeof(int)*n);

        double up = 0.0;
        double down = 0.0;
        for (int i = 0; i < n; i++){
            for (vector<pair<double, int>>::const_iterator it_i = v[i].cbegin(); it_i != v[i].cend(); it_i++){
                int j = it_i->second;
                d[i]++;
                d[j]++;
            }
        }
        int edges_remain = 0;
        for (int i = 0; i < n; i++){
            down += sqrt(d[i]);
            edges_remain += d[i];
        }
        while (iter++ < max_iter && (edges_remain >= 2*(n-1)) && (!disconnected)){
            double delta_x = 1e100;
            int inner = 0;
            double max_diff = 0.0;
            double sum_diff = 0.0;
            vector<int> remove_list;
            for (int i = 0; i < n; i++){
                old_x[i] = x[i];
            }
            int count_down = 0;
            while (delta_x > 1e-2 && (!disconnected) && (edges_remain >= 2*(n-1))){
                inner++;
                max_diff = 0.0;
                up = 0.0;
                for (int i = 0; i < n; i++){
                    int d_i = 0.0;
                    double weighted_sum = 0.0;
                    remove_list.clear();
                    for (vector<pair<double, int>>::const_iterator it_i = v[i].begin(); it_i != v[i].end(); it_i++){
                        int j = it_i->second;
                        double t_ij = it_i->first;
                        //cerr << "i=" << i << ", j=" << j << ", t_ij=" << t_ij << endl;
                        double diff_i = fabs(t_ij - (x[i] - x[j]));
                        if ( diff_i < threshold ){
                            weighted_sum += x[j] + t_ij;
                            d_i++;
                            if (diff_i > max_diff){
                                max_diff = diff_i;
                            }
                        } else {
                            remove_list.push_back(it_i-v[i].cbegin());
                        }
                    }
                    for (int k = remove_list.size()-1; k >= 0; k--){
                        int to_remove = remove_list[k];
                        v[i][to_remove] = v[i][v[i].size()-1];
                        v[i].pop_back();
			            edges_remain--;
                    }
                    if (d_i == 0){
                        disconnected = true;
                        break;
                    }
                    //cerr << "weighted_sum=" << weighted_sum << ", d_i=" << d_i << endl;
                    new_x[i] = weighted_sum / d_i;
                    up += new_x[i] * sqrt(d[i]);
                }
                double shift = up/down;
                delta_x = 0.0;
                for (int i = 0; i < n; i++){
                    new_x[i] -= shift;
                    delta_x += fabs(new_x[i] - x[i]);
                    //cerr << "new_x=" << new_x[i] << ", x=" << x[i] << ", gt=" << graph.x[i] << endl;
                }
                double* temp = new_x; new_x = x; x = temp;
                //cerr << "delta_x=" << delta_x << endl;
            }

            //maintain and output
            //if (iter == 1){
            //    threshold = max_diff;
            //}
            if (threshold > max_diff){
                threshold = max_diff;
            }
            threshold = threshold * decay;

            double loss = linf_loss(n, x, graph.x);
            if (loss < min_loss){
                min_loss = loss;
            }
            double outer_delta_x = 0.0;
            for (int i = 0; i < n; i++){
                outer_delta_x += (x[i] - old_x[i]);
            }
            fout << "iter=" << iter;
            fout << ", #inner=" << inner;
            fout << ", linf_loss=" << loss;
            fout << ", l1_loss=" << l1_loss(n, x, graph.x);
            fout << ", min_loss=" << min_loss;
            fout << ", threshold=" << threshold;
            fout << ", outer_delta_x=" << outer_delta_x;
            fout << ", elapsed_time=" << (omp_get_wtime() - start_time);
            fout << endl;
            if (fabs(outer_delta_x) < 1e-5){
                count_down++;
            } else {
                count_down = 0;
            }
            if (count_down >= 3){
                break;
            }
        }
        fout.close();
        delete d;
        delete index;
        delete x;
        delete new_x;
        return min_loss;
    }

};
#endif
