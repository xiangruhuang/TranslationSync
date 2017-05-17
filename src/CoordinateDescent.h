#ifndef CD_H
#define CD_H

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

using namespace std;

class CoordinateDescent : public Solver{
    public:
    int max_iter;
    CoordinateDescent(Params* params){
        this->max_iter = params->max_iter;
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

        int iter = 0;
        double min_loss = 1e100;

        int* d = new int[n];
        for (int i = 0; i < n; i++){
            d[i] = v[i].size();
        }
        
        double elapsed_time = 0.0;
        double last_update = omp_get_wtime();
        while (iter++ < max_iter){
            double up = 0.0, down = 0.0;
            for (int i = 0; i < n; i++){
                vector<double> c;
                c.clear();
                for (vector<pair<double, int>>::const_iterator it_i = v[i].cbegin(); it_i != v[i].cend(); it_i++){
                    int j = it_i->second;
                    double t_ij = it_i->first;
                    c.push_back(x[j] + t_ij);
                }
                int middle_point = c.size() / 2;
                nth_element(c.begin(), c.begin()+middle_point, c.end());
                new_x[i] = c[middle_point];
                if (c.size() % 2 == 0){
                    nth_element(c.begin(), c.begin()+middle_point-1, c.end());
                    new_x[i] = 0.5*(new_x[i] + c[middle_point-1]);
                }
                up += new_x[i]*sqrt(d[i]);
                down += sqrt(d[i]); 
            }
            double shift = up/down;

            double delta_x = 0.0;
            for (int i = 0; i < n; i++){
                new_x[i] -= shift;
                delta_x += fabs(x[i] - new_x[i]);
                x[i] = new_x[i];
            }
            double loss = linf_loss(n, x, graph.x);

            if (loss < min_loss){
                min_loss = loss;
                last_update = omp_get_wtime();
            }

            elapsed_time = (omp_get_wtime() - start_time);
            fout << "iter=" << iter;
            fout << ", linf_loss=" << loss;
            fout << ", l1_loss=" << l1_loss(n, x, graph.x);
            fout << ", min_loss=" << min_loss;
            fout << ", elapsed_time=" << elapsed_time;
            fout << endl;
            if (fabs(min_loss) < 1e-6 || fabs(delta_x) < 1e-10 || (omp_get_wtime() - last_update > 10.0)){
                break;
            }
        }
        fout.close();
        delete d;
        delete x;
        delete new_x;
        return min_loss;
    }
};

#endif
