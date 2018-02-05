#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;

struct point{
    double f, k1, k2;
    vector<double> r;
    vector<double> t;
};

struct comparison{
    int i, j;
    vector<double> t_ij;
};

inline vector<double> transpose(int n, int m, vector<double> a){
    vector<double> ans;
    for (int j = 0; j < m; j++){
        for (int i = 0; i < n; i++){
            ans.push_back(a[i*m+j]);
        }
    }
    return ans;
}

//compute <a, b>, assume a.size() > b.size()
inline vector<double> multiply(vector<double> a, vector<double> b){
    vector<double> ans;
    int n = a.size(), m = b.size();
    assert(n % m == 0);
    auto it_a = a.cbegin();
    for (int i = 0; i < n; i+=m){
        double dot_product = 0.0;
        auto it_b = b.cbegin();
        for (auto it_b = b.cbegin(); it_b != b.cend(); it_b++, it_a++){
            dot_product += *it_b*(*it_a);
        }
        ans.push_back(dot_product);
    }
    return ans;
}

//compute a - b, assume a.size() > b.size()
inline vector<double> subtract(vector<double>& a, vector<double>& b){
    vector<double> ans;
    int n = a.size(), m = b.size();
    assert(n == m);
    auto it_a = a.cbegin();
    for (auto it_b = b.cbegin(); it_b != b.cend(); it_b++, it_a++){
        ans.push_back((*it_a)-(*it_b));
    }
    return ans;
}

void get_cameras(char* filename, vector<point>& ans){
    ifstream fin(filename);
    int n, m;
    fin >> n >> m;
    for (int i = 0; i < n; i++){
        point camera;
        fin >> camera.f >> camera.k1 >> camera.k2;
        for (int j = 0; j < 9; j++){
            double r_j;
            fin >> r_j;
            camera.r.push_back(r_j);
        }
        for (int j = 0; j < 3; j++){
            double t_j;
            fin >> t_j;
            camera.t.push_back(t_j);
        }
        ans.push_back(camera);
    }
    fin.close();
}

void get_relative_translations(char* filename, vector<comparison>& rel){
    ifstream fin(filename);
    int n, m;
    fin >> n >> m;
    for (int e = 0; e < m; e++){
        comparison comp;
        fin >> comp.i >> comp.j;
        for (int k = 0; k < 9; k++){
            double temp;
            fin >> temp;
        }
        for (int k = 0; k < 3; k++){
            double t_ij_k;
            fin >> t_ij_k;
            comp.t_ij.push_back(t_ij_k);
        }
        for (int k = 0; k < 4; k++){
            double temp;
            fin >> temp;
        }
        rel.push_back(comp);
    }
    fin.close();
}

inline void print_vec(vector<double> v){
    for (auto it = v.cbegin(); it != v.cend(); it++){
        cout << *it << " ";
    }
    cout << endl;
}

int main(int argc, char** argv){
    if (argc <= 1){
        cerr << "./preprocess [camera] [relative translation]" << endl;
        return -1;
    }
    vector<point> camera;
    get_cameras(argv[1], camera);
    int n = camera.size();
    vector<comparison> rel;
    get_relative_translations(argv[2], rel);
    double loss = 0.0;
    double norm = 0.0;
    for (auto it = rel.cbegin(); it != rel.cend(); it++){
        int i = it->i, j = it->j;
        vector<double> t_ij = it->t_ij;
        vector<double>& r_i = camera[i].r;
        vector<double>& t_i = camera[i].t;
        vector<double>& t_j = camera[j].t;
        vector<double> rhs = multiply(transpose(3,3, r_i), subtract(t_j, t_i));
        double norm_l = 0.0, norm_r = 0.0;
        double norm_rhs = 0.0;
        for (int k = 0; k < 3; k++){
            norm_r += rhs[k]*t_ij[k];
            norm_l += t_ij[k]*t_ij[k];
            norm_rhs += rhs[k]*rhs[k];
        }
        norm += sqrt(norm_rhs);
        double lambda = 1.0;
        if (fabs(norm_l) > 1e-20){
            lambda = norm_r/norm_l;
            //print_vec(t_ij);
            //print_vec(rhs);
        } else {
            assert(false);
        }
        double dist = 0.0;
        for (int k = 0; k < 3; k++){
            double diff = t_ij[k]*lambda-rhs[k];
            dist += diff*diff;
        }
        loss += sqrt(dist);
    }
    cout << "loss=" << loss;
    cout << ", avg_loss=" << loss/rel.size();
    cout << ", norm=" << norm;
    cout << ", avg_norm=" << norm/rel.size();
    cout << endl;
    return 0;
}
