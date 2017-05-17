#ifndef UTIL_H
#define UTIL_H

double l1_loss(int n, double* x, double* gt){
    //compute l1 loss to ground truth
    double loss = 0.0;
    for (int i = 0; i < n; i++){
        loss += fabs(gt[i] - x[i]);
    }
    return loss;
}

double linf_loss(int n, double* x, double* gt){
    //compute l1 loss to ground truth
    double loss = 0.0;
    for (int i = 0; i < n; i++){
        loss = max(loss, fabs(gt[i] - x[i]));
    }
    return loss;
}


double mean(vector<double> v){
    double sum = 0.0;
    for (auto it = v.begin(); it != v.end(); it++){
        sum += *it;
    }
    if (v.size() > 0){
        return sum/v.size();
    } else {
        return 0.0;
    }
}

double median(vector<double> v){
    vector<double> dup = v;
    sort(dup.begin(), dup.end());
    if (dup.size() % 2 == 1){
        return dup[dup.size() / 2];
    } else {
        return (dup[dup.size() / 2] + dup[(dup.size() / 2)-1])/2.0;
    }
}

double max(vector<double> v){
    double ans = -1e100;
    for (auto it = v.begin(); it != v.end(); it++){
        if (*it > ans){
            ans = *it;
        }
    }
    return ans;
}

double min(vector<double> v){
    double ans = 1e100;
    for (auto it = v.begin(); it != v.end(); it++){
        if (*it < ans){
            ans = *it;
        }
    }
    return ans;
}

double zero_prob(vector<double> v){
    double sum = 0.0;
    for (auto it = v.begin(); it != v.end(); it++){
        if (*it < 1e-5){
            sum += 1.0;
        }
    }
    if (v.size() > 0){
        return sum/v.size();
    } else {
        return 0.0;
    }
}

#endif
