n = 100;
edge_density = 0.3;
sigma = 1e-3;
error_rate = 0.1; p = 1-error_rate;
sigma_t = sqrt(p*sigma^2 + (1.0-p)/3.0);
num_experiments = 1000;
c = 2; % success prob. >= 1-n^{-c}

B = make_Gnp(n, edge_density);
L = B'*B;
Ld = pinv(L);
max_eig = max(eig(Ld));
max_diag = -1000000;
for k = 1:n
   if max_diag < Ld(k, k)
      max_diag = Ld(k, k) 
   end
end
num_edges = size(B, 1);
delta = sqrt((c+1)*log(n)*max_diag/2.0)*sigma_t;
upper_delta = sqrt((c+1)*log(n)*max_eig/2.0)*sigma_t;

inf_norms = [];
for i = 1:num_experiments
    
    t = make_observations(num_edges, sigma, p);

    x = theirCG(L, B'*t, zeros(n, 1));
    result = norm(x, 'inf');
    inf_norms = [inf_norms result];
end
fprintf('upper_delta=%f, delta=%f, success_prob_delta=%f, success_prob_delta=%f\n', upper_delta, delta, mean(inf_norms <= delta), mean(inf_norms <= upper_delta));
