n = 2000;
sigma = 1e-3;
error_rate = 0.1; p = 1-error_rate;
sigma_t = sqrt(p*sigma^2 + (1.0-p)/3.0);
num_experiments = 100;
c = 2; % success prob. >= 1-n^{-c}
type = 1;

fprintf('generating graph: type %d...', type);
B = make_graph(type);
fprintf('done.\n');
L = B'*B;
% LL = full(L);
% Ld = pinv(LL);
% max_eig = max(eig(Ld));
% max_diag = max(diag(Ld));
% delta = sqrt((c+1)*log(n)*max_diag/2.0)*sigma_t;
% upper_delta = sqrt((c+1)*log(n)*max_eig/2.0)*sigma_t;

num_edges = size(B, 1);

results_us = [];
results_interior = [];
for i = 1:num_experiments
    
    t = make_observations(num_edges, sigma, p);
    
%     y = B'*t;
%     x = theirCG(L, y, zeros(n, 1));
    x = CD(B, L, t, randn(n, 1));
    result = norm(x, 'inf');
    results_us = [results_us result];
end
fprintf('mean_delta=%f', mean(results_us));
% fprintf(', delta=%f (%f)', delta, mean(inf_norms < delta));
% fprintf(', upper_delta=%f (%f)', upper_delta, mean(inf_norms < upper_delta));
fprintf('\n');