n = 2000;
sigma = 1e-3;
error_rate = 0.1; p = 1-error_rate;
sigma_t = sqrt(p*sigma^2 + (1.0-p)/3.0);
num_experiments = 100;
c = 2; % success prob. >= 1-n^{-c}
type = 1;
% 
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

% B = make_graph(type);
% L = B'*B;
num_edges = size(B, 1);

result_CG = [];
result_CD = [];
observations = [];
infs_CG = [];
infs_CD = [];
times_CG = [];
times_CD = [];
for i = 1:num_experiments
    
    t = make_observations(num_edges, sigma, p);
    observations = [observations t];
    x0 = randn(n, 1);

%     [x1, infs, times] = theirCG(L, B'*t, x0);
%     infs_CG = [infs_CG infs];
%     times_CG = [times_CG times];
%     result_CG = [result_CG norm(x1, 'inf')];

%     [x2, infs, times] = CD(B, L, t, x0);
%     infs_CD = [infs_CD infs];
%     times_CD = [times_CD times];
%     result_CD = [result_CD norm(x2, 'inf')];
    
    x3 = IPM(B, t, x0);

%     fprintf('experiment %d: %f (CG), %f (CD)\n', i, norm(x1, 'inf'), norm(x2, 'inf'));
%     save('type1.mat');
end
% save('type1.mat');
% fprintf('mean_delta=%f', mean(results_us));
% fprintf(', delta=%f (%f)', delta, mean(inf_norms < delta));
% fprintf(', upper_delta=%f (%f)', upper_delta, mean(inf_norms < upper_delta));
% fprintf('\n');