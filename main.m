% sigma = 0.00;
% % error_rate = 0.9; p = 1-error_rate;
% % sigma_t = sqrt(p*sigma^2 + (1.0-p)/3.0);
% num_experiments = 1;
% c = 2; % success prob. >= 1-n^{-c}

initialize;

sample_setting.sigma = 0.0;
verbose = false;

num_experiments = 100;
for graph_type = 1:4
    %load graph, process
    graph_name = strcat(['Graph_',num2str(graph_type)]);
    load(strcat([graph_name, '.mat']));
    
    params = graph_params{graph.type};
    
    experiment.graph_name = graph_name;
    experiment.graph_params = params;
    dir = strcat(['./results_symmetric/', graph_name]);
    
    %perform experiments
    for error_rate = 1.0:-0.01:0.0
        fprintf('Error Rate=%f\n', error_rate);
        sample_setting.error_rate = error_rate;
        
        dir2 = strcat([dir, '_p_', num2str(1.0-error_rate)]);
        if exist(dir2, 'dir') ~= 7
            mkdir(dir2);
        else
            continue;
        end
        M = 100;
        data = cell(M, 1);
        TL2 = cell(M, 1);
        CD = cell(M, 1);
        for T = 1:(num_experiments/M)
            parfor (TT = 1:M, 10)
                data{TT} = make_data(graph, sample_setting);
                TL2{TT} = Truncated_L2(graph, data{TT}, verbose);
                CD{TT} = CoordinateDescent(graph, data{TT}, verbose);
            end
            for TT = 1:M
                experiment.TL2 = TL2{TT};
                experiment.CD = CD{TT};
                experiment.sample_setting = sample_setting;
                save(strcat([dir2, '/', num2str(TT+(T-1)*M), '.mat']), 'experiment');
            end
        end
        if (verbose)
            fprintf('\n');
        end
    end
end

% 
% result_CG = [];
% result_CD = [];
% observations = [];
% infs_CG = [];
% infs_CD = [];
% times_CG = [];
% times_CD = [];
% for error_rate=0.7
%     fprintf('error_rate=%f\n', error_rate);
%     par = 1.0 - error_rate;
%     t = make_observations(num_edges, sigma, par);
% %     observations = [observations t];
% %     x0 = randn(n, 1)+1.0;
%     x0 = zeros(n, 1);
% %     hold on;
% %     h1 = figure;
% %     [x1, history1] = theirCG(L, B'*t, x0);
%     [x1] = Truncated_L2(B, t, x0, 1, 0.95);
% %     disp(norm(x1, 'inf'));
% %     plot(history1.time, history1.error_inf);
% 
% %     h2 = figure;
%     [x2, history2] = CD(B, L, t, x0);
% %     plot(history2.time, history2.error_inf);
% 
% %     [x3, history3] = IPM(L, B'*t, x0);
% %     plot(history3.time, history3.error_inf);
%     
% %     hold off;
% %     fprintf('experiment %d: %f (CG), %f (CD)\n', i, norm(x1, 'inf'), norm(x2, 'inf'));
% %     save('type1.mat');
% end
% save('type1.mat');
% fprintf('mean_delta=%f', mean(results_us));
% fprintf(', delta=%f (%f)', delta, mean(inf_norms < delta));
% fprintf(', upper_delta=%f (%f)', upper_delta, mean(inf_norms < upper_delta));
% fprintf('\n');
