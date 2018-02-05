function [TL2, CD, plist] = parse_experiments(dir, graph_type, p_range, num_experiments)
    TL2.min = [];
    TL2.max = [];
    TL2.mean = [];
    CD.min = [];
    CD.max = [];
    CD.mean = [];
    plist = [];
    for p = p_range
        subdir = strcat([dir, '/Graph_', num2str(graph_type), '_p_', num2str(p)]);
        if exist(subdir, 'dir') ~= 7
            %fprintf('data does not exist %s\n', subdir);
            continue;
        end
        plist = [plist; p];
        [TL2_p, CD_p] = parse_subdir(subdir, num_experiments);
        TL2.min = [TL2.min; min(TL2_p)];
        TL2.max = [TL2.max; max(TL2_p)];
        TL2.mean = [TL2.mean; mean(TL2_p)];
        CD.min = [CD.min; min(CD_p)];
        CD.max = [CD.max; max(CD_p)];
        CD.mean = [CD.mean; mean(CD_p)];
    end
end

function [ TL2, CD ] = parse_subdir(dirname, num_experiments)
    TL2 = [];
    %min(experiment.TL2.log.error_inf);
    CD = [];
    %min(experiment.CD.log.error_inf);
    for T = 1:num_experiments
        filename = strcat([dirname, '/', num2str(T), '.mat']);
        if exist(filename, 'file') ~= 2
            continue;
        end
        load(filename);
        TL2 = [TL2; min(experiment.TL2.log.error_inf)];
        CD = [CD; min(experiment.CD.log.error_inf)];
    end
end
