function [ params ] = get_algorithm_params( algorithm )
    if strcmp(algorithm, 'Truncated_L2')
        params.eps_max = 1.0;
        params.c = 0.95;
    end
    if strcmp(algorithm, 'CD')
        params.delta = 1e-5;
    end
end

