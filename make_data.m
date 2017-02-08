function [ data ] = make_data( graph, setting )
    [num_edges, n] = size(graph.B);
    p = 1.0 - setting.error_rate;
    sigma = setting.sigma;
    
    t = rand(num_edges, 1) * 2;
    choice = rand(num_edges, 1);
    dist1 = find(choice <= p);
    n1 = size(dist1, 1);
    t(dist1) = randn(n1, 1) * sigma;
    data.observations = t;
    data.x0 = randn(n, 1);
    data.xstar = ones(n, 1);
end

