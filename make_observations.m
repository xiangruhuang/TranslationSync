function [ t ] = make_observations( num_edges, sigma, p )
    t = rand(num_edges, 1);
    choice = rand(num_edges, 1) * 2 - 1.0;
    dist1 = find(choice <= p);
    n1 = size(dist1, 1);
    t(dist1) = randn(n1, 1) * sigma;
end

