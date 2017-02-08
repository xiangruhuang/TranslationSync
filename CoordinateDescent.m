function CD = CoordinateDescent( graph, data, verbose )
    log.error_inf = [];
    log.time = [];
    log.params = get_algorithm_params('CD');
    tic;
    x = data.x0;
    xstar = data.xstar;
    t = data.observations;
    L = graph.L;
    [num_edges, n] = size(graph.B);
    deg = full(diag(L));
    maxdeg = max(deg);
    adj = zeros(n, maxdeg);
    dist = zeros(n, maxdeg);
    for i = 1:n
        a = L(i, :);
        adj(i, 1:deg(i)) = find(a < 0);
    end
    deg = zeros(n, 1);
    [row, col, v] = find(L == -1);
    count = 0;
    for e = 1:num_edges*2
        i = row(e);
        j = col(e);
        if i > j
            continue;
        end
        count = count + 1;
        deg(i) = deg(i) + 1;
        deg(j) = deg(j) + 1;
        dist(i, deg(i)) = t(count);
        dist(j, deg(j)) = -t(count);
    end
%     fprintf('CD: iter=     , error_inf=               , time=               s');
    line_len = 70;
    if verbose
        for bbb = 1:line_len
            fprintf(' ');
        end
    end
    opt_mean = mean(xstar);
    for iter = 1:10000
        delta = 0.0;
        for i = randperm(n)       
            x_old = x(i);
            x(i) = median(x(adj(i, 1:deg(i))) + dist(i, 1:deg(i))');
            delta = delta + abs(x(i)-x_old);

        end
        x = x - mean(x) + opt_mean;
        error_inf = norm(x-xstar, 'inf');
        log.error_inf = [log.error_inf; error_inf];
        log.time = [log.time; toc];
        if verbose
            for bbb = 1:line_len
                fprintf('\b');
            end
            fprintf('\tCD\t: iter=%5d, min_error_inf=%15g, time=%15gs', iter, error_inf, toc);
        end
        if (delta < log.params.delta)
            break; 
        end
    end
    if verbose
        fprintf('\n');
    end
    CD.x = x;
    CD.log = log;
end