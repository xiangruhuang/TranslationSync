function [ x, infs, times ] = CD( B, L, t, x0 )
    %\|Bx - t\|^2_2 - mu log(1^T x)
    infs = [];
    times = [];
    tic;
    x = x0;
    [num_edges, n] = size(B);
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
    
    dists = zeros(maxdeg, 1);
    for time = 1:1000
        for i = randperm(n)       
            x(i) = median(x(adj(i, 1:deg(i))) + dist(i, 1:deg(i))');
%             if isnan(x(i))
%                disp('==============================='); 
%             end
%             x(adj(i, 1:deg(i))) + dist(i, 1:deg(i));
%             for s = 1:deg(i)
%                 j = adj;
%                 tij = dist{i}(s);
%                 dists = [dists; x(j) + tij];
%             end
%             x(i) = median(dists(1:deg(i)));
        end
        z = mean(x)/n;
        x = x - z;
%         fprintf('iter=%d, inf_norm=%f\n', time, norm(x, 'inf'));
        infs = [infs; norm(x, 'inf')];
        times = [times; toc];
%         disp('==================================================================\n');
    end
end

