function [ B, L ] = make_Gnps( n, p, s )
    ok = false;
    while ~ok
        a = [];
        b = [];
        v = [];
    %     deg = zeros(n, 1);
        num_edges = 0;
        for i = 1:n
            for j = i+1:n
               if rand <= p*s(i)*s(j)
                   num_edges = num_edges + 1;
                   a = [a; num_edges];
                   b = [b; i];
                   v = [v; 1];
                   a = [a; num_edges];
                   b = [b; j];
                   v = [v; -1];
               end
            end
        end

        B = sparse(a, b, v, num_edges, n);
        L = B'*B;
        ok = connected(L);
    end
end

function ok = connected(L)
    n = size(L, 1);
    [row, col, v] = find(L == -1);
    parent = 1:n;
    rank = zeros(n, 1);
    num_edges = size(row);
    count = 0;
    for e = 1:num_edges
        i = row(e);
        j = col(e);
        %find(i)
        pi = i;
        q = [];
        while pi ~= parent(pi)
            q = [q; pi];
            pi = parent(pi); 
        end
        for t = 1:size(q)
            parent(q(t)) = pi;
        end
        %find(j)
        pj = j;
        q = [];
        while pj ~= parent(pj)
            q = [q; pj];
            pj = parent(pj); 
        end
        for t = 1:size(q)
            parent(q(t)) = pj;
        end
        if pi == pj
            continue;
        end
        %merge
        count = count + 1;
        if rank(pi) < rank(pj)
            parent(pi) = pj;
        else
            if rank(pi) > rank(pj)
                parent(pj) = pi;
            else
                parent(pj) = pi;
                rank(pi) = rank(pj) + 1;
            end
        end
    end
    if count == n-1
        ok = true; 
    else
        ok = false;
    end
end
