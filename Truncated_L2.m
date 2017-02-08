function TL2 = Truncated_L2( graph, data, verbose )
%     B0=graph.B;
%     t0=data.observations;
%     x0=data.x0;

    tic;
    log.params = get_algorithm_params('Truncated_L2');
    eps_max = log.params.eps_max;
    c = log.params.c;
    log.error_inf = [];
    log.time = [];
    
    x = data.x0;
    xstar = data.xstar;
    B = graph.B;
    t = data.observations;
    eps = eps_max;
    opt_mean = mean(data.xstar);
    line_len = 92;
    if verbose
        for bbb = 1:line_len
            fprintf(' ');
        end
    end
    iter = 0; min_error_inf = 100000;
    while eps > 1e-5
        iter = iter + 1;
        x = CG(B'*B, B'*t, x);
        x = x - mean(x) + opt_mean;
        nnz_indices = find(abs(B*x-t) <= eps);
        log.time = [log.time; toc];
        error_inf = norm(x-xstar, 'inf');
        log.error_inf = [log.error_inf; error_inf];
        if min_error_inf > error_inf
            min_error_inf = error_inf;
        end
        if verbose
            for bbb = 1:line_len
                fprintf('\b');
            end
            fprintf('Truncated_L2\t: iter=%5d, min_error_inf=%15g, time=%15gs, nnz=%7d', iter, min_error_inf, toc, size(nnz_indices, 1));
        end
        B = B(nnz_indices, :);
        t = t(nnz_indices, :);
        eps = eps*c;
    end
    if verbose
        fprintf('\n');
    end
    TL2.x = x;
    TL2.log = log;
end

function [x] = CG(A, b, x0)
%     history.error_inf = [];
%     history.time = [];
%     history.x = [];
%     global log;
    x = x0;
    r = b - A * x0;
    p = r;
    rsold = r' * r;
    n = size(A, 1);
    
    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        x = x - mean(x) + 1.0;
%         infs = [infs; norm(zz, 'inf')];
%         times = [times; toc];
%         log.error_inf = [log.error_inf; norm(x, 'inf')];
%         log.time = [log.time; toc];
%         history.x = [history.x; x'];
        
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
%     fprintf('\n');
end