function [x, infs, times] = theirCG(A, b, x)
    infs = [];
    times = [];
    tic;
    r = b - A * x;
    p = r;
    rsold = r' * r;
    n = size(A, 1);
    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        zz = x - sum(x)/n;
        infs = [infs; norm(zz, 'inf')];
        times = [times; toc];
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end