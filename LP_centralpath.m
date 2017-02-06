function [ xs ] = LP_centralpath( c, A, b, al, x0, a, be )
    if ~exist('a', 'var')
        a = 0.05;
    end
    if ~exist('be', 'var')
        be = 0.9;
    end
    m = size(A, 1);
    n = size(A, 2);
    for i = 1:m
        if (b(i) - A(i,:)*x0 <= 0)
            disp('bi - Aix <= 0');
            return;
        end
    end
%     for i = 1:n
%         if (x0(i) <= 0)
%             disp('xi <= 0');
%             return;
%         end
%     end
    x = x0;
    t = 0;
    xs = x;
    while t < 1e10
        g = t*c;
        for i = 1:m
            g = g + A(i,:)'/(b(i) - A(i,:)*x);
        end
%         for i = 1:n
%             g(i) = g(i) - 1/x(i);
%         end
        h = zeros(n);
        for i = 1:m
            h = h + A(i,:)'*A(i,:)/(b(i) - A(i,:)*x)^2;
        end
%         for i = 1:n
%             h(i,i) = h(i,i) + 1/x(i)^2; 
%         end
        v = -inv(h)*g;
        eta = 1;
%         [f(x + eta*v) , f(x)]
        while true
            xplus = x+eta*v;
            fplus = t*c'*xplus; f = t*c'*x;
            for i = 1:m
                if (b(i) - A(i,:)*xplus <= 0)
                    fplus = 1000000000;
                    f     =-1000000000;
                    break;
                end
                fplus = fplus - log(b(i) - A(i,:)*xplus);
                f     = f     - log(b(i) - A(i,:)*x    );
            end
%             for i = 1:n
%                 if (xplus(i) <= 0)
%                     fplus = 1000000000;
%                     f     =-1000000000;
%                     break;
%                 end
%                 fplus = fplus - log(xplus(i));
%                 f     = f     - log(x    (i));
%             end
            if fplus < f + a*eta*dot(g,v)
                flag = false;
                if (~flag)
                    break;
                end
            end
            eta = eta * be
        end
        x = x + eta*v;
        norm(x-mean(x), 1)
        if abs(g'*v) < 1e-6
            xs = [xs x];        

            if (t == 0)
                t = 1;
            else
                t = t*(1+al);
            end
        end
    end
end

