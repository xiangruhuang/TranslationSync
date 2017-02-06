function [x, output] = IPM(B, t, x0)
%     x = x0;
    [num_edges, n] = size(B);
    v0 = zeros(num_edges, 1) + max(abs(B*x0 - t)) + 0.01;
    f = ones(n+num_edges, 1);
    f(1:n) = zeros(n, 1);
    b = [t; -t];
    A = [B, -speye(num_edges); -B, -speye(num_edges)];
%     options = optimoptions('linprog','Algorithm','interior-point');
%     func = @(x) norm(x, 1);
    options = optimoptions(@fmincon);
%     options.OutputFcn = @(x) func(x, n);
    options.OutputFcn = @outfun;
    options.Display = 'iter';
    options.TolCon = 1e-5;
    options.TolFun = 1e-5;
    options.MaxIter = 10000;
    options.Diagnostics = 'on';
    options.MaxFunEvals = 10000000;
% %     options.outputFcn = @(x) norm(x(1:n) - sum(x(1:n))/n, 1);
%     [x, fval, exitflag, output] = linprog(f, A, b, [], [], [], [], [x0; v0], options);
    xv0 = [x0; v0];
%     history.x = [];
%     history.fval = [];
    [x, fval, exitflag, output] = fmincon(@(x) f'*x, xv0, A, b, [], [], [], [], [], options);
    x = x(1:n);
end

function stop = outfun(x,optimValues,state)
stop = false;
   
   switch state
       case 'init'
           hold on
       case 'iter'
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
%            history.fval = [history.fval; optimValues.fval];
            disp(optimValues.fval);
            disp(norm(x(1:2000)-mean(x(1:2000)), 'inf'));
%            history.x = [history.x; x];
           % Concatenate current search direction with 
           % searchdir.
%            searchdir = [searchdir;...
%                         optimValues.searchdirection'];
%            plot(x(1),x(2),'o');
           % Label points with iteration number.
           % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),num2str(optimValues.iteration));
       case 'done'
           hold off
       otherwise
   end
end
% function [ x ] = IPM( B, t, x0 )
%     x = x0;
%     [num_edges, n] = size(B);
%     v = zeros(num_edges, 1) + max(abs(B*x - t)) + 0.01;
%     eta = 0.000001;
%     mu = 1.0;
% %     x0
% %     v
% %     max(abs(B*x - t))
%     for time = 1:10000
%         gradv = ones(num_edges, 1);
%         gradx = zeros(n, 1);
%         bottom1 = B*x - t + v;
%         bottom2 = v - B*x + t;
%         for e = 1:num_edges
%             gradv(e) = gradv(e) - mu*(1.0/bottom1(e) - 1.0/bottom2(e));
%             gradx = gradx - B(e,:)'*mu*(1.0/bottom1(e) - 1.0/bottom2(e));
%         end
% %         gradx
% %         gradv
%         x = x + eta*gradx;
%         v = v + eta*gradv;
%         [norm(gradx, 2), norm(gradv, 2)]
%         if (norm(gradx, 2) + norm(gradv, 2) < 0.1)
%             mu = mu / 1.5;
%         end
%         disp(sum(v));
%     end
% end

