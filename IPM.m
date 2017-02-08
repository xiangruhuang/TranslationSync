function [x, history] = IPM(B, t, x0)
%     x = x0;
    tic;
    [num_edges, n] = size(B);
    v0 = zeros(num_edges, 1) + max(abs(B*x0 - t)) + 0.01;
    xv0 = [x0; v0];
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
%     options.Diagnostics = 'on';
    options.MaxFunEvals = 10000000;
% %     options.outputFcn = @(x) norm(x(1:n) - sum(x(1:n))/n, 1);
%     [x, fval, exitflag, output] = linprog(f, A, b, [], [], [], [], [x0; v0], options);
    
    history_x = [];
    history_time = [];
    history_error = [];
    fprintf('IPM: error_inf =        ');
    [x, fval, exitflag, output] = fmincon(@(x) f'*x, xv0, A, b, [], [], [], [], [], options);
    x = x(1:n);
    history.x = history_x;
    history.error_inf = history_error;
    history.time = history_time;
    fprintf('\n');
end

function stop = outfun(x,optimValues,state)
stop = false;
   global history_x;
   global history_time;
   global history_error;
   global n;
   switch state
       case 'init'
%            hold on
       case 'iter'
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
           
%            fval_collection = [fval_collection; optimValues.fval];
%            disp(optimValues.fval);
            fprintf('\b\b\b\b\b\b\b%7g', norm(x(1:n)-mean(x(1:n)), 'inf'));
%             disp(size(x));
            history_x = [history_x; x(1:n)'];
            history_error = [history_error; norm(x(1:n)-mean(x(1:n)), 'inf')];
            history_time = [history_time; toc];
           % Concatenate current search direction with 
           % searchdir.
%            searchdir = [searchdir;...
%                         optimValues.searchdirection'];
%            plot(x(1),x(2),'o');
           % Label points with iteration number.
           % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),num2str(optimValues.iteration));
       case 'done'
%            hold off
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

