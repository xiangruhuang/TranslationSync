function [ B ] = make_Gnp( n, p )
    B = [];
    for i = 1:n 
        for j = i+1:n
           if rand <= p
               row = zeros(1, n);
               row(i) = 1;
               row(j) = -1;
               B = [B; row];
           end
        end
    end
end

