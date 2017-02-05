function [ B ] = make_graph( i )
    if i == 1
       s = ones(2000, 1);
       B = make_Gnps(2000, 0.2, s); 
    end
    if i == 2
       s = zeros(2000, 1);
       for j = 1:2000
           s(j) = 0.05 + 0.15*(j-1)/1999;
       end
       B = make_Gnps(2000, 0.2, s); 
    end
    if i == 3
       s = ones(20000, 1);
       B = make_Gnps(20000, 0.002, s); 
    end
    if i == 4
       s = zeros(20000, 1);
       for j = 1:20000
           s(j) = 0.005 + 0.015*(j-1)/19999;
       end
       B = make_Gnps(20000, 0.2, s); 
    end
end

