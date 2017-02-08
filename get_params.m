function [ params ] = get_params( )
    params = cell(4, 1);

    for i = 1:4
        params{i}.type = i;
    end
    
    params{1}.n = 2000;
    params{1}.p = 0.2;
    params{1}.s = ones(params{1}.n, 1);
    
    
    params{2}.n = 2000;
    params{2}.p = 0.2;
    params{2}.s = zeros(params{2}.n, 1);
    for j = 1:params{2}.n
        params{2}.s(j) = 0.1 + 0.15*(j-1)/(params{2}.n-1);
    end
    
    params{3}.n = 20000;
    params{3}.p = 0.002;
    params{3}.s = ones(params{3}.n, 1);
    
    params{4}.n = 20000;
    params{4}.p = 0.2;
    params{4}.s = zeros(params{4}.n, 1);
    for j = 1:params{4}.n
        params{4}.s(j) = 0.05 + 0.075*(j-1)/(params{4}.n-1);
    end
    
end

