graph_params = get_params();

for graph_type = 1:4
    filename = strcat(['Graph_',num2str(graph_type),'.mat']);
    if (exist(filename, 'file') ~= 2)
        graph = make_Gnps(graph_params{graph_type});
        save(filename, 'graph');
    end
end