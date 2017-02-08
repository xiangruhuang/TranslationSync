function [experiment] = run(graph, data )
%     data = make_data(graph, sample_setting);
    experiment.TL2 = Truncated_L2(graph, data);
    experiment.CD = CD(graph, data);
%     experiment.sample_setting = sample_setting;
%     save(strcat([dir2, '/', num2str(T), '.mat']), 'experiment');
end

