function draw_figure(graph_types, folder)
    for graph_type = graph_types
        [TL2, CD, plist] = parse_experiments(strcat([folder, '/']), graph_type, 0.0:0.01:1.0, 100);
        TL2.error = 0; %[TL2.mean - TL2.min;TL2.max - TL2.mean];
        CD.error = 0; %[CD.mean - CD.min; CD.max - CD.mean];
        h = figure;
        set(h, 'Visible', 'off');
        plot(plist, TL2.mean, 'b-'); hold on;
        plot(plist, CD.mean, 'k-'); hold off;
        xlabel('p');
        ylabel('$\|x-x^*\|_{\infty}$', 'Interpreter', 'latex');
        title(strcat(['Graph Type ', num2str(graph_type)]));
        legend({'TL2', 'CD'}, 'FontSize', 18);
        saveas(h, strcat([folder, '/', 'Graph_', num2str(graph_type)]), 'epsc');
    end
end
