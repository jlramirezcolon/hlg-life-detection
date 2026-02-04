function f = PlotFeatures2d(F,f1,f2)
    % Get the data for this feature
    class = F.Class;
    abiotic = (class==0);
    environmental = (class==1);
    
    X = F.(f1);
    X_txt = strrep(f1,'_',' ');

    Y = F.(f2);
    Y_txt = strrep(f2,'_',' ');
    
    % Custom colors
    color_abiotic = hex2rgb('#e9963e');  % orange
    color_biotic = hex2rgb('#b2d49a');   % green

    % Plot the two histograms
    f = figure('Color', [1 1 1]);
    hold on;

    % Plot data
    scatter(X(abiotic), Y(abiotic), 80, color_abiotic, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Extraterrestrial');
    scatter(X(environmental), Y(environmental), 80, color_biotic, 's', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Environmental');

    % Add labels, legend, and title
    xlabel(sprintf('%s', X_txt), 'FontSize', 18, 'FontName', 'Helvetica');
    ylabel(sprintf('%s', Y_txt), 'FontSize', 18, 'FontName', 'Helvetica');

    legend({'Abiotic', 'Environmental'}, 'FontSize', 12, 'FontName', 'Helvetica', 'Location', 'northwest');

    % Set axis properties
    set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
    set(gca, 'Box', 'on', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');
    %grid on;

    set(gca,'xscale','log','yscale','log');
    axis square
end
