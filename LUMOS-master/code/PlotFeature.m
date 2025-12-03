% PlotFeature plots a feature by plotting a histogram of the feature
% and returns the figure handle.
%
% f = PlotFeature(F,feature) plots a feature (column) of F with the
% column label given by the contents of feature.
%
% f = PlotFeature(F,feature,N) uses N bins when plotting a histogram.
%
% f = PlotFeature(F,feature,BINEDGES) uses the bin edges BINEDGES.
%
% f = PlotFeature(F,feature,...,COLORS) uses the colors COLORS for
% meteoritic [COLORS(1,:)] and environmental [COLORS(2,:)] samples.

function f = PlotFeature(F, feature, bins, colors)
    if nargin < 4
        colors(1,:) = hex2rgb('#e9963e'); % orange, for abiotic
        colors(2,:) = hex2rgb('#b2d49a'); % green, for environmental
    end
    
    % Get the data for this feature
    class = F.Class;
    feature_values = F.(feature); 
    feature_txt = strrep(feature, '_', ' ');
    
    % Define colors
    c_abiotic = colors(1,:);
    c_life = colors(2,:);
    
    % Determine if we should use log or linear scaling
    min_val = min(feature_values);
    max_val = max(feature_values);
    
    % Decision criteria for log vs linear scaling
    use_log_scale = false;
    
    % Will use log scale if:
    % 1. All values are positive
    % 2. The range spans at least 2 orders of magnitude
    % 3. The minimum positive value is significantly above zero
    if min_val > 0
        range_orders_of_magnitude = log10(max_val) - log10(min_val);
        if range_orders_of_magnitude >= 2
            use_log_scale = true;
        end
    end
    
    % Filter and prepare data based on scaling type
    if use_log_scale
        % For log scale, use a small positive minimum to avoid log(0)
        min_threshold = max(min_val, eps);
        
        abiotic_values = feature_values(class == 0);
        abiotic_values = abiotic_values(abiotic_values >= min_threshold);
        
        environmental_values = feature_values(class == 1);
        environmental_values = environmental_values(environmental_values >= min_threshold);
        
        plot_min = min_threshold;
        plot_max = max_val;
    else
        % For linear scale, use all values including negatives
        abiotic_values = feature_values(class == 0);
        environmental_values = feature_values(class == 1);
        
        plot_min = min_val;
        plot_max = max_val;
        
        % Add small padding for linear scale
        range_padding = (plot_max - plot_min) * 0.05;
        plot_min = plot_min - range_padding;
        plot_max = plot_max + range_padding;
    end
    
    % Create figure with 2 subplots (histogram on top, boxplot on bottom)
    f = figure('Color', [1 1 1]);
    
    % Create a 2-row grid for the plots
    subplot(2, 1, 1); % Top subplot for histogram
    hold on;
    
    % Create bin edges based on scaling type
    if use_log_scale
        edges = logspace(log10(plot_min), log10(plot_max), bins);
        edges = unique(edges);
    else
        edges = linspace(plot_min, plot_max, bins);
    end
    
    % Plot abiotic histogram
    histogram(abiotic_values, 'BinEdges', edges, 'FaceColor', c_abiotic, ...
        'EdgeColor', c_abiotic*0.5, 'LineWidth', 0.5);
    
    % Plot environmental histogram (bring to front)
    histogram(environmental_values, 'BinEdges', edges, 'FaceColor', c_life, ...
        'EdgeColor', c_life*0.5, 'LineWidth', 0.5);
    
    % Add labels, legend, and title for histogram
    xlabel(sprintf('%s', feature_txt), 'FontSize', 12, 'FontName', 'Helvetica');
    ylabel('Number of Samples', 'FontSize', 12, 'FontName', 'Helvetica');
    legend({'Abiotic', 'Biotic'}, 'FontSize', 12, 'FontName', 'Helvetica', 'Location', 'northwest');
    
    % Set axis properties for histogram
    set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
    set(gca, 'Box', 'on', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');
    
    % Apply appropriate scaling
    if use_log_scale
        set(gca, 'XScale', 'log');
    end
    
    grid on;
    xlim([plot_min, plot_max]);
    
    % Create bottom subplot for horizontal boxplot
    subplot(2, 1, 2);
    
    % For horizontal boxplot, we need to create a combined dataset with a grouping vector
    all_data = [abiotic_values; environmental_values];
    group = [ones(size(abiotic_values)); 2*ones(size(environmental_values))];
    
    % Create horizontal boxplot
    h = boxplot(all_data, group, 'orientation', 'horizontal', ...
        'Labels', {'Abiotic', 'Biotic'}, ...
        'Width', 0.6, 'Symbol', '.');
    
    % Customize boxplot appearance
    set(h, 'LineWidth', 1.5);
    
    % Get all the different parts of the boxplot
    h_boxes = findobj(gca, 'Tag', 'Box');
    h_medians = findobj(gca, 'Tag', 'Median');
    h_outliers = findobj(gca, 'Tag', 'Outliers');
    h_whiskers = findobj(gca, 'Tag', 'Whisker');
    h_caps = findobj(gca, 'Tag', 'Cap');
    
    % Color the boxes correctly - IMPORTANT: The order in the findobj results
    % may not be what we expect, so we identify boxes by their Y position
    for j = 1:length(h_boxes)
        y_position = mean(get(h_boxes(j), 'YData'));
        if abs(y_position - 1) < 0.5  % This is the Abiotic box (position 1)
            patch(get(h_boxes(j), 'XData'), get(h_boxes(j), 'YData'), c_abiotic, 'FaceAlpha', 0.5);
        else  % This is the Environmental box (position 2)
            patch(get(h_boxes(j), 'XData'), get(h_boxes(j), 'YData'), c_life, 'FaceAlpha', 0.5);
        end
    end
    
    % Color the whiskers and caps to match
    for j = 1:length(h_whiskers)
        y_position = mean(get(h_whiskers(j), 'YData'));
        if abs(y_position - 1) < 0.5  % Abiotic whiskers
            set(h_whiskers(j), 'Color', c_abiotic*0.7);
        else  % Environmental whiskers
            set(h_whiskers(j), 'Color', c_life*0.7);
        end
    end
    
    for j = 1:length(h_caps)
        y_position = mean(get(h_caps(j), 'YData'));
        if abs(y_position - 1) < 0.5  % Abiotic caps
            set(h_caps(j), 'Color', c_abiotic*0.7);
        else  % Environmental caps
            set(h_caps(j), 'Color', c_life*0.7);
        end
    end
    
    % Color the medians with a darker version of the fill color
    for j = 1:length(h_medians)
        y_position = mean(get(h_medians(j), 'YData'));
        if abs(y_position - 1) < 0.5  % Abiotic median
            set(h_medians(j), 'Color', c_abiotic*0.4);
        else  % Environmental median
            set(h_medians(j), 'Color', c_life*0.4);
        end
    end
    
    % Color the outliers
    for j = 1:length(h_outliers)
        y_position = get(h_outliers(j), 'YData');
        if abs(y_position - 1) < 0.5  % Abiotic outliers
            set(h_outliers(j), 'MarkerEdgeColor', c_abiotic*0.7);
        else  % Environmental outliers
            set(h_outliers(j), 'MarkerEdgeColor', c_life*0.7);
        end
    end
    
    % Set axis properties for boxplot
    xlabel(sprintf('%s', feature_txt), 'FontSize', 12, 'FontName', 'Helvetica');
    ylabel('Category', 'FontSize', 12, 'FontName', 'Helvetica');
    set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
    set(gca, 'Box', 'on', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');
    
    % Apply appropriate scaling to boxplot too
    if use_log_scale
        set(gca, 'XScale', 'log');
    end
    
    grid on;
    xlim([plot_min, plot_max]); % Set x limits to match histogram
    
    % Adjust spacing between subplots
    set(f, 'Position', get(f, 'Position') .* [1 1 1 1.2]); % Make figure taller
    set(f, 'Units', 'normalized');
    
    % Adjust subplot spacing
    p = get(f, 'Position');
    set(f, 'Position', [p(1), p(2), p(3), p(4)*1.2]);
    
    % Optional: Display scaling information in title or as text
    if use_log_scale
        sgtitle(sprintf('%s (Log Scale)', feature_txt), 'FontSize', 14, 'FontName', 'Helvetica');
    else
        sgtitle(sprintf('%s (Linear Scale)', feature_txt), 'FontSize', 14, 'FontName', 'Helvetica');
    end
end