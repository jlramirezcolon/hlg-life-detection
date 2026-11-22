% export_distribution_results.m
% Runs molecular_descriptor_distribution.m (from the repository root) and saves its statistics to out/:
%   metrics_distribution.csv   class statistics of every descriptor, and its relative entropy (nats) in the
%                              RE_Biotic_vs_Combined_Abiotic rows: mean = biotic vs combined abiotic,
%                              median = biotic vs abiotic, min = biotic vs abiotic simulations
%   ranksum_p_values.csv       pairwise Wilcoxon rank-sum p values between subcategories (Fig. 1c descriptor)
%   subcategory_ranges_<descriptor>.csv     min, max and range of the Fig. 1c descriptor in each subcategory
%   statistics_table4_<descriptor>.csv      Supplemental Dataset, Statistics Tables, Table 4: subcategory
%                                           statistics, each amino acid counted once per subcategory
run('molecular_descriptor_distribution.m');
rows = {};
for v = 1:numel(variables)
    for c = 1:numel(class)
        s = stats.(variables{v}).(class{c});
        rows(end+1, :) = {variables{v}, class{c}, s.n, s.mean, s.median, s.min, s.max, s.range, s.std}; %#ok<AGROW>
    end
    re = relative_entropy.(variables{v});
    rows(end+1, :) = {variables{v}, 'RE_Biotic_vs_Combined_Abiotic', NaN, re.Biotic_vs_Combined_Abiotic, ...
        re.Biotic_vs_Abiotic_D, re.Biotic_vs_Abiotic_S, NaN, NaN, NaN}; %#ok<AGROW>
end
writetable(cell2table(rows, 'VariableNames', {'variable', 'class', 'n', 'mean', 'median', 'min', 'max', ...
    'range', 'std'}), fullfile('out', 'metrics_distribution.csv'));
if exist('p_values', 'var'), writematrix(p_values, fullfile('out', 'ranksum_p_values.csv')); end

names = cellstr(desired_order);
rng = NaN(numel(names), 3);
for i = 1:numel(names)
    v = y_data_ordered{i};
    if ~isempty(v), rng(i, :) = [min(v), max(v), max(v) - min(v)]; end
end
writetable(table(names(:), rng(:, 1), rng(:, 2), rng(:, 3), 'VariableNames', {'subcategory', 'min', 'max', 'range'}), ...
    fullfile('out', sprintf('subcategory_ranges_%s.csv', property_name)));

order4 = {'Asteroid Parent Body', 'Meteorite', 'Asteroid', 'Hadean Earth', 'Noachian Mars', 'Lunar', 'Titan', ...
          'Europa', 'Hydrothermal', 'Organism', 'Terrestrial biome'};
T4 = cell(numel(order4), 9);
for k = 1:numel(order4)
    i = find(strcmp(names, order4{k}));
    [~, ui] = unique(y_name_ordered{i}, 'stable');    % each amino acid once (by name; D-/L- forms merged)
    u = y_data_ordered{i}(ui);
    nsamp = sum(categories(:) == order4{k});         % all samples in the subcategory
    T4(k, :) = {order4{k}, numel(u), mean(u), std(u), median(u), min(u), max(u), max(u) - min(u), nsamp};
end
writetable(cell2table(T4, 'VariableNames', {'Category', 'AAs_unique', 'Mean', 'Std', 'Median', 'Min', 'Max', ...
    'Range', 'Samples'}), fullfile('out', sprintf('statistics_table4_%s.csv', property_name)));
