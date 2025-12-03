% molecular_descriptor_distribution.m
%
% This script was used for the following article: Distinguishing life from 
% non-life via molecular frontier orbital energy gaps
%
% Authors: Jose L. Ramirez-Colon (jcol6@gatech.edu), Ziqin Ni (zni47@gatech.edu), 
% and Christopher E. Carr (cecarr@gatech.edu)
%
% Amino Acid Distribution Analysis: This script analyzes the 
% distribution of amino acids across different environments
% (biotic and abiotic) based on 10 molecular descriptors. 
% 
% Inputs:
%  - Amino acid database
%  - Properties data
%
% Summary of sections:
% 1. Imports and organizes amino acid detection data from Excel database
% 2. Processes molecular property data for detected amino acids
% 3. Visualizes the distribution of molecular descriptors across environmental categories
% 4. Calculates statistical summaries for each descriptor by environment class
% 5. Computes and plots relative entropy between biotic and abiotic environments
% 6. Computes and plots p-values pairwise comparisons between categories
% 7. Extended histogram analysis for all molecular descriptors
% 8. Summary table
% 
% Generates three figures published in the manuscript:
%    Figure 1c. Box plots showing amino acid's HOMO-LUMO gap distributions by class
%    Figure 1d. Bar chart of relative entropy values for different molecular descriptors
%    Figure 1e. Statistical significance matrix showing pairwise comparisons between categories
%

%% 1. Organizing Database Data

% Start fresh
clear all; close all; clc;

% Importing data from Excel file
[num_data, txt_data, raw_data] = xlsread('Amino.Acid.Database.v1.3.Release.2025-07-19.xlsx', 'MATLAB_Distributions');
    
% Getting dimensions of the data
[rows, cols] = size(raw_data);
    
% Extracting metadata
id = raw_data(2:end, 1);              % Column 1: Sample ID
class = raw_data(2:end, 2);           % Column 2: Category
subcategories = raw_data(2:end, 3);   % Column 3: Subcategory
conditions = raw_data(2:end, 4);      % Column 4: Experimental conditions  

% Getting amino acid names 
amino_acid_names = raw_data(1, 5:125);
    
% Extracting amino acid detection data 
detection_data = num_data;
    
% Creating categorical arrays 
id_cat = categorical(string(id));
categories_cat = categorical(string(class));
subcategories_cat = categorical(string(subcategories));
conditions_cat = categorical(string(conditions));

% Storing organized data in a structure
organized_data = struct();
organized_data.detection_matrix = detection_data;
organized_data.amino_acids = amino_acid_names;
    
% Store metadata in a separate structure
metadata = struct();
metadata.id = id_cat;
metadata.categories = categories_cat;
metadata.subcategories = subcategories_cat;
metadata.conditions = conditions_cat;
metadata.unique_id = unique(id_cat);
metadata.unique_categories = unique(categories_cat);
metadata.unique_subcategories = unique(subcategories_cat);
metadata.unique_conditions = unique(conditions_cat);
    
% Calculate summary statistics
detected = sum(detection_data(:) == 1);
not_detected = sum(detection_data(:) == 2);
not_reported = sum(detection_data(:) == 3);
low_confidence = sum(detection_data(:) == 4);
    
% Store summary in structure
detection_summary = struct();
detection_summary.detected = detected;
detection_summary.not_detected = not_detected;
detection_summary.not_reported = not_reported;
detection_summary.low_confidence = low_confidence;
    
metadata.detection_summary = detection_summary;

%% 2. Organizing Properties Data

% Import data from Properties sheet
[num_data, txt_data, raw_data] = xlsread('Amino.Acid.Database.v1.3.Release.2025-07-19.xlsx', 'Properties_Distributions');
    
% Extract property names (first column, skipping header)
property_names = raw_data(2:end, 1);
    
% Get the numeric data
n_properties = length(property_names);
n_amino_acids = length(amino_acid_names);
properties_matrix = zeros(n_properties, n_amino_acids);
    
% Fill properties matrix using num_data
if ~isempty(num_data)
    properties_matrix = num_data;
end
    
% Create structure to hold organized data
properties_data = struct();
properties_data.amino_acid_names = amino_acid_names;
properties_data.property_names = property_names;
properties_data.properties_matrix = properties_matrix;
    
% Create a properties table for easier access
properties_table = array2table(properties_matrix, ...
    'RowNames', property_names, ...
    'VariableNames', matlab.lang.makeValidName(amino_acid_names));
    
properties_data.properties_table = properties_table;
    
% Check for any NaN values
nan_count = sum(isnan(properties_matrix(:)));
    if nan_count > 0
        fprintf('\nWarning: Found %d NaN values in the properties data\n', nan_count);
    end

%% 3. Distribution of amino acids by molecular descriptor (Yields Figure 1c.)

% Define the three classes
biotic_categories = {'Terrestrial biome', 'Organism', 'Hydrothermal'};
abiotic_sim_categories = {'Europa', 'Noachian Mars', 'Titan', 'Hadean Earth', 'Interstellar Clouds'};
abiotic_categories = {'Lunar', 'Asteroid', 'Meteorite'};

% Define custom colors for each plot 
% Biotic class plot: Green gradient
color_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
color_map('Terrestrial biome') = [0.0, 0.6, 0.3];   % Forest green
color_map('Organism') = [0.1, 0.7, 0.4];            % Sea green
color_map('Hydrothermal') = [0.2, 0.8, 0.5];        % Light sea green

% Abiotic (Simulations) class plot: Blue gradient
color_map('Europa') = [0.0, 0.3, 0.7];              % Dark blue
color_map('Noachian Mars') = [0.1, 0.4, 0.8];       % Steel blue
color_map('Titan') = [0.2, 0.5, 0.9];               % Medium blue
color_map('Hadean Earth') = [0.3, 0.6, 1.0];        % Sky blue
color_map('Interstellar Clouds') = [0.4, 0.7, 1.0]; % Maya blue

% Abiotic class plot: Orange gradient
color_map('Lunar') = [0.8, 0.4, 0.0];               % Dark orange
color_map('Asteroid') = [0.9, 0.5, 0.1];            % Tiger orange
color_map('Meteorite') = [1.0, 0.6, 0.2];           % Saffron

% Get detection matrix and categorical data
detection_matrix = detection_data;
id = metadata.id;
categories = metadata.categories;

% Select property to plot
property_name = 'HOMO_LUMO_Gap'; % Modify depending on descriptor of interest  
property_idx = find(strcmp(properties_data.property_names, property_name));
if isempty(property_idx)
    error('Property name not found in properties data');
end
property_values = properties_data.properties_matrix(property_idx, :);

% Define desired category order
desired_order = categorical({'Interstellar Clouds', 'Meteorite', ...
    'Asteroid', 'Hadean Earth', 'Noachian Mars', 'Lunar', 'Titan', ...
    'Europa', 'Hydrothermal', 'Terrestrial biome', 'Organism'});

% Initialize data storage for property values by category
y_data_ordered = cell(length(desired_order), 1);
for i = 1:length(desired_order)
    y_data_ordered{i} = [];
end

% Find environments with detected amino acids
has_detected = false(size(id));
for i = 1:size(detection_matrix, 1)
    if any(detection_matrix(i, :) == 1)
        has_detected(i) = true;
    end
end

% Collect data points by category
for i = 1:size(detection_matrix, 1)
    if has_detected(i)
        curr_cat = categories(i);
        cat_idx = find(desired_order == curr_cat);
        if ~isempty(cat_idx)
            detected = detection_matrix(i, :) == 1;
            if any(detected)
                values = property_values(detected);
                y_data_ordered{cat_idx} = [y_data_ordered{cat_idx}, values];
            end
        end
    end
end

% Create figure
figure('Position', [1665, -84, 1627, 929]);

% Calculate heights based on number of categories
total_plot_height = 0.80; 
gap = 0.02;

% Calculate proportional heights based on category counts
biotic_height = total_plot_height * (length(biotic_categories) / (length(biotic_categories) + length(abiotic_sim_categories) + length(abiotic_categories)));
abiotic_simulated_height = total_plot_height * (length(abiotic_sim_categories) / (length(biotic_categories) + length(abiotic_sim_categories) + length(abiotic_categories)));
abiotic_height = total_plot_height * (length(abiotic_categories) / (length(biotic_categories) + length(abiotic_sim_categories) + length(abiotic_categories)));

% Calculate positions of subplots
abiotic_bottom = 0.08;
abiotic_simulated_bottom = abiotic_bottom + abiotic_height + gap;
biotic_bottom = abiotic_simulated_bottom + abiotic_simulated_height + gap;

% Biotic plot: 
subplot('Position', [0.1, biotic_bottom, 0.8, biotic_height]);
hold on;

% Plot each category
for i = 1:length(biotic_categories)
    cat_idx = find(strcmp(cellstr(desired_order), biotic_categories{i}));
    
    if ~isempty(cat_idx) && ~isempty(y_data_ordered{cat_idx})
        current_values = y_data_ordered{cat_idx}(:);
        y_positions = repmat(i, size(current_values));
        
        bx = boxchart(y_positions, current_values, 'Orientation', 'horizontal');
        
        color = color_map(biotic_categories{i});
        bx.BoxFaceColor = color;
        bx.MarkerColor = color;
        bx.WhiskerLineColor = color;
        bx.BoxFaceAlpha = 0.7;
        bx.LineWidth = 1.5;
    end
end

% Set up plot appearance
yticks(1:length(biotic_categories));
yticklabels(biotic_categories);
set(gca, 'XTickLabel', []);
xlim([8.5, 11.5]);
grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
set(gca, 'FontSize', 18);
set(gca, 'Box', 'on');

% Abiotic (Simulated) plot 
subplot('Position', [0.1, abiotic_simulated_bottom, 0.8, abiotic_simulated_height]);
hold on;

% Plot each category
for i = 1:length(abiotic_sim_categories)
    cat_idx = find(strcmp(cellstr(desired_order), abiotic_sim_categories{i}));
    
    if ~isempty(cat_idx) && ~isempty(y_data_ordered{cat_idx})
        current_values = y_data_ordered{cat_idx}(:);
        y_positions = repmat(i, size(current_values));
        
        bx = boxchart(y_positions, current_values, 'Orientation', 'horizontal');
        
        color = color_map(abiotic_sim_categories{i});
        bx.BoxFaceColor = color;
        bx.MarkerColor = color;
        bx.WhiskerLineColor = color;
        bx.BoxFaceAlpha = 0.7;
        bx.LineWidth = 1.5;
    end
end

% Set up plot appearance
yticks(1:length(abiotic_sim_categories));
yticklabels(abiotic_sim_categories);
set(gca, 'XTickLabel', []);
xlim([8.5, 11.5]);
grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
set(gca, 'FontSize', 18);
set(gca, 'Box', 'on');

% Abiotic plot
subplot('Position', [0.1, abiotic_bottom, 0.8, abiotic_height]);
hold on;

% Plot each category
for i = 1:length(abiotic_categories)
    cat_idx = find(strcmp(cellstr(desired_order), abiotic_categories{i}));
    
    if ~isempty(cat_idx) && ~isempty(y_data_ordered{cat_idx})
        current_values = y_data_ordered{cat_idx}(:);
        y_positions = repmat(i, size(current_values));
        
        bx = boxchart(y_positions, current_values, 'Orientation', 'horizontal');
        
        color = color_map(abiotic_categories{i});
        bx.BoxFaceColor = color;
        bx.MarkerColor = color;
        bx.WhiskerLineColor = color;
        bx.BoxFaceAlpha = 0.7;
        bx.LineWidth = 1.5;
    end
end

% Set up plot appearance
yticks(1:length(abiotic_categories));
yticklabels(abiotic_categories);
xlabel(property_name, 'FontSize', 20);
xlim([8.5, 11.5]);
grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
set(gca, 'FontSize', 18);
set(gca, 'Box', 'on');

% Add overall title
sgtitle(['Distribution of ', property_name, ' by Sample Category'], 'FontSize', 24);

%% 4. Calculate statistics for each variable by class

% Reorient properties table
P = rows2vars(properties_table,'VariableNamingRule','preserve');

% Get indices for each class
abiotic_d_idx = find(P.Abiotic_D == 1);                % Abiotic 
abiotic_s_idx = find(P.Abiotic_S == 1);                % Abiotic (Simulated)
biotic_idx = find(P.Biotic == 1);                      % Biotic
combined_abiotic_idx = unique([abiotic_d_idx; abiotic_s_idx]); % Abiotic + Abiotic (Simulated)

% Variables to analyze
variables = {'HOMO_LUMO_Gap', 'Molar_mass', 'MA_index', 'Carbon_number', 'Heat_Formation_kJ_mol', 'Gibbs_Energy_kJ_mol', 'LogP','LogS','TopoPSA','Dipole_moment'};
class = {'Abiotic_D', 'Abiotic_S', 'Biotic', 'Combined_Abiotic'};
indices = {abiotic_d_idx, abiotic_s_idx, biotic_idx, combined_abiotic_idx};

% Initialize results structure
stats = struct();

% Calculate statistics for each variable and category
for v = 1:length(variables)
    var_name = variables{v};
    
    for c = 1:length(class)
        class_name = class{c};
        idx = indices{c};
        
        if ~isempty(idx)
            % Extract values for this variable and category
            values = P.(var_name)(idx);
            
            % Calculate statistics
            stats.(var_name).(class_name).n = length(values);
            stats.(var_name).(class_name).mean = mean(values);
            stats.(var_name).(class_name).median = median(values);
            stats.(var_name).(class_name).min = min(values);
            stats.(var_name).(class_name).max = max(values);
            stats.(var_name).(class_name).range = max(values) - min(values);
            stats.(var_name).(class_name).std = std(values);
            stats.(var_name).(class_name).values = values;
        end
    end
end

%% 5. Calculate Relative Entropy 

relative_entropy = struct();

for v = 1:length(variables)
    var_name = variables{v};
    
    % Get the values for each class
    biotic_values = stats.(var_name).Biotic.values;
    abiotic_d_values = stats.(var_name).Abiotic_D.values;
    abiotic_s_values = stats.(var_name).Abiotic_S.values;
    combined_abiotic_values = stats.(var_name).Combined_Abiotic.values;
    
    % Calculate relative entropy for Biotic vs Abiotic_D
    X_bd = [biotic_values; abiotic_d_values];
    I_bd = [true(size(biotic_values)); false(size(abiotic_d_values))];
    relative_entropy.(var_name).Biotic_vs_Abiotic_D = relativeEntropy(X_bd, I_bd);
    
    % Calculate relative entropy for Biotic vs Abiotic_S
    X_bs = [biotic_values; abiotic_s_values];
    I_bs = [true(size(biotic_values)); false(size(abiotic_s_values))];
    relative_entropy.(var_name).Biotic_vs_Abiotic_S = relativeEntropy(X_bs, I_bs);
    
    % Calculate relative entropy for Biotic vs Combined Abiotic
    X_bc = [biotic_values; combined_abiotic_values];
    I_bc = [true(size(biotic_values)); false(size(combined_abiotic_values))];
    relative_entropy.(var_name).Biotic_vs_Combined_Abiotic = relativeEntropy(X_bc, I_bc);
end

%% 5. Relative entropy plot using Biotic vs. Combined Abiotic values (Yields Figure 1d.)

% Create figure
figure;

% Extract relative entropy values from the structure
rel_entropy_values = zeros(length(variables), 1);
for v = 1:length(variables)
    var_name = variables{v};
    rel_entropy_values(v) = relative_entropy.(var_name).Biotic_vs_Combined_Abiotic;
end

% Sort the values for better visualization
[sorted_values_nats, sort_idx] = sort(rel_entropy_values);
sorted_values = sorted_values_nats / log(2);  % Convert to bits
sorted_variables = variables(sort_idx);

% Create an horizontal bar plot
h = barh(sorted_values);
set(h, 'FaceColor', [244/255, 194/255, 69/255]);
yticks(1:length(variables))
xlim([1e-3 2.3])
yticklabels(sorted_variables)
xlabel("Relative entropy (bits)")
ylabel("Molecular descriptor")
set(gca, 'Xscale','log', 'FontSize', 18)

%% 6. Statistical Significance Matrix (Yields Figure 1e.)

% Create figure with square shape
figure('Position', [100, 100, 450, 450]);

% Define mapping from full names to abbreviations
name_abbr_map = containers.Map({'Hydrothermal', 'Organism', 'Terrestrial biome', ...
                               'Interstellar Clouds', 'Hadean Earth', 'Titan', ...
                               'Noachian Mars', 'Europa', 'Meteorite', 'Asteroid', 'Lunar'}, ...
                               {'Hy', 'O', 'Tb', 'I', 'H', 'T', 'N', 'E', 'M', 'A', 'L'});

% Define custom category order (based on user request)
custom_order_names = {'Hydrothermal', 'Organism', 'Terrestrial biome', ...
                               'Interstellar Clouds', 'Hadean Earth', 'Titan', ...
                               'Noachian Mars', 'Europa', 'Meteorite', 'Asteroid', 'Lunar'};

% Find indices of existing categories that match the custom order
custom_order_indices = [];
for i = 1:length(custom_order_names)
    desired_order_str = cellstr(desired_order);
    idx = find(strcmp(desired_order_str, custom_order_names{i}));
    if ~isempty(idx)
        custom_order_indices = [custom_order_indices, idx];
    end
end

% Find matching indices and get abbreviated labels
[~, ~, idx] = intersect(custom_order_names, desired_order_str, 'stable');
custom_order_indices = idx(idx > 0);
num_valid = length(custom_order_indices);

% Get abbreviated labels
valid_abbr = cell(num_valid, 1);
for i = 1:num_valid
    orig_name = desired_order_str{custom_order_indices(i)};
    valid_abbr{i} = name_abbr_map(orig_name);
end

% Compute p-values matrix
p_values = ones(num_valid);  % Initialize with 1 (not significant)
for i = 1:num_valid
    for j = 1:num_valid
        if i ~= j
            idx_i = custom_order_indices(i);
            idx_j = custom_order_indices(j);
            [p_values(i,j), ~] = ranksum(y_data_ordered{idx_i}, y_data_ordered{idx_j});
        end
    end
end

% Set up axes with more centered position
mainAx = axes('Position', [0.25, 0.25, 0.65, 0.65]);
hold on;
grid_size = num_valid;
xlim([0.5, grid_size+0.5]); ylim([0.5, grid_size+0.5]);
set(gca, 'XTick', 1:grid_size, 'YTick', 1:grid_size, 'XTickLabel', {}, ...
    'YTickLabel', {}, 'FontSize', 18, 'YDir', 'reverse');

% Define colors
colors = struct('p001', [244/255, 194/255, 69/255], ...   % f4c245
                'p05', [158/255, 179/255, 189/255], ...   % 9eb3bd
                'ns', [231/255, 232/255, 232/255]);       % e7e8e8

% Add labels
arrayfun(@(i) text(-0.5, i, valid_abbr{i}, 'HorizontalAlignment', 'right', 'FontSize', 18), 1:num_valid);
% Fix bottom labels alignment
arrayfun(@(i) text(i, grid_size+1.5, valid_abbr{i}, 'HorizontalAlignment', 'center', 'FontSize', 18), 1:num_valid);

% Plot squares
square_size = 0.8;
for i = 1:num_valid
    for j = 1:num_valid
        x_pos = j - square_size/2;
        y_pos = i - square_size/2;
        
        if i ~= j  % Skip diagonal
            % Determine color based on p-value
            if p_values(i, j) < 0.001
                square_color = colors.p001;
            elseif p_values(i, j) < 0.05
                square_color = colors.p05;
            else
                square_color = colors.ns;
            end
        else
            square_color = 'white';  % Diagonal squares are white
        end
        
        % Plot the square
        rectangle('Position', [x_pos, y_pos, square_size, square_size], ...
                 'FaceColor', square_color, 'EdgeColor', 'none');
    end
end

% Finalize plot
grid on;
set(gca, 'GridColor', [0.9, 0.9, 0.9], 'GridAlpha', 0.5);
axis square;
box on;

%% 7. Extended histogram analysis for all molecular descriptors

% Create comprehensive histogram comparison for all descriptors
figure('Position', [100, 100, 1600, 1200]);

% Fixed number of bins for all descriptors
fixed_num_bins = 40; 

for v = 1:length(variables)
    var_name = variables{v};
    
    % Extract values
    biotic_vals = stats.(var_name).Biotic.values;
    abiotic_vals = stats.(var_name).Combined_Abiotic.values;
    all_vals = [biotic_vals; abiotic_vals];
    
    subplot(2, 5, v);
    hold on;
    
    bin_edges = linspace(min(all_vals), max(all_vals), fixed_num_bins + 1);
        
    % Create overlaid histograms with fixed bin edges
    histogram(biotic_vals, 'BinEdges', bin_edges, ...
        'FaceColor', colors(2,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(abiotic_vals, 'BinEdges', bin_edges, ...
       'FaceColor', colors(1,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    
    xlabel(strrep(var_name, '_', ' '), 'FontSize', 10);
    ylabel('Probability', 'FontSize', 10);
    title(sprintf('%s\n(RE: %.3f nats)', strrep(var_name, '_', ' '), ...
        relative_entropy.(var_name).Biotic_vs_Combined_Abiotic), 'FontSize', 11);
    
    if v == 1
        legend({'Biotic', 'Combined Abiotic'}, 'Location', 'best', 'FontSize', 9);
    end
    
    grid on;
    set(gca, 'GridAlpha', 0.3);
    set(gca, 'FontSize', 18);
end

sgtitle('Distribution Comparison: Biotic vs Combined Abiotic', 'FontSize', 16);

%% 8. Generate summary table

% Create summary table
summary_table = table();
summary_table.Descriptor = variables';

% Add statistics for each group
for v = 1:length(variables)
    var_name = variables{v};
    summary_table.Biotic_Mean(v) = stats.(var_name).Biotic.mean;
    summary_table.Biotic_Std(v) = stats.(var_name).Biotic.std;
    summary_table.Abiotic_Mean(v) = stats.(var_name).Combined_Abiotic.mean;
    summary_table.Abiotic_Std(v) = stats.(var_name).Combined_Abiotic.std;
    summary_table.Relative_Entropy(v) = relative_entropy.(var_name).Biotic_vs_Combined_Abiotic;
end

% Display the table
disp(summary_table);