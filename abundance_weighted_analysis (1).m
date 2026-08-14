% abundance_weighted_analysis.m
%
% This script was used for the following article: Distinguishing life from 
% non-life via molecular frontier orbital energy gaps
%
% Authors: Jose L. Ramirez-Colon (jcol6@gatech.edu), Ziqin Ni (zni47@gatech.edu), 
% and Christopher E. Carr (cecarr@gatech.edu)
% 
% Abundance-weighted Amino Acid Distribution Analysis: This script 
% amino acid data, computes weighted molecular descriptor statistics, 
% evaluates multiple class separation metrics, and performs machine-learning 
% classification to identify the most discriminative features. 
% 
% Summary of sections:
% 1. Imports amino acid data from an Excel database
% 2. Calculates weighted statistical metrics (mean, variance, Gini coefficient)
%    for molecular descriptors across samples
% 3. Evaluates several class separation methods (relative entropy, AUC, Chi-square, MRMR)
%    to identify the best descriptors for class separation
% 4. Machine learning classification analysis of features
% 
% Generates the following figures published in the manuscript:
%    Figure 2b. Abundance-weighted molecular descriptors improve classification of
%               biotic and abiotic amino acid samples
%    Figure 3.  Performance of abundance-weighted molecular descriptors in
%               biotic–abiotic classification
%    Figure S1. Abundance of Amino Acids Across Biotic and Abiotic Classes
%    Figure S2. Amino Acid Abundance Patterns by Biotic and Abiotic Subcategories
%    Figure S3. Evaluation of Alternative Class Separation Methods Using
%               Abundance-Weighted Amino Acid Descriptors to Differentiate Biotic and
%               Abiotic Samples
%

% Copyright (C) 2026 Planetary eXploration Lab (PXL) - Georgia Institute of Technology
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% COLLABORATION NOTICE: For flight hardware integration, NASA mission 
% proposals, or alternative licensing, please contact PXL at 
% https://www.pxl.earth/.
%
% Includes 3rd part code for calculating gini coefficient:
% Copyright (c) 2010, Yvan Lengwiler. All rights reserved.

%% Subsection 1a. Preliminaries

% Start fresh
clear all; close all; clc;

% Add custom code to the path
addpath('./code');
addpath('./code/3rdparty/gini');

% Define output directory for figures and create it if needed
out=fullfile('.','out','figures');
if ~exist(out,'dir'), mkdir(out); end

%% Subsection 1b. Read in Amino Acid Abundance Data

% Add file path
datafile = fullfile('.','data','Amino.Acid.Database.Release.2026-07-20.xlsx');
% Set up import options for our data file; ignore first 3 lines
opts = detectImportOptions(datafile,'Sheet','MATLAB','NumHeaderLines',2);
% Force input of column 7 (Pub_chem_ID) using char type
opts.VariableTypes(7) = {'char'};
% Force input of columns 8 onward using double data type
opts.VariableTypes(8:end) = {'double'};
% Read in table data, avoid warning about converting column names
warning off; T = readtable(datafile,opts,'Sheet','MATLAB'); warning on;

% Define start column for Class and Class_Label
Start_Col = 9;

% Get header rows data
[~,~,raw]=xlsread(datafile,'MATLAB');
Class = [raw{2,Start_Col:end}]';
Class_Label = raw(1,Start_Col:end)';

% Extract sample identifiers: First five columns are Amino_Acid_Label,
% Amino_Acid, Chirality, Symbol, Molar_Mass, and others are the sample names.
Samples = string(T.Properties.VariableNames(Start_Col:end))';

%% Subsection 1c. Read in properties

% Read in properties
P = readtable(datafile,'Sheet','Properties');
P = removevars(P,{'Biotic','Abiotic'});

%% Section 2. Make Feature Table with weighted statistical metrics

% Make initial table
Class = logical(Class);
F = table(Samples,Class,Class_Label);

% Establish the metric
% Available metrics: 
    % Basic metrics:
        % 'sum': Total abundance of amino acid in each profile.
        % 'wsum'(weighted sum): Total abundance of amino acid in each profile, weighted by the property value.
        % 'wmean'(weighted mean): The average of property values in the profiles, weighted by amino acid abundance.
    % Dispersion metric:
        % 'wvar' (weighted variance): Measures the spread of property values in the profiles, weighted by amino acid abundance.
    % Inequality metric:
        % 'gini': Metric of statistical dispersion measuring the inequality among the property values of the amino acid distribution. 

% Add abundance
F = addFeature(F, T, P, '', 'sum');

% If you want to analyze all metrics
% metrics = {'wsum','wmean','wvar','gini'};

% To analyze non-wsum metrics as discussed in the methods
metrics = {'wmean','wvar','gini'};

% Iterate through each column starting from column 7 in the P table
for k=1:numel(metrics)
    metric_k = metrics{k};
    for i = 6:width(P)
        p = P.Properties.VariableNames{i}; 
        F = addFeature(F, T, P, p, metric_k);
    end
end

%% Subsection 3a. Calculate relative entropy

% Initialize relative entropy values
relative_entropy = zeros(1, width(F) - 4);

% Calculate relative entropy for each column starting from the fifth column
% (ignores first 4 columns which contain metadata)
for i = 5:width(F)
    relative_entropy(i - 4) = relativeEntropy(F{:, i}, F.Class);
end

% Rank features by relative entropy in descending order
[sorted_rel_entropy, sorted_indices_re] = sort(relative_entropy, 'descend');
sorted_values = sorted_rel_entropy / log(2);  % Convert to bits

%% Figure 2b. Generate plots for n top features 

% Histogram displaying distribution of top performing metric/molecular
% descriptor as determined by relative entropy 
% The current setup will yield plots in Figure 2b.
% To visualize results from any other class separation method modify
% "sorted_indices_re" by the respective variable for the method.

varnames = F.Properties.VariableNames(5:width(F));

% Loop for visualization of many features at once
%N_top = 1; % Modify depending on number of top features wanting to be visualized
%top_features = varnames(sorted_indices_re(1:N_top)); 
%for k=1:N_top
    % Plot the feature
    %PlotFeature(F,top_features{k}, 40);
%end

% Specify the base property you want to plot
base_property = 'HLG_wB97XD_TZVP_pH_7_eV'; 
suffixes = {'gini', 'wmean', 'wvar'};  % variants to look for

for k = 1:numel(suffixes)
    feat_name = sprintf('%s_%s', base_property, suffixes{k});

    % Find this feature's position in varnames
    varname_idx = find(strcmp(varnames, feat_name));
    if isempty(varname_idx)
        warning('Feature "%s" not found in varnames.', feat_name);
        continue
    end

    fig = PlotFeature(F, varnames{varname_idx}, 40);
    set(fig, 'Renderer', 'painters');
    fn = fullfile(out, sprintf('Fig2b.%s.pdf', suffixes{k}));
    exportgraphics(gcf, fn, 'ContentType', 'vector');
end

%% Figure 3a. Relative entropy rank of all features

figure('Color', [1 1 1],'Position', [200 100 600 700]); 

% Get variable names for categorization
varnames = F.Properties.VariableNames(5:width(F));
sorted_varnames = varnames(sorted_indices_re);

% Get the number of variables
num_vars = length(sorted_varnames);

% Define colors for each metric type 
wvar_color = [175/255, 48/255, 52/255];    
wmean_color = [236/255, 171/255, 152/255];
gini_color = [218/255, 56/255, 50/255];    
default_color = [0.5, 0.5, 0.5];           

% Process each variable and create cleaned labels
cleaned_varnames = cell(num_vars, 1);
bar_colors = zeros(num_vars, 3);

for i = 1:num_vars
    varname = sorted_varnames{i};
    
    % Determine color based on variable name suffix
    if contains(varname, '_wvar')
        bar_colors(i,:) = wvar_color;
        cleaned_name = strrep(varname, '_wvar', '');
    elseif contains(varname, '_wmean')
        bar_colors(i,:) = wmean_color;
        cleaned_name = strrep(varname, '_wmean', '');
    elseif contains(varname, '_gini')
        bar_colors(i,:) = gini_color;
        cleaned_name = strrep(varname, '_gini', '');
    else
        bar_colors(i,:) = default_color;
        cleaned_name = varname;
    end
    
    % Fix underscores for display
    cleaned_name = strrep(cleaned_name, '_', '\_');
    cleaned_varnames{i} = cleaned_name;
end

% Create reversed list of values and labels for plotting
rev_auc = flip(sorted_values);
rev_colors = flip(bar_colors, 1);
rev_labels = flip(cleaned_varnames);

% Create the horizontal bar chart
b = barh(1:num_vars, rev_auc);

% Apply colors to bars
b.FaceColor = 'flat';
b.CData = rev_colors;
b.FaceAlpha = 0.7;

% Add labels and formatting
ylabel('Relative entropy rank', 'FontSize', 18);
xlabel('Relative entropy (bits)', 'FontSize', 18);
set(gca, 'ytick', 1:num_vars);
set(gca, 'yticklabels', rev_labels);

% Set log scale and font size
set(gca, 'xscale', 'log', 'FontSize', 16);
xlim([1e-1 1e4]);

% Create a custom legend
legend_h = zeros(3,1);
legend_h(1) = patch([0 0 0 0], [0 0 0 0], wvar_color, 'FaceAlpha', 0.7);
legend_h(2) = patch([0 0 0 0], [0 0 0 0], wmean_color, 'FaceAlpha', 0.7);
legend_h(3) = patch([0 0 0 0], [0 0 0 0], gini_color, 'FaceAlpha', 0.7);
set(legend_h, 'Visible', 'off'); 
legend(legend_h, {'Weighted variance', 'Weighted mean', 'Gini coefficient'}, ...
       'Location', 'southeast', 'FontSize', 14);

% Export Figure 3a
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'Fig3a.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure 3b. Two-dimensional scatter plot of weighted variance for HLG (MNDO method) versus MAI

% Yields a 2D scatter plot of any two distributions of results from table F
f1 = 'HLG_MNDO_pH_7_eV_wvar'; % 1st molecular descriptor for a specific metric
f2 = 'MA_index_wvar';        % 2nd molecular descriptor for a specific metric
PlotFeatures2d(F,f1,f2);

% Export Figure 3a
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'Fig3b.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Subsection 3b. Calculate AUC of each feature

% Initialize auc values
auc = zeros(1, width(F) - 4); 

% Calculate ROC for each column starting from the fourth column
for i = 5:width(F)
    % find whether property has better AUC with positive
    [x,y,t,auc_i,optrocpt,suby]=perfcurve(Class,F{:,i},1);
    % or negative relationship
    [xn,yn,tn,auc_in,optrocptn,subyn]=perfcurve(Class,-F{:,i},1);
    if auc_i>auc_in
        auc(i - 4) = auc_i;
        relationship(i-4) = 1;
    else
        auc(i - 4) = auc_in;
        relationship(i-4) = -1;
    end
end

% Rank using relative entropy
[sorted_auc,sorted_indices_auc]=sort(auc,'descend');

%% Figure 3c. Receiver operator characteristic (ROC) curves comparing classification performance of wvar for HLG, MAI, molar mass, and carbon number

% Create figure first
figure('color',[1 1 1]);
hold on  

% Set up custom colors
hlg_color = [78/255, 178/255, 229/255];      % Light blue
ma_color = [60/255, 145/255, 201/255];       % Blue
mass_color = [55/255, 100/255, 170/255];     % Navy
cn_color = [26/255, 37/255, 72/255];         % Deep blue
perfect_color = [205/255, 223/255, 236/255]; % Pastel blue
random_color = [0, 0, 0];                    % Black

% Plot HOMO-LUMO Gap (AUC = 0.968)
field1 = 'HLG_wB97XD_TZVP_pH_7_eV_wvar';
[x1,y1,t1,auc1,optrocpt1,suby1] = perfcurve(Class,F.(field1),1);
plot(x1,y1,'--','LineWidth',2,'Color',hlg_color,'DisplayName', 'HOMO-LUMO Gap')

% Plot MA index (AUC = 0.93)
field2 = 'MA_index_wvar';
[x2,y2,t2,auc2,optrocpt2,suby2] = perfcurve(Class,F.(field2),1);
plot(x2,y2,'--','LineWidth',2,'Color',ma_color,'DisplayName', 'MA index')

% Plot Molar mass (AUC = 0.869)
field3 = 'Molar_mass_g_mol_wvar';
[x3,y3,t3,auc3,optrocpt3,suby3] = perfcurve(Class,F.(field3),1);
plot(x3,y3,'-.','LineWidth',2,'Color',mass_color,'DisplayName','Molar mass')

% Plot Carbon number (AUC = 0.863)
field4 = 'Carbon_number_wvar';
[x4,y4,t4,auc4,optrocpt4,suby4] = perfcurve(Class,F.(field4),1);
plot(x4,y4,'--','LineWidth',2,'Color',cn_color,'DisplayName','Carbon number')

% Add perfect classifier (AUC = 1.0)
x_perfect = [0 0 1];
y_perfect = [0 1 1];
plot(x_perfect,y_perfect,'LineWidth',2,'Color',perfect_color,'DisplayName','Perfect Classifier')

% Add random classifier line (AUC = 0.5)
x_random = [0 1];
y_random = [0 1];
plot(x_random,y_random,':','LineWidth',2,'Color',random_color,'DisplayName','Random Classifier')

% Add labels and title
xlabel('False positive rate', 'FontSize', 18)
ylabel('True positive rate', 'FontSize', 18)
set(gca, 'FontSize', 16);

% Create a nice legend with box
leg = legend('show', 'Location', 'southeast');  
leg.FontSize = 14;
leg.Box = 'on';
leg.EdgeColor = [0.5 0.5 0.5];

% Add grid and customize appearance
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);

% Set axis limits 
xlim([0 1])
ylim([0 1])

hold off

% Export Figure 3c
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'Fig3c.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');


%% Figure S3b. Area Under the Receiver Operating Characteristic Curve rank of all features

figure('Color', [1 1 1],'Position', [200 100 600 700]);

% Get variable names for categorization
varnames = F.Properties.VariableNames(5:width(F));
sorted_varnames = varnames(sorted_indices_auc);

% Get the number of variables
num_vars = length(sorted_varnames);         

% Process each variable and create cleaned labels
cleaned_varnames = cell(num_vars, 1);
bar_colors = zeros(num_vars, 3);

for i = 1:num_vars
    varname = sorted_varnames{i};
    
    % Determine color based on variable name suffix
    if contains(varname, '_wvar')
        bar_colors(i,:) = wvar_color;
        cleaned_name = strrep(varname, '_wvar', '');
    elseif contains(varname, '_wmean')
        bar_colors(i,:) = wmean_color;
        cleaned_name = strrep(varname, '_wmean', '');
    elseif contains(varname, '_gini')
        bar_colors(i,:) = gini_color;
        cleaned_name = strrep(varname, '_gini', '');
    else
        bar_colors(i,:) = default_color;
        cleaned_name = varname;
    end
    
    % Fix underscores for display
    cleaned_name = strrep(cleaned_name, '_', '\_');
    cleaned_varnames{i} = cleaned_name;
end

% Create reversed list of values and labels for plotting
rev_auc = flip(sorted_auc);
rev_colors = flip(bar_colors, 1);
rev_labels = flip(cleaned_varnames);

% Create the horizontal bar chart
b = barh(1:num_vars, rev_auc);

% Apply colors to bars
b.FaceColor = 'flat';
b.CData = rev_colors;
b.FaceAlpha = 0.7;

% Add labels and formatting
ylabel('AUC rank', 'FontSize', 18);
xlabel('Area Under the ROC Curve', 'FontSize', 18);
set(gca, 'ytick', 1:num_vars);
set(gca, 'yticklabels', rev_labels);

% Set scale limits and font size
set(gca, 'xlim',[0 1], 'FontSize', 16);

% Create a custom legend
legend_h = zeros(3,1);
legend_h(1) = patch([0 0 0 0], [0 0 0 0], wvar_color, 'FaceAlpha', 0.7);
legend_h(2) = patch([0 0 0 0], [0 0 0 0], wmean_color, 'FaceAlpha', 0.7);
legend_h(3) = patch([0 0 0 0], [0 0 0 0], gini_color, 'FaceAlpha', 0.7);
set(legend_h, 'Visible', 'off'); 
legend(legend_h, {'Weighted variance', 'Weighted mean', 'Gini coefficient'}, ...
       'Location', 'southeast', 'FontSize', 14);

% Export Figure S3b
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS3b.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Subsection 3c. Rank the features using Chi2
% https://www.mathworks.com/help/stats/fscmrmr.html

% Remove sample description and class lable to avoid leakage
F_filt = removevars(F,["Samples","Class_Label",'_sum']);
% Move Class (response) to final column
F_filt = movevars(F_filt,"Class",'After',width(F_filt));

[indices_chi,scores_chi] = fscchi2(F_filt,'Class');

%% Figure S3c. Chi2 rank of all features

figure('Color', [1 1 1],'Position', [200 100 600 700]);

% Get variable names for categorization
varnames = F_filt.Properties.VariableNames(1:end-1);
sorted_varnames = varnames(indices_chi);

% Get the number of variables
num_vars = length(sorted_varnames);         

% Process each variable and create cleaned labels
cleaned_varnames = cell(num_vars, 1);
bar_colors = zeros(num_vars, 3);

for i = 1:num_vars
    varname = sorted_varnames{i};
    
    % Determine color based on variable name suffix
    if contains(varname, '_wvar')
        bar_colors(i,:) = wvar_color;
        cleaned_name = strrep(varname, '_wvar', '');
    elseif contains(varname, '_wmean')
        bar_colors(i,:) = wmean_color;
        cleaned_name = strrep(varname, '_wmean', '');
    elseif contains(varname, '_gini')
        bar_colors(i,:) = gini_color;
        cleaned_name = strrep(varname, '_gini', '');
    else
        bar_colors(i,:) = default_color;
        cleaned_name = varname;
    end
    
    % Fix underscores for display
    cleaned_name = strrep(cleaned_name, '_', '\_');
    cleaned_varnames{i} = cleaned_name;
end

% Create reversed list of values and labels for plotting
rev_chi = flip(scores_chi(indices_chi));
rev_colors = flip(bar_colors, 1);
rev_labels = flip(cleaned_varnames);

% Create the horizontal bar chart
b = barh(1:num_vars, rev_chi);

% Apply colors to bars
b.FaceColor = 'flat';
b.CData = rev_colors;
b.FaceAlpha = 0.7;

% Add labels and formatting
ylabel('Chi-squared rank', 'FontSize', 18);
xlabel('Predictor importance score', 'FontSize', 18);
set(gca, 'ytick', 1:num_vars);
set(gca, 'yticklabels', rev_labels);

% Set font size
set(gca, 'FontSize', 16);

% Create a custom legend
legend_h = zeros(3,1);
legend_h(1) = patch([0 0 0 0], [0 0 0 0], wvar_color, 'FaceAlpha', 0.7);
legend_h(2) = patch([0 0 0 0], [0 0 0 0], wmean_color, 'FaceAlpha', 0.7);
legend_h(3) = patch([0 0 0 0], [0 0 0 0], gini_color, 'FaceAlpha', 0.7);
set(legend_h, 'Visible', 'off'); 
legend(legend_h, {'Weighted variance', 'Weighted mean', 'Gini coefficient'}, ...
       'Location', 'southeast', 'FontSize', 14);

% Export Figure S3c
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS3c.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Subsection 3d. Rank the features using the minimum redundancy maximum relevance (MRMR) algorithm
% https://www.mathworks.com/help/stats/fscmrmr.html

% Fit using MRMR for classification
[indices_mrmr,scores_mrmr] = fscmrmr(F_filt,'Class','Verbose',1);

%% Figure S3d. MRMR rank of all features

figure('Color', [1 1 1], 'Position', [200 100 600 700]);

% Get variable names for categorization
varnames = F_filt.Properties.VariableNames(1:end-1);
sorted_varnames = varnames(indices_mrmr);

% Get the number of variables
num_vars = length(sorted_varnames);         

% Process each variable and create cleaned labels
cleaned_varnames = cell(num_vars, 1);
bar_colors = zeros(num_vars, 3);

for i = 1:num_vars
    varname = sorted_varnames{i};
    
    % Determine color based on variable name suffix
    if contains(varname, '_wvar')
        bar_colors(i,:) = wvar_color;
        cleaned_name = strrep(varname, '_wvar', '');
    elseif contains(varname, '_wmean')
        bar_colors(i,:) = wmean_color;
        cleaned_name = strrep(varname, '_wmean', '');
    elseif contains(varname, '_gini')
        bar_colors(i,:) = gini_color;
        cleaned_name = strrep(varname, '_gini', '');
    else
        bar_colors(i,:) = default_color;
        cleaned_name = varname;
    end
    
    % Fix underscores for display
    cleaned_name = strrep(cleaned_name, '_', '\_');
    cleaned_varnames{i} = cleaned_name;
end

% Create reversed list of values and labels for plotting
rev_mrmr = flip(scores_mrmr(indices_mrmr));
rev_colors = flip(bar_colors, 1);
rev_labels = flip(cleaned_varnames);

% Create the horizontal bar chart
b = barh(1:num_vars, rev_mrmr);

% Apply colors to bars
b.FaceColor = 'flat';
b.CData = rev_colors;
b.FaceAlpha = 0.7;

% Add labels and formatting
ylabel('MRMR rank', 'FontSize', 18);
xlabel('Predictor importance score', 'FontSize', 18);
set(gca, 'ytick', 1:num_vars);
set(gca, 'yticklabels', rev_labels);

% Set font size
set(gca, 'FontSize', 16);

% Create a custom legend
legend_h = zeros(3,1);
legend_h(1) = patch([0 0 0 0], [0 0 0 0], wvar_color, 'FaceAlpha', 0.7);
legend_h(2) = patch([0 0 0 0], [0 0 0 0], wmean_color, 'FaceAlpha', 0.7);
legend_h(3) = patch([0 0 0 0], [0 0 0 0], gini_color, 'FaceAlpha', 0.7);
set(legend_h, 'Visible', 'off'); 
legend(legend_h, {'Weighted variance', 'Weighted mean', 'Gini coefficient'}, ...
       'Location', 'southeast', 'FontSize', 14);

% Export Figure S3d
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS3d.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure S3a. Class separation methods ranking heatmap 

% Rank relative entropy
rank_re = zeros(size(sorted_rel_entropy));
rank_re(sorted_indices_re) = 1:numel(sorted_rel_entropy);

% Rank using AUC
rank_auc = zeros(size(sorted_auc));
rank_auc(sorted_indices_auc) = 1:numel(sorted_auc);

% Rank using Chi-Square
rank_chi2 = zeros(size(scores_chi));
rank_chi2(indices_chi) = 1:numel(scores_chi);

% Rank using MRMR
rank_mrmr = zeros(size(scores_mrmr));
rank_mrmr(indices_mrmr) = 1:numel(scores_mrmr);

% Combine ranks into a matrix
rank_matrix = [rank_re(:), rank_auc(:), rank_chi2(:), rank_mrmr(:)];
avg_rank = mean(rank_matrix, 2);

% Sort features by average rank (ascending: best-ranked first)
[~, sorted_idx] = sort(avg_rank, 'descend');
rank_matrix = rank_matrix(sorted_idx, :);

% Get sorted feature names
varnames = F.Properties.VariableNames(5:width(F));
sorted_varnames = varnames(sorted_idx);

% Define ranking method names
methodNames = {'Relative Entropy','AUC', 'Chi-Square', 'MRMR'};

% Create heatmap with sorted features
figure('Color', [1 1 1],'Position', [200 100 600 700]); 
h = heatmap(methodNames, strrep(sorted_varnames, '_', '\_'), rank_matrix, ...
            'Colormap', flipud(bone), ...
            'ColorbarVisible', 'on', ...
            'Title', 'Feature Rank Aggregation');

% Flip the y-axis so top-ranked features appear at the top
h.YDisplayData = flipud(h.YDisplayData);

% Customize labels
h.XLabel = 'Ranking Methods';
h.YLabel = 'Features';

% Move the color bar 
h.ColorbarVisible = 'on';

% Set font size
set(gca, 'FontSize', 16);

% Export Figure S3a
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS3a.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure S1a. Total concentration of amino acids (nmol/g) per sample across biotic and abiotic classes

% Apply colors
color_biotic = hex2rgb('#b2d49a');     % Green
color_abiotic = hex2rgb('#e9963e');    % Orange

% Create histogram
conc = F.('_sum');
log10_conc = log10(conc);
BinEdges = linspace(-2,8,51);
f = figure('color',[1 1 1]);
histogram(log10_conc(Class==0),'BinEdges',BinEdges,'FaceColor',color_abiotic,'EdgeColor',color_abiotic*0.5);
hold on;
histogram(log10_conc(Class==1),'BinEdges',BinEdges,'FaceColor',color_biotic,'EdgeColor',color_biotic*0.5);
xlabel('Concentration log_{10} (nmol/g)');
ylabel('# Samples');
legend('Abiotic','Biotic','Location','northwest');
set(gca, 'FontSize', 16);

% Export Figure S1a
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS1a.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure S1b. Scatter plot of AA richness versus total abundances

% First, calculate the amino acid count per sample
num_amino_acids = height(T);
num_samples = width(T) - (Start_Col - 1);
aa_count_per_sample = zeros(num_samples, 1);

% For each sample, count non-NaN amino acid values
for i = 1:num_samples
    sample_col = T{:, i + (Start_Col - 1)};
    aa_count_per_sample(i) = sum(~isnan(sample_col));
end

% Create the scatter plot with your specifications
f = figure('color', [1 1 1]);

% Set font size to 16
set(groot, 'DefaultAxesFontSize', 16, 'DefaultTextFontSize', 16);

% Plot abiotic samples (Class == 0)
scatter(aa_count_per_sample(Class==0), log10_conc(Class==0), 100, ...
    'MarkerFaceColor', color_abiotic, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceAlpha', 0.7, 'Marker', 'o', 'LineWidth', 0.5);
hold on;

% Plot environmental/biotic samples (Class == 1)
scatter(aa_count_per_sample(Class==1), log10_conc(Class==1), 100, ...
    'MarkerFaceColor', color_biotic, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceAlpha', 0.6, 'Marker', 'o', 'LineWidth', 0.5);

% Add labels and legend
xlabel('Number of Amino Acids per Sample', 'FontSize', 16);
ylabel('Total Abundance log_{10} (nmol/g)', 'FontSize', 16);
legend('Abiotic', 'Biotic', 'FontSize', 16);

% Add more transparent grid
grid on;
set(gca, 'GridAlpha', 0.15); % Make grid more transparent

% Export Figure S1b
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS1b.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure S2. Heatmaps of AA abundances across abiotic and biotic subcategories

% Process amino acid data
[~,ia] = unique(T.Amino_Acid);
G = groupsummary(T,'Amino_Acid','sum',Samples);
G = addvars(G,T.Symbol(ia),'After','Amino_Acid','NewVariableNames','Symbol');

% Split and sort by properties
[unique_classes, idx_first] = unique(Class_Label, 'stable');
class_values = Class(idx_first);
life_idx = class_values == 1;
life_data = zeros(size(G,1), sum(life_idx));
nonlife_data = zeros(size(G,1), sum(~life_idx));

% Calculate life/nonlife data
for i = 1:length(unique_classes)
    samples = find(strcmp(Class_Label, unique_classes{i}));
    if class_values(i)
        life_data(:,sum(life_idx(1:i))) = sum(G{:,samples+3}, 2);
    else
        nonlife_data(:,sum(~life_idx(1:i))) = sum(G{:,samples+3}, 2);
    end
end
    
% Sort amino acids by life abundance only
mean_life = mean(life_data, 2);
mean_nonlife = mean(nonlife_data, 2);
[~, sort_idx_life] = sort(mean_life, 'descend');

% Sort classes by abundance (this stays the same)
life_class_sums = sum(life_data, 1);
nonlife_class_sums = sum(nonlife_data, 1);
[~, life_class_sort] = sort(life_class_sums, 'descend');
[~, nonlife_class_sort] = sort(nonlife_class_sums, 'descend');

% Get sorted class labels
life_classes = unique_classes(life_idx);
nonlife_classes = unique_classes(~life_idx);

% Calculate log values
life_log = log10(life_data);
nonlife_log = log10(nonlife_data);

% Calculate global color limits
all_log_data = [life_log(:); nonlife_log(:)];
valid_data = all_log_data(isfinite(all_log_data));
color_min = min(valid_data);
color_max = max(valid_data);

figure('color', [1 1 1], 'Position', [150 150 1200 600]);

% Biotic subplot
subplot(2,1,1);
sorted_life = life_log(sort_idx_life, life_class_sort);
h1 = heatmap(G.Amino_Acid(sort_idx_life), life_classes(life_class_sort), sorted_life', ...
    'ColorLimits', [color_min, color_max]);
xlabel('Amino Acid');
ylabel('Class');
title('Biotic');

% Abiotic subplot 
subplot(2,1,2);
sorted_nonlife = nonlife_log(sort_idx_life, nonlife_class_sort);  % Changed to use sort_idx_life
h2 = heatmap(G.Amino_Acid(sort_idx_life), nonlife_classes(nonlife_class_sort), sorted_nonlife', ...
    'ColorLimits', [color_min, color_max]);
xlabel('Amino Acid');
ylabel('Class');
title('Abiotic');

colormap(turbo);

% Export Figure S2
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS2.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Subsection 4. Machine Learning Classification

% Define predictors and response
labels = F.Class;
%selectedPredictors = 5:width(F); % All predictors
selectedPredictors = {'HLG_MNDO_pH_7_eV_wvar'}; 
features = F{:, selectedPredictors};

% Define models
models = {'Fine Tree', 'Medium Tree', 'Coarse Tree', 'Quadratic Discriminant', 'Binary GLM Logistic Regression'};
numHoldouts = 10; % Number of different holdouts subsets of interest

% Storage for results
results = struct();

for modelIdx = 1:length(models)
    modelType = models{modelIdx};
    savedModels = {}; % Store all the models for each specific type

    % Initialize arrays for all metrics 
    metrics = struct();
    metricNames = {'Accuracy_CrossVal', 'Total_Cost_CrossVal', 'Accuracy_Test', 'Total_Cost_Test', ...
                   'Error_Rate_CrossVal', 'Macro_Precision_CrossVal', 'Micro_Precision_CrossVal', ...
                   'Weighted_Precision_CrossVal', 'Macro_Recall_CrossVal', 'Micro_Recall_CrossVal', ...
                   'Weighted_Recall_CrossVal', 'Macro_F1_Score_CrossVal', 'Micro_F1_Score_CrossVal', ...
                   'Weighted_F1_Score_CrossVal', 'Error_Rate_Test', 'Macro_Precision_Test', ...
                   'Micro_Precision_Test', 'Weighted_Precision_Test', 'Macro_Recall_Test', ...
                   'Micro_Recall_Test', 'Weighted_Recall_Test', 'Macro_F1_Score_Test', ...
                   'Micro_F1_Score_Test', 'Weighted_F1_Score_Test'};
    
    for i = 1:length(metricNames)
        metrics.(metricNames{i}) = zeros(numHoldouts, 1);
    end

    rng('default')
    holdout_partitions = cell(numHoldouts, 1);
    kfold_partitions = cell(numHoldouts, 1);

    for holdoutIdx = 1:numHoldouts
        holdout_partitions{holdoutIdx} = cvpartition(size(features,1), 'HoldOut', 0.2);
        trainValIdx = training(holdout_partitions{holdoutIdx});
        kfold_partitions{holdoutIdx} = cvpartition(sum(trainValIdx), 'KFold', 10);
    end
    
    for holdoutIdx = 1:numHoldouts
        cv_holdout = holdout_partitions{holdoutIdx};
        cv_kfold = kfold_partitions{holdoutIdx};
        trainValIdx = training(cv_holdout);
        testIdx = test(cv_holdout);

        % Cross validations
        cvAccuracies = zeros(10, 1);
        cvPredictions = [];
        cvTrueLabels = [];
        
        % Get indices for the training/validation portion
        trainValIndices = find(trainValIdx);
        
        % Run cross-validation
        for fold = 1:10
            % Get fold indices
            foldTestLogical = test(cv_kfold, fold);
            foldTrainLogical = training(cv_kfold, fold);
            
            foldTestIdx = trainValIndices(foldTestLogical);
            foldTrainIdx = trainValIndices(foldTrainLogical);
            
            % Train model on this fold
            foldModel = trainModel(features(foldTrainIdx,:), labels(foldTrainIdx), modelType);
            
            % Test on fold validation set
            foldPred = predict(foldModel, features(foldTestIdx,:));
            cvAccuracies(fold) = sum(foldPred == labels(foldTestIdx)) / length(foldTestIdx);
            
            % Collect predictions for overall CV metrics
            cvPredictions = [cvPredictions; foldPred];
            cvTrueLabels = [cvTrueLabels; labels(foldTestIdx)];
        end
        
        % Train final model on entire 80% 
        finalModel = trainModel(features(trainValIdx,:), labels(trainValIdx), modelType);
        
        % Save the final model
        savedModels{holdoutIdx} = finalModel;
        
        % Test on the held-out 20%
        testPred = predict(finalModel, features(testIdx,:));
        
        % Calculate all metrics for cross-validation
        [acc_cv, cost_cv, err_cv, prec_cv, rec_cv, f1_cv] = calculateMetrics(cvTrueLabels, cvPredictions);
        metrics.Accuracy_CrossVal(holdoutIdx) = acc_cv;
        metrics.Total_Cost_CrossVal(holdoutIdx) = cost_cv;
        metrics.Error_Rate_CrossVal(holdoutIdx) = err_cv;
        metrics.Macro_Precision_CrossVal(holdoutIdx) = prec_cv.macro;
        metrics.Micro_Precision_CrossVal(holdoutIdx) = prec_cv.micro;
        metrics.Weighted_Precision_CrossVal(holdoutIdx) = prec_cv.weighted;
        metrics.Macro_Recall_CrossVal(holdoutIdx) = rec_cv.macro;
        metrics.Micro_Recall_CrossVal(holdoutIdx) = rec_cv.micro;
        metrics.Weighted_Recall_CrossVal(holdoutIdx) = rec_cv.weighted;
        metrics.Macro_F1_Score_CrossVal(holdoutIdx) = f1_cv.macro;
        metrics.Micro_F1_Score_CrossVal(holdoutIdx) = f1_cv.micro;
        metrics.Weighted_F1_Score_CrossVal(holdoutIdx) = f1_cv.weighted;
        
        % Calculate all metrics for test set
        [acc_test, cost_test, err_test, prec_test, rec_test, f1_test] = calculateMetrics(labels(testIdx), testPred);
        metrics.Accuracy_Test(holdoutIdx) = acc_test;
        metrics.Total_Cost_Test(holdoutIdx) = cost_test;
        metrics.Error_Rate_Test(holdoutIdx) = err_test;
        metrics.Macro_Precision_Test(holdoutIdx) = prec_test.macro;
        metrics.Micro_Precision_Test(holdoutIdx) = prec_test.micro;
        metrics.Weighted_Precision_Test(holdoutIdx) = prec_test.weighted;
        metrics.Macro_Recall_Test(holdoutIdx) = rec_test.macro;
        metrics.Micro_Recall_Test(holdoutIdx) = rec_test.micro;
        metrics.Weighted_Recall_Test(holdoutIdx) = rec_test.weighted;
        metrics.Macro_F1_Score_Test(holdoutIdx) = f1_test.macro;
        metrics.Micro_F1_Score_Test(holdoutIdx) = f1_test.micro;
        metrics.Weighted_F1_Score_Test(holdoutIdx) = f1_test.weighted;
        
    end
    
    % Save all models for this type
    modelName = strrep(modelType, ' ', '_');
    out_path = fullfile('.','out','ml');
    if ~exist(out_path,'dir'), mkdir(out_path); end
    out_fn = fullfile(out_path,[modelName '_models.mat']);
    save(out_fn, 'savedModels');
    
    % Store results
    results.(modelName) = metrics;
    
end

% Save the overall results
overall_fn = fullfile('.','out','ml','holdout_results.mat');
save(overall_fn, 'results', 'F', 'features', 'labels');

% Create Excel file with results

% Get all model names
modelNames = fieldnames(results);

% Initialize master table
allData = table();

for i = 1:length(modelNames)
    modelName = modelNames{i};
    modelData = results.(modelName);
    
    % Get all metric names for this model
    metricNames = fieldnames(modelData);
    
    % Create a table for this model
    for holdoutIdx = 1:numHoldouts
        % Create a row for each holdout iteration
        row = table();
        row.Model = {strrep(modelName, '_', ' ')}; 
        row.Holdout_Iteration = holdoutIdx;
        
        % Add all metrics
        for j = 1:length(metricNames)
            metricName = metricNames{j};
            row.(metricName) = modelData.(metricName)(holdoutIdx);
        end
        
        % Append to master table
        allData = [allData; row];
    end
end

% Export detailed results
MLM_All_fn = fullfile('.','out','ml','MLM_All_results.xlsx');
writetable(allData, MLM_All_fn, 'Sheet', 'Results');

% Create summary table
summaryData = table();
for i = 1:length(modelNames)
    modelName = modelNames{i};
    modelData = results.(modelName);
    
    summaryRow = table();
    summaryRow.Model = {strrep(modelName, '_', ' ')};
    summaryRow.CrossVal_Accuracy_Mean = mean(modelData.Accuracy_CrossVal);
    summaryRow.CrossVal_Accuracy_Std = std(modelData.Accuracy_CrossVal);
    summaryRow.Test_Accuracy_Mean = mean(modelData.Accuracy_Test);
    summaryRow.Test_Accuracy_Std = std(modelData.Accuracy_Test);
    summaryRow.CrossVal_F1_Macro_Mean = mean(modelData.Macro_F1_Score_CrossVal);
    summaryRow.CrossVal_F1_Macro_Std = std(modelData.Macro_F1_Score_CrossVal);
    summaryRow.Test_F1_Macro_Mean = mean(modelData.Macro_F1_Score_Test);
    summaryRow.Test_F1_Macro_Std = std(modelData.Macro_F1_Score_Test);
    
    summaryData = [summaryData; summaryRow];
end

% Export summary
writetable(summaryData, MLM_All_fn, 'Sheet', 'Summary');

%% Helper Function

function z = relativeEntropy(x, labels)
    % Replicates Predictive Maintenance Toolbox's relativeEntropy
    % Assumes Gaussian-distributed data; labels is logical
    x1 = x(labels);
    x2 = x(~labels);
    m1 = mean(x1); v1 = var(x1);
    m2 = mean(x2); v2 = var(x2);
    z = 0.5*((v2/v1)+(v1/v2)-2) + 0.5*(v1+v2)*((m1-m2)^2)/(v1*v2);
end

% Helper functions
function model = trainModel(features, labels, modelType)
  
    switch modelType
        case 'Fine Tree'
            model = fitctree(features, labels, 'MinLeafSize', 4);
        case 'Medium Tree'
            model = fitctree(features, labels, 'MinLeafSize', 12);
        case 'Coarse Tree'
            model = fitctree(features, labels, 'MinLeafSize', 36);
        case 'Quadratic Discriminant'
            model = fitcdiscr(features, labels, 'DiscrimType', 'quadratic');
        case 'Binary GLM Logistic Regression'
            model = fitclinear(features, labels, 'Learner', 'logistic');
    end
end

function [accuracy, totalCost, errorRate, precision, recall, f1] = calculateMetrics(trueLabels, predictions)
    
    % Handle case where predictions might be empty
    if isempty(predictions)
        accuracy = 0; totalCost = length(trueLabels); errorRate = 1;
        precision = struct('macro', 0, 'micro', 0, 'weighted', 0);
        recall = struct('macro', 0, 'micro', 0, 'weighted', 0);
        f1 = struct('macro', 0, 'micro', 0, 'weighted', 0);
        return;
    end
    
    % Confusion matrix
    try
        [C, order] = confusionmat(trueLabels, predictions);
    catch
        accuracy = sum(predictions == trueLabels) / length(trueLabels);
        errorRate = 1 - accuracy;
        totalCost = sum(predictions ~= trueLabels);
        precision = struct('macro', accuracy, 'micro', accuracy, 'weighted', accuracy);
        recall = struct('macro', accuracy, 'micro', accuracy, 'weighted', accuracy);
        f1 = struct('macro', accuracy, 'micro', accuracy, 'weighted', accuracy);
        return;
    end
    
    % Metrics
    accuracy = sum(predictions == trueLabels) / length(trueLabels);
    errorRate = 1 - accuracy;
    totalCost = sum(predictions ~= trueLabels);
    
    % For binary classification
    if size(C, 1) == 2
        TP = C(2,2); TN = C(1,1); FP = C(1,2); FN = C(2,1);
        
        % Handle division by zero
        if (TP + FP) == 0
            prec_pos = 0;
        else
            prec_pos = TP / (TP + FP);
        end
        
        if (TN + FN) == 0
            prec_neg = 0;
        else
            prec_neg = TN / (TN + FN);
        end
        
        if (TP + FN) == 0
            rec_pos = 0;
        else
            rec_pos = TP / (TP + FN);
        end
        
        if (TN + FP) == 0
            rec_neg = 0;
        else
            rec_neg = TN / (TN + FP);
        end
        
        % Macro averages
        precision.macro = (prec_pos + prec_neg) / 2;
        recall.macro = (rec_pos + rec_neg) / 2;
        
        if (precision.macro + recall.macro) == 0
            f1.macro = 0;
        else
            f1.macro = 2 * (precision.macro * recall.macro) / (precision.macro + recall.macro);
        end
        
        % Micro averages 
        precision.micro = accuracy;
        recall.micro = accuracy;
        f1.micro = accuracy;
        
        % Weighted averages
        total = sum(C(:));
        if total > 0
            weight_neg = sum(C(1,:)) / total;
            weight_pos = sum(C(2,:)) / total;
            precision.weighted = prec_neg * weight_neg + prec_pos * weight_pos;
            recall.weighted = rec_neg * weight_neg + rec_pos * weight_pos;
            
            if (precision.weighted + recall.weighted) == 0
                f1.weighted = 0;
            else
                f1.weighted = 2 * (precision.weighted * recall.weighted) / (precision.weighted + recall.weighted);
            end
        else
            precision.weighted = 0;
            recall.weighted = 0;
            f1.weighted = 0;
        end
    else
        % Multiclass case
        n_classes = size(C, 1);
        class_precision = zeros(n_classes, 1);
        class_recall = zeros(n_classes, 1);
        
        for i = 1:n_classes
            if sum(C(:, i)) > 0
                class_precision(i) = C(i, i) / sum(C(:, i));
            end
            if sum(C(i, :)) > 0
                class_recall(i) = C(i, i) / sum(C(i, :));
            end
        end
        
        precision.macro = mean(class_precision);
        recall.macro = mean(class_recall);
        precision.micro = accuracy;
        recall.micro = accuracy;
        
        if (precision.macro + recall.macro) == 0
            f1.macro = 0;
        else
            f1.macro = 2 * (precision.macro * recall.macro) / (precision.macro + recall.macro);
        end
        f1.micro = accuracy;
        
        % Weighted averages
        total = sum(C(:));
        if total > 0
            class_weights = sum(C, 2) / total;
            precision.weighted = sum(class_precision .* class_weights);
            recall.weighted = sum(class_recall .* class_weights);
            
            if (precision.weighted + recall.weighted) == 0
                f1.weighted = 0;
            else
                f1.weighted = 2 * (precision.weighted * recall.weighted) / (precision.weighted + recall.weighted);
            end
        else
            precision.weighted = 0;
            recall.weighted = 0;
            f1.weighted = 0;
        end
    end
end