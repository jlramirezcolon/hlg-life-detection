% simulation_plots.m
%
% This script was used for the following article: Distinguishing life from 
% non-life via molecular frontier orbital energy gaps
%
% Authors: Jose L. Ramirez-Colon (jcol6@gatech.edu), Ziqin Ni (zni47@gatech.edu), 
% and Christopher E. Carr (cecarr@gatech.edu)
% 
% Simulations Plots: This script generates plots from the Bayesian simulations
% obtained in simulations.m
%
% Generates the following figures published in the manuscript:
%    Figure 4b. Heatmap of P (B|E) as function of the number of amino acids
%               measured.
%    Figure S4. Probability of evidence given abiotic and biotic datasets
%    Figure S5. Confidence of biogenicity at prior of 0.5 across metrics.
%    Figure S6. Effect of prior on evidence distribution and posterior 
%               confidence in biogenicity
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

%% Setup 

clear all; close all; clc;

% Plotting options
fs = 12;                    % Font size
out = fullfile('.','out','figures');      % output folder

letter_width = 8.5;         % Figure export width in inches
letter_height = 11.0;       % Figure export height in inches
contenttype = 'vector';     % Slow but high quality
contenttype = 'image';      % Fast but low quality
yscalelog = false;          % Change y-scale to log

% Sources for plots
base_dir = fullfile('.', 'out', 'simulations', 'HLG_MNDO_pH_7_eV');
suffix_dir = fullfile('1000000_class_freq_uniform_wor', 'sim.mat');

% Construct the cell array
files = {fullfile(base_dir, 'gini', suffix_dir), ...
         fullfile(base_dir, 'wmean', suffix_dir), ...
         fullfile(base_dir, 'wvar', suffix_dir)};

%% Configuration
if ~exist(out,'dir'), mkdir(out); end;

%% Figure S4: Probability of evidence given abiotic and biotic datasets
% Left column: P(E|A), Right Column: P(E|B)
% Rows: Gini, Mean, Weighted Variance

figure('Color',[1 1 1]);
ht = tiledlayout(3,2);

for row=1:3
    fn = files{row};
    load(fn);

    % get label
    label = [op.metric_label ' of ' op.field_label];

    % Plot left column: P(E|A)
    nexttile;
    x = op.N_AA; y = (op.metric_edges(1:end-1)+op.metric_edges(2:end))/2;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,sim.P_of_E_given_A.cdf','FaceAlpha', 0.7, 'EdgeColor','none');
    set(gca,'ylim',op.metric_range_abiotic); set(gca,'FontSize',fs);
    xlabel('Number of Amino Acids');
    ylabel(label);
    colormap jet;
    view(2);
    colorbar;
    %clim([0 1]);
    hold on
    [C1,h1] = contour(X, Y, sim.P_of_E_given_A.cdf', [0.01, 0.10, 0.5, 0.9 0.99],'black','LineWidth',0.5,'LineStyle','-','ShowText','on');
    clabel(C1,h1,'FontSize', fs);
    if row==1, title('P(E|A) (CDF)'); end;
    if yscalelog, set(gca,'yscale','log'); end;

    % Plot right column: P(E|B)
    nexttile;
    x = op.N_AA; y = (op.metric_edges(1:end-1)+op.metric_edges(2:end))/2;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,sim.P_of_E_given_B.cdf','FaceAlpha', 0.7, 'EdgeColor','none');
    set(gca,'xlim',[min(op.N_AA) max(op.N_AA)]); set(gca,'ylim',op.metric_range_biotic); set(gca,'FontSize',fs);
    xlabel('Number of Amino Acids');
    ylabel(label);
    colormap jet;
    view(2);
    colorbar;
    %clim([0 1]);
    hold on
    [C1,h1] = contour(X, Y, sim.P_of_E_given_B.cdf', [0.01, 0.10, 0.5, 0.9 0.99],'black','LineWidth',0.5,'LineStyle','-','ShowText','on');
    clabel(C1,h1,'FontSize', fs);
    if row==1, title('P(E|B) (CDF)'); end;
    if yscalelog, set(gca,'yscale','log'); end;
end

% save
fn = fullfile(out,'FigS4_P_of_E_across_metrics_raster.pdf');
hFig = gcf; set(hFig, 'Units', 'inches'); set(hFig, 'OuterPosition', [0 0 letter_width letter_height]);
exportgraphics(gcf,fn,'contenttype',contenttype);

% save pdf
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS4_P_of_E_across_metrics_vector.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure S5: Confidence of biogenicity at prior of 0.5 across metrics
% Left column: P(E) for P(B)=0.5, Right Column: P(B|E) for P(B)=0.5
% Rows: Gini, Mean, Weighted Variance.

% Specify desired prior P(B)
P_of_B = 0.5;

figure('Color',[1 1 1]);
ht = tiledlayout(3,2);

for row=1:3
    fn = files{row};
    load(fn);

    % get label
    label = [op.metric_label ' of ' op.field_label];

    % Plot left column: P(E) for desired prior
    nexttile;

    % Find the entry for our desired P(B)
    idx_prior = find(sim.P_of_B==P_of_B);

    x = op.N_AA; y = (op.metric_edges(1:end-1)+op.metric_edges(2:end))/2;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,sim.P_of_E(idx_prior).cdf','FaceAlpha', 0.7,'EdgeColor','none');
    set(gca,'xlim',[min(op.N_AA) max(op.N_AA)]); set(gca,'ylim',op.metric_range_evidence); set(gca,'FontSize',fs);
    xlabel('Number of Amino Acids');
    ylabel(label);
    colormap jet;
    view(2);
    colorbar;
    %clim([0 1]);
    hold on
    [C1,h1] = contour(X,Y,sim.P_of_E(idx_prior).cdf',[0.01 0.1 0.5, 0.9, 0.99],'black','LineWidth',1,'LineStyle','-','ShowText','on','ZLocation','zmax');
    clabel(C1,h1,'FontSize', fs);
    if row==1, title(sprintf('P(E) (CDF) for P(B)=%0.3f',P_of_B)); end;
    if yscalelog, set(gca,'yscale','log'); end;

    % Plot right column: P(B|E) for P(B)=0.5
    nexttile;
    
    x = op.N_AA; y = (op.metric_edges(1:end-1)+op.metric_edges(2:end))/2;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,sim.P_of_B_given_E(idx_prior).pdf','FaceAlpha', 0.7, 'EdgeColor','none');
    set(gca,'xlim',[min(op.N_AA) max(op.N_AA)]); set(gca,'ylim',op.metric_range_posterior); set(gca,'FontSize',fs);
    xlabel('Number of Amino Acids');
    ylabel(label);
    colormap jet;
    view(2);
    colorbar;
    %clim([0 1]);
    %hold on
    %[C1,h1] = contour(X,Y,sim.P_of_B_given_E.pdf',[0.001 0.01 0.1 0.5, 0.9, 0.99],'black','LineWidth',2,'LineStyle',':','ShowText','on','ZLocation','zmax');
    %clabel(C1,h1,'FontSize', fs);
    if row==1, title(sprintf('Posterior P(B|E) for P(B)=%0.3f',P_of_B)); end;
    if yscalelog, set(gca,'yscale','log'); end;
end

% save
fn = fullfile(out,'FigS5_Equal_priors_across_metrics_raster.pdf');
hFig = gcf; set(hFig, 'Units', 'inches'); set(hFig, 'OuterPosition', [0 0 letter_width letter_height]);
exportgraphics(gcf,fn,'contenttype',contenttype);

% save pdf
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS5_Equal_priors_across_metrics_vector.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure S6: Effect of prior

% Left column: P(E) for given P(B)
% Right Column: P(B|E) for given P(B)
% Rows: P(B) = 0.1, 0.01, 0.001

% Specify desired prior P(B)
P_of_B = [0.1,0.01,0.001];

figure('Color',[1 1 1]);
ht = tiledlayout(3,2);

% Same file for all 3 rows
fn = files{3};
load(fn);

% get label
label = [op.metric_label ' of ' op.field_label];

% override default plot options
op.metric_range_evidence = [0 0.15];
op.metric_range_posterior = [0 0.15];

for row=1:3
    
    % Plot left column: P(E) for desired prior
    nexttile;

    % Find the entry for our desired P(B)
    idx_prior = find(sim.P_of_B==P_of_B(row));

    x = op.N_AA; y = (op.metric_edges(1:end-1)+op.metric_edges(2:end))/2;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,sim.P_of_E(idx_prior).cdf','FaceAlpha', 0.7,'EdgeColor','none');
    set(gca,'xlim',[min(op.N_AA) max(op.N_AA)]); set(gca,'ylim',op.metric_range_evidence); set(gca,'FontSize',fs);
    xlabel('Number of Amino Acids');
    %ylabel(label);
    ylabel(sprintf('\\textbf{P(B)=%0.3f}\n%s',P_of_B(row),label),'Interpreter','latex')
    colormap jet;
    view(2);
    colorbar;
    %clim([0 1]);
    hold on
    [C1,h1] = contour(X,Y,sim.P_of_E(idx_prior).cdf',[0.01 0.1 0.5, 0.9, 0.99],'black','LineWidth',1,'LineStyle','-','ShowText','on','ZLocation','zmax');
    clabel(C1,h1,'FontSize', fs);
    if row==1, title('P(E) (CDF)'); end;
    %title(sprintf('P(E) for P(B)=%0.3f (CDF)',P_of_B));
    if yscalelog, set(gca,'yscale','log'); end;

    % Plot right column: P(B|E) for P(B)=0.5
    nexttile;
    
    x = op.N_AA; y = (op.metric_edges(1:end-1)+op.metric_edges(2:end))/2;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,sim.P_of_B_given_E(idx_prior).pdf','FaceAlpha', 0.7, 'EdgeColor','none');
    set(gca,'xlim',[min(op.N_AA) max(op.N_AA)]); set(gca,'ylim',op.metric_range_posterior); set(gca,'FontSize',fs);
    xlabel('Number of Amino Acids');
    ylabel(label);
    colormap jet;
    view(2);
    colorbar;
    %clim([0 1]);
    %hold on
    %[C1,h1] = contour(X,Y,sim.P_of_B_given_E.pdf',[0.001 0.01 0.1 0.5, 0.9, 0.99],'black','LineWidth',2,'LineStyle',':','ShowText','on','ZLocation','zmax');
    %clabel(C1,h1,'FontSize', fs);
    if row==1, title('Posterior, P(B|E)'); end;
    % title(sprintf('Posterior P(B|E) for P(B)=%0.3f',P_of_B(row)));
    if yscalelog, set(gca,'yscale','log'); end;
end

% save
fn = fullfile(out,'FigS6_Effect_of_prior_raster.pdf');
hFig = gcf; set(hFig, 'Units', 'inches'); set(hFig, 'OuterPosition', [0 0 letter_width letter_height]);
exportgraphics(gcf,fn,'contenttype',contenttype);

% save pdf
set(gcf, 'Renderer', 'painters');
fn = fullfile(out,'FigS6_Effect_of_prior_vector.pdf');
exportgraphics(gcf, fn, 'ContentType', 'vector');

%% Figure 4b: Posterior P(B|E) individual plots
% P(B|E) for given P(B), with P(B) = 0.1, 0.01, 0.001

% Specify desired prior P(B)
P_of_B = 0.01;

% Same file for all 3 rows
fn = files{3};
load(fn);

% get label
label = [op.metric_label ' of ' op.field_label];

% override default plot options
op.metric_range_evidence = [0 0.15];
op.metric_range_posterior = [0 0.15];

for row=length(P_of_B)
    
    % New figure for each row
    figure('Color',[1 1 1]);

    % Find the entry for our desired P(B)
    idx_prior = find(sim.P_of_B==P_of_B(row));

    % P(B|E) for each P(B)
    x = op.N_AA; y = (op.metric_edges(1:end-1)+op.metric_edges(2:end))/2;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,sim.P_of_B_given_E(idx_prior).pdf','FaceAlpha', 0.7, 'EdgeColor','none');
    set(gca,'xlim',[min(op.N_AA) max(op.N_AA)]); set(gca,'ylim',op.metric_range_posterior); set(gca,'FontSize',fs);
    xlabel('Number of Amino Acids');
    ylabel(label);
    colormap jet;
    view(2);
    colorbar;
    if yscalelog set(gca,'yscale','log');
    end

    % save
    fn = fullfile(out,sprintf('Fig4b_Posterior_for_prior_of_P_of_B_%0.3f.pdf',P_of_B(row)));
    hFig = gcf; set(hFig, 'Units', 'inches'); set(hFig, 'OuterPosition', [0 0 4 4]);
    exportgraphics(gcf,fn,'contenttype',contenttype);
end
