% simulations.m
%
% This script was used for the following article: Distinguishing life from 
% non-life via molecular frontier orbital energy gaps
%
% Authors: Jose L. Ramirez-Colon (jcol6@gatech.edu), Ziqin Ni (zni47@gatech.edu), 
% and Christopher E. Carr (cecarr@gatech.edu)
% 
% Simulations: This script implements a Bayesian framework to calculate:
%   P(E|A): Probability of observing evidence E given an abiotic sample
%   P(E|B): Probability of observing evidence E given a biotic sample
%   P(B|E): Probability that a sample is biotic given observed evidence (the posterior)
% 
% Summary of sections:
% 1. Imports simulation settings from Excel file, sets a random seed for
%    reproducibility, and loads amino acid (AA) property database.
% 2. Creates output folder, and calculates derived parameters (e.g., number
%    of amino acids to simulate).
% 3. Filters to abiotic AAs, performs P(E|A) Monte Carlo simulations for each 
%    AA count and metric, and builds probability distributions. 
% 4. Filters to biotic AAs, performs P(E|B) Monte Carlo simulations for each 
%    AA count and metric, and builds probability distributions. 
% 5. Calculates marginal probability of evidence for each prior probability
%    P(B) and computes PDFs and CDFs. Applied Bayes' theorem and calculates
%    posterior probability P(B|E).
% 6. Saves results.
%
% Includes 3rd part code for calculating gini coefficient:
% Copyright (c) 2010, Yvan Lengwiler. All rights reserved.

%% Start fresh
clear all; close all; clc;

% Add custom code to the path
addpath('./code');
addpath('./code/3rdparty/gini');

% Run wvar simulation
s = readsettings('sim_settings.xlsx','wvar'); simulate(s);
s = readsettings('sim_settings.xlsx','wmean'); simulate(s);
s = readsettings('sim_settings.xlsx','gini'); simulate(s);

%% Simulate 
function simulate(op)

    % Tell user
    fprintf('Running simulation with the following settings:\n');
    disp(op); fprintf('\n');

    % Time the simulation
    t_start = tic;

    % Initialize random number seed for deterministic performance
    rng(op.seed);

    % Get a properties table from the datafile (fn) with occurrence and
    % frequency information calculated and added.
    P = GetPropertiesTable(op.dataset);

    %% Derived values 

    subfolder = sprintf('%d_%s_%s',op.N_reps,op.aa_sel,op.aa_abm);
    folder = sprintf('./figures-manuscript/%s/%s/%s',op.field,op.metric,subfolder);
    N_N_AA = numel(op.N_AA);       % Number of different AA counts to test
    N_bins = numel(op.metric_edges)-1;    % Number of bins for probability distributions

    % Create output folder if it doesn't exist
    if ~exist(folder), mkdir(folder); end

    %% Compute P(E|A) as function of number of amino acids (nAA)
    
    % Allocate memory for P(E|A)
    sim.P_of_E_given_A.values = NaN([N_N_AA op.N_reps]);    % replicate results
    sim.P_of_E_given_A.count = op.N_AA;                     % Number of AAs (x)
    sim.P_of_E_given_A.edges = op.metric_edges;             % PDF edges (y)
    sim.P_of_E_given_A.hist = NaN([N_N_AA N_bins]);         % Counts
    sim.P_of_E_given_A.pdf = NaN([N_N_AA N_bins]);          % Density
    sim.P_of_E_given_A.cdf = NaN([N_N_AA N_bins]);          % Cumulative distribution
    
    % Limit the amino acids considered to those that are ABIOTIC (A)
    Pf = P(P.Abiotic,:);
        
    % Limit the amino acids to those represented in ABIOTIC category
    Pf = Pf(Pf.Occurrence_Abiotic>0,:);
    
    % Determine number of amino acids from which we can sample
    N_available_AA = height(Pf);
    
    % perform computation for each # of AAs
    for i_AA=1:N_N_AA
        % Get number of amino acids for this computation
        N_AA_i = op.N_AA(i_AA);    
        
        % Initialize temporary storage of feature values in a way compatible
        % with full variable slicing to allow parallel computation
        values = NaN(1,op.N_reps);
    
        % Perform N_reps calculations of P(E|A)
        parfor r=1:op.N_reps
            % Initialize variables
            AA_r = NaN(N_AA_i,1);
            ABD_r = NaN(N_AA_i,1);
            feature_r = NaN;
            % Select which amino acids to include in this computation
            switch op.aa_sel
                case 'uniform'
                    % Randomly choose N_AA_i amino acids from N_available_AA
                    AA_r = randsample(N_available_AA,N_AA_i);
                case 'class_freq'
                    % Randomly choose N_AA_i amino acids from N_available_AA
                    % weighted by occurrence frequency in ABIOTIC category
                    AA_r = randsample(N_available_AA,N_AA_i,true,Pf.Frequency_Abiotic);
                case 'overall_freq'
                    % Randomly choose N_AA_i amino acids from N_available_AA
                    % weighted by overall dataset occurrence frequency
                    AA_r = randsample(N_available_AA,N_AA_i,true,Pf.Frequency);
            end
            
            % Get the values of the target field for the selected amino acids
            val_r = Pf{AA_r,op.field};
            
            % Assign the abundance (ABD_r) for the selected amino acids (AA_r)
            switch op.aa_abm
                case 'uniform'
                    % Abundances should be uniform over the range from the 
                    % minimum to maximum abundances observed in the database
                    ABD_max = Pf.Max_abiotic(AA_r);
                    ABD_min = Pf.Min_abiotic(AA_r);
                    ABD_r = (ABD_max-ABD_min)*rand(1,1)+ABD_min;
                case 'expand'
                    % Abundances should be uniform over the range from scaled 
                    % minimum to maximum abundances observed in the database.
                    % Max is scaled up by a factor of aa_exp
                    ABD_max = Pf.Max_abiotic(AA_r)*aa_exp;
                    % Min is scaled down by a factor of aa_exp
                    ABD_min = Pf.Min_abiotic(AA_r)*1/aa_exp;
                    ABD_r = (ABD_max-ABD_min)*rand(1,1)+ABD_min;
            end
    
            % Compute the feature based on the specified metric
            switch op.metric
                case 'wmean'
                    feature_r = sum(ABD_r .* val_r,'omitnan') / sum(ABD_r,'omitnan');
                case 'wvar'
                    % Calculate weighted mean
                    weighted_mean = sum(ABD_r .* val_r,'omitnan') / sum(ABD_r,'omitnan');
                    % Calculate the weighted variance
                    squared_diff = (val_r - weighted_mean).^2;
                    feature_r = sum(ABD_r .* squared_diff,'omitnan') / sum(ABD_r,'omitnan');
                case 'gini'
                    % Population is equal to amino acid abundance
                    pop = ABD_r;
                    % Calculate gini coefficient with 'Income' equal to 
                    % property value, constrained to be positive
                    feature_r = gini(pop,abs(val_r));
            end
            % Store the feature value
            values(1,r)=feature_r;
        end
        % store the values
        sim.P_of_E_given_A.values(i_AA,:)=values;
        % compute the histogram (counts)
        sim.P_of_E_given_A.hist(i_AA,:)=histcounts(values,op.metric_edges);
        % compute the histogram (frequency or PDF)
        sim.P_of_E_given_A.pdf(i_AA,:)=histcounts(values,op.metric_edges,'Normalization','pdf');
        % compute the cumulative distribution (CDF)
        sim.P_of_E_given_A.cdf(i_AA,:)=histcounts(values,op.metric_edges,'Normalization','cdf');
        % Tell user
        fprintf('Calculated P(E|A) for %d amino acids\n',N_AA_i);
    end
    
    %% Compute P(E|B) as function of number of amino acids (nAA)
    
    % Allocate memory for P(E|B)
    sim.P_of_E_given_B.values = NaN([N_N_AA op.N_reps]);  % replicate results
    sim.P_of_E_given_B.count = op.N_AA;                % Number of AAs (x)
    sim.P_of_E_given_B.edges = op.metric_edges;        % PDF edges (y)
    sim.P_of_E_given_B.hist = NaN([N_N_AA N_bins]);    % Counts
    sim.P_of_E_given_B.pdf = NaN([N_N_AA N_bins]);     % Density
    sim.P_of_E_given_B.cdf = NaN([N_N_AA N_bins]);     % Cumulative distribution
    
    % Limit the amino acids considered to those that are BIOTIC (B)
    Pf = P(P.Biotic,:);
    
    % Limit the amino acids to those represented in BIOTIC category
    Pf = Pf(Pf.Occurrence_Biotic>0,:);
    
    % Determine number of amino acids from which we can sample
    N_available_AA = height(Pf);
    
    % perform computation for each # of AAs
    for i_AA=1:N_N_AA
        % Get number of amino acids for this computation
        N_AA_i = op.N_AA(i_AA);    
        
        % Initialize temporary storage of feature values in a way compatible
        % with full variable slicing to allow parallel computation
        values = NaN(1,op.N_reps);
    
        % Perform N_reps calculations of P(E|B)
        parfor r=1:op.N_reps
            % Initialize variables
            AA_r = NaN(N_AA_i,1);
            ABD_r = NaN(N_AA_i,1);
            feature_r = NaN;
            % Select which amino acids to include in this computation
            switch op.aa_sel
                case 'uniform'
                    % Randomly choose N_AA_i amino acids from N_available_AA
                    AA_r = randsample(N_available_AA,N_AA_i);
                case 'class_freq'
                    % Randomly choose N_AA_i amino acids from N_available_AA
                    % weighted by occurrence frequency in BIOTIC category
                    AA_r = randsample(N_available_AA,N_AA_i,true,Pf.Frequency_Biotic);
                case 'overall_freq'
                    % Randomly choose N_AA_i amino acids from N_available_AA
                    % weighted by overall dataset occurrence frequency
                    AA_r = randsample(N_available_AA,N_AA_i,true,Pf.Frequency);
            end
            
            % Get the values of the target field for the selected amino acids
            val_r = Pf{AA_r,op.field};
            
            % Assign the abundance (ABD_r) for the selected amino acids (AA_r)
            switch op.aa_abm
                case 'uniform'
                    % Abundances should be uniform over the range from the 
                    % minimum to maximum abundances observed in the database
                    ABD_max = Pf.Max_biotic(AA_r);
                    ABD_min = Pf.Min_biotic(AA_r);
                    ABD_r = (ABD_max-ABD_min)*rand(1,1)+ABD_min;
                case 'expand'
                    % Abundances should be uniform over the range from scaled 
                    % minimum to maximum abundances observed in the database.
                    % Max is scaled up by a factor of aa_exp
                    ABD_max = Pf.Max_biotic(AA_r)*aa_exp;
                    % Min is scaled down by a factor of aa_exp
                    ABD_min = Pf.Min_biotic(AA_r)*1/aa_exp;
                    ABD_r = (ABD_max-ABD_min)*rand(1,1)+ABD_min;
            end
    
            % Compute the feature based on the specified metric
            switch op.metric
                case 'wmean'
                    feature_r = sum(ABD_r .* val_r,'omitnan') / sum(ABD_r,'omitnan');
                case 'wvar'
                    % Calculate weighted mean
                    weighted_mean = sum(ABD_r .* val_r,'omitnan') / sum(ABD_r,'omitnan');
                    % Calculate the weighted variance
                    squared_diff = (val_r - weighted_mean).^2;
                    feature_r = sum(ABD_r .* squared_diff,'omitnan') / sum(ABD_r,'omitnan');
                case 'gini'
                    % Population is equal to amino acid abundance
                    pop = ABD_r;
                    % Calculate gini coefficient with 'Income' equal to 
                    % property value, constrained to be positive
                    feature_r = gini(pop,abs(val_r));
            end
            % Store the feature value
            values(1,r)=feature_r;
        end
        % store the values
        sim.P_of_E_given_B.values(i_AA,:)=values;
        % compute the histogram (counts)
        sim.P_of_E_given_B.hist(i_AA,:)=histcounts(values,op.metric_edges);
        % compute the histogram (frequency or PDF)
        sim.P_of_E_given_B.pdf(i_AA,:)=histcounts(values,op.metric_edges,'Normalization','pdf');
        % compute the cumulative distribution (CDF)
        sim.P_of_E_given_B.cdf(i_AA,:)=histcounts(values,op.metric_edges,'Normalization','cdf');
        
        % Tell user
        fprintf('Calculated P(E|B) for %d amino acids\n',N_AA_i);
    end

    % Analyze P(E) and posterior P(B|E) for each prior P(B) desired
    for k=1:numel(op.biotic_priors)
        %% Assume a prior P(B) = 1 - P(A)
        P_of_B = op.biotic_priors(k);

        %% Compute P(E) = P(E|A)*P(A) + P(E|B)*P(B)

        % Allocate memory and/or directly set some results
        sim.P_of_B(k,1) = P_of_B;
        sim.P_of_E(k,1).count = op.N_AA;                % Number of AAs (x)
        sim.P_of_E(k,1).edges = op.metric_edges;        % PDF edges (y)
        sim.P_of_E(k,1).pdf = NaN([N_N_AA N_bins]);     % Density
        sim.P_of_E(k,1).cdf = NaN([N_N_AA N_bins]);     % Cumulative distribution
        % Compute counts, density (PDF), and cumulative (CDF)
        for i_AA=1:N_N_AA
            % compute the frequency or PDF
            sim.P_of_E(k,1).pdf(i_AA,:)= sim.P_of_E_given_A.pdf(i_AA,:)*(1-P_of_B) + ...
                                    sim.P_of_E_given_B.pdf(i_AA,:)*(P_of_B);
            % compute the cumulative distribution (CDF)
            bin_width = op.metric_edges(2)-op.metric_edges(1);
            bin_probs = sim.P_of_E(k,1).pdf(i_AA,:)*bin_width;
            sim.P_of_E(k,1).cdf(i_AA,:)= cumsum(bin_probs);
        end
        
        % Compute P(B|E) = (P(E|B)*P(B))/P(E)
        
        % Allocate memory and/or directly set some results
        sim.P_of_B_given_E(k,1).count = op.N_AA;                % Number of AAs (x)
        sim.P_of_B_given_E(k,1).edges = op.metric_edges;        % PDF edges (y)
        sim.P_of_B_given_E(k,1).pdf = NaN([N_N_AA N_bins]);     % Density
        sim.P_of_B_given_E(k,1).cdf = NaN([N_N_AA N_bins]);     % Cumulative distribution
        % Compute counts, density (PDF), and cumulative (CDF)
        for i_AA=1:N_N_AA
            % compute the PDF using Bayes theorem
            % P(B|E) = (P(E|B)*P(B))/P(E)
            sim.P_of_B_given_E(k,1).pdf(i_AA,:)=sim.P_of_E_given_B.pdf(i_AA,:)*P_of_B./sim.P_of_E(k,1).pdf(i_AA,:);
            % compute the cumulative distribution (CDF)
            bin_width = op.metric_edges(2)-op.metric_edges(1);
            bin_probs = sim.P_of_B_given_E(k,1).pdf(i_AA,:)*bin_width;
            sim.P_of_B_given_E(k,1).cdf(i_AA,:)=cumsum(bin_probs);
        end
    end

    % Save simulation results in 'folder' and include settings 'op' and
    % simulation results 'sim'
    save(fullfile(folder,'sim.mat'),'op','sim');

    %% Calculate simulation time
    t_elapsed = toc(t_start);
    fprintf('Elapsed time %0.2f s\n',t_elapsed);
end

%% Helper functions

% GetPropertiesTable imports sample data and properties from the database
% and identifies which amino acids are biotic and abiotic and determines
% the occurrence and frequency of samples in the database with each 
% amino acid, overall and within the biotic and abiotic classes. A table
% of amino acid properties containing this information is returned so it
% can be used to perform the simulation.
function P = GetPropertiesTable(datafile)
    % Set up import options for our data file; ignore first 3 lines
    opts = detectImportOptions(datafile,'Sheet','MATLAB','NumHeaderLines',2);
    % Force input of column 7 (Pub_chem_ID) using char type (db v1.3+)
    opts.VariableTypes(7) = {'char'};
    % Force input of columns 8 onward using double data type (db v1.3+)
    opts.VariableTypes(8:end) = {'double'};
    % Read in table data, avoid warning about converting column names
    warning off; T = readtable(datafile,opts,'Sheet','MATLAB'); warning on;

    % Get class label from header rows, biotic = 1, abiotic = 0
    [~,~,raw]=xlsread(datafile,'MATLAB');
    bBiotic = logical([raw{2,9:end}]); % Class starts at column 9 (db v1.3+)

    % Get number of samples in the database
    N_Samples = numel(bBiotic);

    % Extract just the sample abundance values
    % Sample abundances start at column 9 (db v1.3+)
    % We only need abundance values, not sample IDs
    Abundance = table2array(T(:, 9:end)); 

    % Calculate frequency of samples with indicated AA measured across db
    T.('Frequency') = T.N/N_Samples; % Frequency of samples with that AA
    
    % Identify min and max abundance measurements for each 
    % amino acid (row) in the table (compute row-wise, ignoring NaNs)
    T.Min_biotic = min(Abundance(:,bBiotic), [], 2, 'omitnan');
    T.Max_biotic = max(Abundance(:,bBiotic), [], 2, 'omitnan');
    T.Min_abiotic = min(Abundance(:,~bBiotic), [], 2, 'omitnan');
    T.Max_abiotic = max(Abundance(:,~bBiotic), [], 2, 'omitnan');

    % Count occurrence separately in Biotic and Abiotic samples
    T.Occurrence_Biotic = sum(~isnan(Abundance(:,bBiotic)),2);
    T.Occurrence_Abiotic = sum(~isnan(Abundance(:,~bBiotic)),2);
    % Compute frequency seprately in Biotic and Abiotic samples
    T.Frequency_Biotic = T.Occurrence_Biotic/sum(bBiotic);
    T.Frequency_Abiotic = T.Occurrence_Abiotic/sum(~bBiotic);

    % Read in properties
    P = readtable(datafile,'Sheet','Properties');

    % Table merging
    % Because different numbers of amino acids were reported in the abundance 
    % table (T) and property table (P), we merged the two tables based on the 
    % symbol of amino acid in each row. 
    newTabCol=zeros(height(P), 1);
    P.("Occurrence") = newTabCol;
    P.("Frequency") = newTabCol;
    P.("Occurrence_Biotic") = newTabCol;
    P.("Frequency_Biotic") = newTabCol;
    P.("Occurrence_Abiotic") = newTabCol;
    P.("Frequency_Abiotic") = newTabCol;
    for i = 1:height(T)
        AA=T.Symbol{i};
        ind = find(strcmp(P.Symbol, AA),1);
        P.Min_abiotic(ind) = T.Min_abiotic(i);
        P.Max_abiotic(ind) = T.Max_abiotic(i);
        P.Min_biotic(ind) = T.Min_biotic(i);
        P.Max_biotic(ind) = T.Max_biotic(i);
        P.Occurrence(ind) = T.N(i);
        P.Frequency(ind) = T.Frequency(i);
        P.Occurrence_Biotic(ind) = T.Occurrence_Biotic(i);
        P.Occurrence_Abiotic(ind) = T.Occurrence_Abiotic(i);
        P.Frequency_Biotic(ind) = T.Frequency_Biotic(i);
        P.Frequency_Abiotic(ind) = T.Frequency_Abiotic(i);
    end

    % Remove entries of amino acids that are not represented in the database. 
    %
    % These are D-cysteine, D-asparagine, D-proline, D-ornithine, L-ornithine, 
    % D-pipecolic acid, L-pipecolic acid, D-alloisoleucine, L-alloisoleucine, 
    % D-5-hydroxylysine, L-5-hydroxylysine, and D-histidine.
    %
    bRemove = (P.Frequency == 0);
    P(bRemove,:)=[];

    % Set the "Biotic" and "Abiotic" column data types to logical.
    P.Biotic = logical(P.Biotic);
    P.Abiotic = logical(P.Abiotic);
end

function settings = readsettings(xlsx,sheet)
    % READSETTINGS Reads settings from a structured Excel spreadsheet
    % assumed to have the values in column 2 and variable names in column 3

    [~,~,raw]=xlsread(xlsx,sheet);
    Value = raw(2:end,2);
    VariableName = cellstr(raw(2:end,3));

    % assign values (we forgo error checking here)
    tmp = [VariableName'; Value']; tmp = reshape(tmp,1,numel(tmp));
    settings = struct(tmp{:});
    % Evaluate fields where this is necessary
    settings.N_AA = eval(settings.N_AA);
    settings.metric_edges = eval(settings.metric_edges);
    settings.metric_range_biotic = eval(settings.metric_range_biotic);
    settings.metric_range_abiotic = eval(settings.metric_range_abiotic);
    settings.metric_range_evidence = eval(settings.metric_range_evidence);
    settings.metric_range_posterior = eval(settings.metric_range_posterior);
    settings.biotic_priors = eval(settings.biotic_priors);
end