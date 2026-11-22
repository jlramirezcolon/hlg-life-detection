% export_rankings.m
% Chi-square and MRMR scores of every feature (Suppl. Fig. S3c,d) from out/features_per_sample.csv
% (written by export_weighted_results.m), with the same calls as abundance_weighted_analysis.m.
% Saves out/rankings_chi2_mrmr.csv.
F = readtable(fullfile('out', 'features_per_sample.csv'), 'VariableNamingRule', 'preserve');
X = F(:, 5:width(F)); X.Class = F.Class;
[~, sc] = fscchi2(X, 'Class');
[~, sm] = fscmrmr(X, 'Class');
names = X.Properties.VariableNames(1:end-1)';
writetable(table(names, sc(:), sm(:), 'VariableNames', {'feature', 'chi2_score', 'mrmr_score'}), ...
    fullfile('out', 'rankings_chi2_mrmr.csv'));
