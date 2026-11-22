% export_weighted_results.m
% Runs abundance_weighted_analysis.m (from the repository root) and saves to out/:
%   metrics_weighted.csv      relative entropy (nats), AUC and direction of every abundance-weighted feature
%   features_per_sample.csv   the abundance-weighted features of every sample
% With include_S1 = true in abundance_weighted_analysis.m, these are the files in out_with_S1/.
run('abundance_weighted_analysis.m');
names = F.Properties.VariableNames(5:width(F))';
writetable(table(names, relative_entropy', auc', relationship', 'VariableNames', ...
    {'feature', 'relent_nats', 'auc', 'direction'}), fullfile('out', 'metrics_weighted.csv'));
writetable(F, fullfile('out', 'features_per_sample.csv'));
