% addFeature.m: Computes abundance-weighted statistical features from molecular properties
%
% Calculates various metrics (sum, weighted sum, weighted mean, weighted variance, 
% or Gini coefficient) for a specified property across sample datasets, using 
% molecular abundance as weights. Adds the computed feature as a new column to 
% the feature table.
%
% Syntax: FT = addFeature(FT,DT,PT,property,metric)
%
% Inputs:
%   FT       Feature table containing sample names
%   DT       Data table with abundance values for each sample
%   PT       Property table containing molecular property values
%   property Name of the property to analyze (string)
%   metric   Statistical metric to compute: 'sum', 'wsum', 'wmean', 'wvar', or 'gini'
%
% Output:
%   FT: Updated feature table with new column named 'property_metric'

function FT = addFeature(FT,DT,PT,property,metric)
    N_Samples = size(FT.Samples,1);
    feature = zeros(N_Samples,1);
    for k=1:N_Samples
        % Get abundance
        abundance = DT.(FT.Samples{k});
        % Calculate appropriate metric
        switch metric
            case 'sum'
                feature(k) = sum(abundance,'omitnan');
            case 'wsum'
                % Get property values for each data table entry
                val = LookupProperties(DT,PT,property);
                val = cell2mat(val);
                % Compute abundance weighted sum
                feature(k) = sum(abundance.*val,'omitnan');
            case 'wmean'
                % Retrieve property values
                val = LookupProperties(DT,PT,property);
                val = cell2mat(val);
                % Calculate weighted mean
                feature(k) = sum(abundance .* val,'omitnan') / sum(abundance,'omitnan');
            case 'wvar'
                % Retrieve property values
                val = LookupProperties(DT,PT,property);
                val = cell2mat(val);
                % Calculate weighted mean
                weighted_mean = sum(abundance .* val,'omitnan') / sum(abundance,'omitnan');
                % Calculate the weighted variance
                squared_diff = (val - weighted_mean).^2;
                feature(k) = sum(abundance .* squared_diff,'omitnan') / sum(abundance,'omitnan');
            case 'gini'
                % Population is equal to amino acid abundance
                pop = abundance;
                % Get property values for each data table entry
                val = LookupProperties(DT,PT,property);
                val = cell2mat(val);
                % Calculate gini coefficient with 'Income' equal to 
                % property value, constrained to be positive
                feature(k) = gini(pop,abs(val));
            otherwise
                error('unknown metric')
        end
    end
    % Add feature to table
    FT.(sprintf('%s_%s',property,metric))=feature;
end