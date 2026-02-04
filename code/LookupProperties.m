% LookupProperties.m: Retrieves property values for molecules in a data table
%
% Maps molecular identifiers from a data table to their corresponding property 
% values in a property table using case-insensitive matching.
%
% Syntax: val = LookupProperties(DT,PT,property,field)
%
% Inputs:
%   DT       Data table containing molecular identifiers
%   PT       Property table containing molecular properties
%   property Name of the property to retrieve (string)
%   field    Field name for matching identifiers (default: 'Symbol')
%
% Output:
%   val: Cell array of property values corresponding to each entry in DT

function val = LookupProperties(DT,PT,property,field)
    if nargin<4, field = 'Symbol'; end % Support any text field
    % Get property values for each data table entry
    val = cellfun(@(x)(PT.(property)(find(strcmpi(PT.(field),x)))),DT.(field),'UniformOutput', false);
end