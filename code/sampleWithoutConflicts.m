% sampleWithoutConflicts.m: Weighted sampling of amino acids without replacement,
% treating an amino acid's unresolved and resolved entries as mutually exclusive
%
% Draws n entries one at a time with probability proportional to their
% weights. After each draw, the drawn entry is removed from the pool, along
% with any entry that would represent the same amino acid twice: drawing an
% unresolved ('Both') entry removes the D- and L- entries of that amino
% acid, and drawing a D- or L- entry removes its 'Both' entry. D- and L-
% entries of the same amino acid can both be drawn.
%
% Syntax: idx = sampleWithoutConflicts(w,grp,isBoth,n)
%
% Inputs:
%   w        Nonnegative selection weight of each entry
%   grp      Group index of each entry; entries of the same amino acid share a group
%   isBoth   Logical, true for unresolved ('Both') entries
%   n        Number of entries to draw
%
% Output:
%   idx: n-by-1 indices of the drawn entries

function idx = sampleWithoutConflicts(w,grp,isBoth,n)
    w = w(:); grp = grp(:); isBoth = isBoth(:);
    idx = zeros(n,1);
    for k = 1:n
        c = cumsum(w);
        if c(end) <= 0
            error('sampleWithoutConflicts:exhausted', ...
                'Only %d of the %d requested entries could be drawn.', k-1, n);
        end
        % Draw one entry with probability proportional to its weight
        i = find(rand*c(end) < c, 1);
        idx(k) = i;
        % Remove the drawn entry and any conflicting entries of the same amino acid
        w(i) = 0;
        if isBoth(i)
            w(grp == grp(i)) = 0;
        else
            w(grp == grp(i) & isBoth) = 0;
        end
    end
end
