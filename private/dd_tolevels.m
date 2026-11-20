function D = dd_tolevels(C, V, where)
%DD_TOLEVELS Convert a design in level values to level indices. Internal.
%   V is a numeric S-by-K matrix in the column order of C.colName, or a table
%   whose variables are named like the columns ('alt.attr', 'alt_attr', or the
%   attribute name for scenario variables). Extra table variables are ignored.

if nargin < 3
    where = 'the design';
end
if istable(V)
    vn = V.Properties.VariableNames;
    M = zeros(height(V), C.K);
    for k = 1:C.K
        want = {C.colName{k}, strrep(C.colName{k}, '.', '_')};
        hit = find(ismember(lower(vn), lower(want)), 1);
        if isempty(hit)
            dd_err('design:columns', '%s has no column ''%s''. Columns expected: %s.', where, C.colName{k}, strjoin(C.colName, ', '));
        end
        M(:, k) = double(V{:, hit});
    end
    V = M;
end
if ~isnumeric(V) || size(V, 2) ~= C.K
    dd_err('design:columns', '%s must have %d columns: %s.', where, C.K, strjoin(C.colName, ', '));
end
D = zeros(size(V));
for k = 1:C.K
    [tf, loc] = ismember(V(:, k), C.levels{k});
    if ~all(tf)
        bad = unique(V(~tf, k))';
        dd_err('design:levels', '%s: column %s contains %s, not among its levels (%s).', where, C.colName{k}, mat2str(bad), mat2str(C.levels{k}));
    end
    D(:, k) = loc;
end
end
