function dup = dd_dups(C, D)
%DD_DUPS Flag choice tasks that repeat an earlier task. Internal.
%   With C.mirror, a task that equals an earlier task after exchanging two
%   interchangeable alternatives also counts as a repeat.

S = size(D, 1);
dup = false(S, 1);
if ~C.unique
    return
end
[~, ia] = unique(D, 'rows', 'stable');
dup(setdiff(1:S, ia)) = true;
if C.mirror
    for i = 1:numel(C.mirrorPairs)
        Dm = D;
        Dm(:, [C.mirrorPairs(i).ka C.mirrorPairs(i).kb]) = D(:, [C.mirrorPairs(i).kb C.mirrorPairs(i).ka]);
        [tf, loc] = ismember(Dm, D, 'rows');
        dup = dup | (tf & loc < (1:S)');
    end
end
end
