function [A, isOA, maxDev] = dd_oa(S, Lv, tries, rs)
%DD_OA Orthogonal array of strength 2 with mixed levels, by search. Internal.
%   A is S-by-K with entries 1..Lv(k), every level equally often in a column.
%   The array is built column by column with the J2 algorithm of Xu (2002,
%   Technometrics 44, 356-368): each new column starts as a random balanced
%   column and pairs of its entries are exchanged while that lowers
%   J2 = sum_{i<j} delta_ij^2, where delta_ij counts (with weights Lv) the
%   columns in which rows i and j coincide. J2 reaches its lower bound exactly
%   when the array is orthogonal.
%   isOA is true if every pair of columns shows every pair of levels equally
%   often. Otherwise A is the most nearly orthogonal array found in 'tries'
%   attempts, and maxDev the largest deviation of a pair count from its
%   expected value, relative to that value.

K = numel(Lv);
if any(mod(S, Lv) ~= 0)
    error('dce:orthogonal:rows', 'The number of rows (%d) must be a multiple of every number of levels.', S);
end
[~, order] = sort(Lv, 'descend');
Ls = Lv(order);
best = [];
bestDev = Inf;
for t = 1:tries
    A = zeros(S, K);
    Delta = zeros(S);
    for c = 1:K
        L = Ls(c);
        base = repelem((1:L)', S / L);
        bestCol = [];
        bestG = Inf;
        % the smallest possible g for this column is reached when it is
        % orthogonal to all earlier ones: then every earlier column contributes
        % w * (number of pairs coinciding in both) = w * S * (S/(L*Lk) - 1) / 2
        gMin = 0;
        for q = 1:c - 1
            gMin = gMin + Ls(q) * S * (S / (L * Ls(q)) - 1) / 2;
        end
        for t1 = 1:max(10, 3 * tries)
            col = base(randperm(rs, S));
            [col, g] = improveColumn(col, Delta, L);
            if g < bestG
                bestG = g;
                bestCol = col;
            end
            if bestG <= gMin + 1e-9
                break
            end
        end
        A(:, c) = bestCol;
        Delta = Delta + L * (bestCol == bestCol');
        Delta(1:S + 1:end) = 0;
    end
    dev = orthDeviation(A, Ls);
    if dev < bestDev
        bestDev = dev;
        best = A;
    end
    if bestDev < 1e-12
        break
    end
end
A = zeros(S, K);
A(:, order) = best;
maxDev = bestDev;
isOA = bestDev < 1e-12;
end

function [col, g] = improveColumn(col, Delta, L)
S = numel(col);
for it = 1:10 * S
    H = double(col == (1:L));
    Sm = Delta * H;
    own = Sm((1:S)' + (col - 1) * S);
    M1 = Sm(:, col);
    dg = M1 + M1' - 2 * Delta - own - own';
    dg(col == col') = 0;
    [m, ix] = min(dg(:));
    if m >= -1e-9
        break
    end
    [a, b] = ind2sub([S S], ix);
    tmp = col(a);
    col(a) = col(b);
    col(b) = tmp;
end
same = col == col';
g = sum(Delta(same), 'all') / 2;
end

function dev = orthDeviation(A, Lv)
[S, K] = size(A);
dev = 0;
for a = 1:K - 1
    for b = a + 1:K
        tab = accumarray([A(:, a) A(:, b)], 1, [Lv(a) Lv(b)]);
        e = S / (Lv(a) * Lv(b));
        dev = max(dev, max(abs(tab(:) - e)) / e);
    end
end
end
