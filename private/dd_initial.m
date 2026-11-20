function D = dd_initial(C, rs)
%DD_INITIAL Random feasible starting design. Internal.
%   Columns are drawn within the level-frequency bounds, then rows that break
%   a constraint or repeat another row are repaired by swaps within columns
%   (and by level changes where the bounds leave room). A repair move is kept
%   when it does not increase the number of offending rows.

S = C.S;
K = C.K;
for attempt = 1:8
    D = zeros(S, K);
    for k = 1:K
        cnt = C.lo{k};
        left = S - sum(cnt);
        while left > 0
            room = find(cnt < C.hi{k});
            l = room(randi(rs, numel(room)));
            cnt(l) = cnt(l) + 1;
            left = left - 1;
        end
        col = repelem((1:C.L(k))', cnt);
        D(:, k) = col(randperm(rs, S));
    end
    okRow = dd_rowok(C, D);
    dup = dd_dups(C, D);
    bad = ~okRow | dup;
    cnt = cell(1, K);
    for k = 1:K
        cnt{k} = accumarray(D(:, k), 1, [C.L(k) 1]);
    end
    it = 0;
    lastGain = 0;
    maxIt = 400 * S * K;
    while any(bad) && it < maxIt && it - lastGain < 60 * S * K
        it = it + 1;
        ib = find(bad);
        s = ib(randi(rs, numel(ib)));
        k = randi(rs, K);
        cur = D(s, k);
        free = find(cnt{k} + 1 <= C.hi{k} & (1:C.L(k))' ~= cur);
        if ~isempty(free) && cnt{k}(cur) - 1 >= C.lo{k}(cur) && rand(rs) < 0.5
            l = free(randi(rs, numel(free)));
            rows = s;
            newRows = D(s, :);
            newRows(k) = l;
            changeLevel = true;
        else
            cand = find(D(:, k) ~= cur);
            if isempty(cand)
                continue
            end
            s2 = cand(randi(rs, numel(cand)));
            rows = [s; s2];
            newRows = D(rows, :);
            newRows(:, k) = D([s2; s], k);
            changeLevel = false;
        end
        okNew = dd_rowok(C, newRows);
        % cheap screen on the changed rows; duplicates are recounted only when it passes
        if nnz(~okNew) > nnz(~okRow(rows))
            continue
        end
        Dn = D;
        Dn(rows, :) = newRows;
        okN = okRow;
        okN(rows) = okNew;
        dupN = dd_dups(C, Dn);
        badN = ~okN | dupN;
        if nnz(badN) <= nnz(bad)
            if nnz(badN) < nnz(bad)
                lastGain = it;
            end
            D = Dn;
            okRow = okN;
            bad = badN;
            if changeLevel
                cnt{k}(cur) = cnt{k}(cur) - 1;
                cnt{k}(l) = cnt{k}(l) + 1;
            end
        end
    end
    if ~any(bad)
        return
    end
end
n = 4000;
Dr = zeros(n, K);
for k = 1:K
    Dr(:, k) = randi(rs, C.L(k), n, 1);
end
[ok, why] = dd_rowok(C, Dr);
rule = {'a reject expression', 'a require expression', 'dominance', 'identical alternatives', 'maxProb'};
share = arrayfun(@(i) mean(why == i), 1:5);
txt = strjoin(arrayfun(@(i) sprintf('%s %.0f%%', rule{i}, 100 * share(i)), find(share > 0), 'UniformOutput', false), ', ');
dd_err('search:infeasible', ['No design was found that satisfies the constraints together with the level frequencies. ' ...
    'Of random choice tasks, %.1f%% are feasible (excluded by: %s). Relax the constraints, use spec.balance = ''none'' ' ...
    'or spec.frequency bounds, or increase the number of rows.'], 100 * mean(ok), txt);
end
