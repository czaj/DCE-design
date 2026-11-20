function dd_precheck(C)
%DD_PRECHECK Stop early when the constraints leave too few feasible choice tasks. Internal.
%   Raises dce:search:infeasible with a diagnosis when (a) no feasible task
%   is found among random tasks, or (b) the feasible tasks can be enumerated
%   and there are fewer distinct ones than the design has rows.

rs = RandStream('mt19937ar', 'Seed', 4711);
n = 4000;
K = C.K;
Dr = zeros(n, K);
for k = 1:K
    Dr(:, k) = randi(rs, C.L(k), n, 1);
end
[ok, why] = dd_rowok(C, Dr);
rule = {'a reject expression', 'a require expression', 'dominance', 'identical alternatives', 'maxProb'};
share = arrayfun(@(i) mean(why == i), 1:5);
txt = strjoin(arrayfun(@(i) sprintf('%s %.0f%%', rule{i}, 100 * share(i)), find(share > 0), 'UniformOutput', false), ', ');
nFull = prod(C.L);
if nFull <= 2e5
    g = cell(1, K);
    for k = 1:K
        g{k} = 1:C.L(k);
    end
    [g{:}] = ndgrid(g{:});
    Dall = cell2mat(cellfun(@(v) v(:), g, 'UniformOutput', false));
    Dall = Dall(dd_rowok(C, Dall), :);
    nFeasible = size(Dall, 1);
    if C.unique
        nFeasible = nFeasible - nnz(dd_dups(C, Dall));
    end
    if nFeasible == 0 || (C.unique && nFeasible < C.S)
        mirrorNote = '';
        if C.unique && C.mirror && ~isempty(C.mirrorPairs)
            mirrorNote = ' (mirror images counted once)';
        end
        dd_err('search:infeasible', ['Only %d distinct choice tasks satisfy the constraints%s, but the design needs %d rows ' ...
            '(random tasks are excluded by: %s). Relax the constraints (spec.reject, spec.require, spec.dominance), ' ...
            'add levels, reduce spec.rows, or allow repeated tasks with spec.unique = false.'], nFeasible, mirrorNote, C.S, txt);
    end
elseif ~any(ok)
    dd_err('search:infeasible', ['None of %d random choice tasks satisfies the constraints (excluded by: %s). ' ...
        'Check spec.reject, spec.require and spec.dominance.'], n, txt);
end
end
