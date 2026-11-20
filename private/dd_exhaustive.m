function [D, f, nTried, nValid] = dd_exhaustive(C, DR, o, pool)
%DD_EXHAUSTIVE Evaluate every design and return the best. Internal.
%   A design is a set of C.S distinct feasible choice tasks whose level
%   frequencies are within their bounds. Possible only for small problems:
%   all feasible tasks must be in the candidate pool and the number of sets
%   must not exceed o.ExhaustiveMax. The result is the global optimum of the
%   search objective (for Bayesian priors: on the search draws).

S = C.S;
nF = size(pool.D, 1);
act = find(C.w > 0)';
total = exp(gammaln(nF + 1) - gammaln(S + 1) - gammaln(max(nF - S, 0) + 1));
if ~pool.all || ~all(pool.full(act))
    dd_err('search:exhaustive', ['Method ''exhaustive'' needs every feasible choice task in the candidate pool, with its information for all ' ...
        'prior draws; this problem has too many. Use the default method.']);
end
if nF < S
    dd_err('search:exhaustive', 'Only %d feasible choice tasks exist for %d rows.', nF, S);
end
if total > o.ExhaustiveMax
    dd_err('search:exhaustive', ['Method ''exhaustive'' would have to evaluate about 10^%.0f designs (%d feasible tasks, %d rows); ' ...
        'the limit is %g (option ExhaustiveMax). Use the default method.'], log10(total), nF, S, o.ExhaustiveMax);
end
combos = nchoosek(1:nF, S);
nTried = size(combos, 1);
% level frequencies
ok = true(nTried, 1);
for k = 1:C.K
    lev = reshape(pool.D(combos, k), nTried, S);
    for l = 1:C.L(k)
        c = sum(lev == l, 2);
        ok = ok & c >= C.lo{k}(l) & c <= C.hi{k}(l);
    end
end
% mirror images of a task inside the same design
if C.unique && C.mirror
    for i = 1:numel(C.mirrorPairs)
        pm = 1:C.K;
        pm([C.mirrorPairs(i).ka C.mirrorPairs(i).kb]) = [C.mirrorPairs(i).kb C.mirrorPairs(i).ka];
        hm = pool.D * C.hashW(pm);
        H = reshape(pool.hash(combos), nTried, S);
        Hm = reshape(hm(combos), nTried, S);
        clash = false(nTried, 1);
        for a = 1:S
            clash = clash | any(Hm(:, a) == H(:, [1:a - 1, a + 1:S]), 2);
        end
        ok = ok & ~clash;
    end
end
combos = combos(ok, :);
nValid = size(combos, 1);
if nValid == 0
    dd_err('search:infeasible', 'No set of %d feasible tasks keeps the level frequencies within their bounds.', S);
end
f = Inf;
best = 1;
chunk = 20000;
for c0 = 1:chunk:nValid
    ii = c0:min(nValid, c0 + chunk - 1);
    fc = zeros(numel(ii), 1);
    for m = act
        Mm = C.model(m);
        It = zeros(Mm.Pt, Mm.Pt, numel(ii), DR(m).R);
        for t = 1:S
            It = It + pool.I{m}(:, :, combos(ii, t), :);
        end
        fc = fc + C.w(m) * dd_crit(It, Mm, C.criterion);
    end
    [fb, ib] = min(fc);
    if fb < f
        f = fb;
        best = ii(ib);
    end
    if posixtime(datetime('now')) > o.deadline
        dd_err('search:exhaustive', 'MaxTime was reached after %d of %d designs; an exhaustive search cannot be cut short.', ii(end), nValid);
    end
end
if ~isfinite(f)
    dd_err('search:notIdentified', 'None of the %d designs identifies all parameters.', nValid);
end
D = pool.D(combos(best, :), :);
end
