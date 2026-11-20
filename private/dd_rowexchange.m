function [D, obj, st] = dd_rowexchange(C, D, DR, o, rs, cand)
%DD_ROWEXCHANGE Improve a design by replacing whole choice tasks. Internal.
%   For every row in random order, each candidate task (dd_candidates) that
%   keeps the level frequencies within their bounds and does not repeat
%   another row is tried in its place; the best replacement is made if it
%   lowers the objective. This is a modified Fedorov exchange. Unlike a swap
%   within a column it changes several attributes of a task at once, which
%   matters when constraints tie attributes together, and it changes the
%   level frequencies, so it only has room to work where they are not fixed.
%   Passes repeat until one makes no replacement (at most o.MaxPasses).
%
%   With many prior draws the candidates are first ranked on the draws stored
%   in cand, and the best are evaluated on all draws.

S = C.S;
K = C.K;
M = C.M;
w = C.w;
crit = C.criterion;
act = find(w > 0)';
X = cell(1, M);
Is = cell(1, M);
It = cell(1, M);
fm = zeros(M, 1);
for m = act
    X{m} = dd_X(C.model(m), D);
    Is{m} = dd_rowinfo(C.model(m), X{m}, DR(m));
    It{m} = sum(Is{m}, 3);
    fm(m) = dd_crit(It{m}, C.model(m), crit);
end
obj = sum(w(act) .* fm(act));
cnt = cell(1, K);
for k = 1:K
    cnt{k} = accumarray(D(:, k), 1, [C.L(k) 1]);
end
hw = C.hashW;
h = D * hw;
mw = {};
if C.unique && C.mirror
    for i = 1:numel(C.mirrorPairs)
        pm = 1:K;
        pm([C.mirrorPairs(i).ka C.mirrorPairs(i).kb]) = [C.mirrorPairs(i).kb C.mirrorPairs(i).ka];
        mw{i} = hw(pm); %#ok<AGROW>
    end
end
candM = cell(1, numel(mw));
for i = 1:numel(mw)
    candM{i} = cand.D * mw{i};
end
st = struct('passes', 0, 'evaluations', 0, 'moves', 0, 'timedOut', false);
nCand = size(cand.D, 1);
screen = ~all(cand.full(act));
nTop = 2 * o.ScreenTop;

for pass = 1:o.MaxPasses
    improved = false;
    for s = randperm(rs, S)
        % candidates that keep every column within its frequency bounds
        ok = true(nCand, 1);
        for k = 1:K
            cur = D(s, k);
            ck = cand.D(:, k);
            ok = ok & (ck == cur | (cnt{k}(ck) + 1 <= C.hi{k}(ck) & cnt{k}(cur) - 1 >= C.lo{k}(cur)));
        end
        ok = ok & cand.hash ~= h(s);
        if C.unique
            hOther = h;
            hOther(s) = NaN;
            ok = ok & ~ismember(cand.hash, hOther);
            for i = 1:numel(mw)
                ok = ok & ~ismember(candM{i}, hOther);
            end
        end
        idx = find(ok);
        if isempty(idx)
            continue
        end
        if screen && numel(idx) > nTop
            fs = zeros(numel(idx), 1);
            for m = act
                rsub = cand.rsub{m};
                fs = fs + w(m) * dd_crit(It{m}(:, :, 1, rsub) - Is{m}(:, :, s, rsub) + cand.I{m}(:, :, idx, :), C.model(m), crit);
            end
            st.evaluations = st.evaluations + numel(idx);
            [~, os] = sort(fs);
            idx = idx(os(1:nTop));
        end
        f = zeros(numel(idx), 1);
        In = cell(1, M);
        for m = act
            if cand.full(m)
                In{m} = cand.I{m}(:, :, idx, :);
            else
                In{m} = dd_rowinfo(C.model(m), cand.X{m}(:, :, idx), DR(m));
            end
            f = f + w(m) * dd_crit(It{m} - Is{m}(:, :, s, :) + In{m}, C.model(m), crit);
        end
        st.evaluations = st.evaluations + numel(idx);
        [fb, ib] = min(f);
        if fb < obj * (1 - o.Tol) || (~isfinite(obj) && isfinite(fb))
            c = idx(ib);
            for k = 1:K
                cnt{k}(D(s, k)) = cnt{k}(D(s, k)) - 1;
                cnt{k}(cand.D(c, k)) = cnt{k}(cand.D(c, k)) + 1;
            end
            D(s, :) = cand.D(c, :);
            h(s) = cand.hash(c);
            for m = act
                It{m} = It{m} - Is{m}(:, :, s, :) + In{m}(:, :, ib, :);
                Is{m}(:, :, s, :) = In{m}(:, :, ib, :);
                X{m}(:, :, s) = cand.X{m}(:, :, c);
                fm(m) = dd_crit(It{m}, C.model(m), crit);
            end
            obj = sum(w(act) .* fm(act));
            improved = true;
            st.moves = st.moves + 1;
        end
    end
    st.passes = pass;
    for m = act
        It{m} = sum(Is{m}, 3);
        fm(m) = dd_crit(It{m}, C.model(m), crit);
    end
    obj = sum(w(act) .* fm(act));
    if ~improved
        break
    end
    if posixtime(datetime('now')) > o.deadline
        st.timedOut = true;
        break
    end
end
end
