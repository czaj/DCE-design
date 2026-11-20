function [D, obj, st] = dd_localsearch(C, D, DR, o, rs)
%DD_LOCALSEARCH Improve a design by swaps within columns and level changes. Internal.
%   Moves: exchange the levels of two tasks in one column (keeps the level
%   frequencies), or change one level where the frequency bounds leave room.
%   Columns are visited in random order; within a column the rows are taken in
%   chunks. For a chunk, every move that involves one of its rows is evaluated
%   in one vectorised step. The best move is applied if it lowers the
%   objective; the remaining moves of the chunk that touch other rows are then
%   re-evaluated against the updated information matrix (their own row
%   contributions are unchanged) and applied while they keep improving.
%   Passes repeat until one finds no improving move, which makes the design a
%   local optimum with respect to all such moves.
%
%   The models are additive over tasks (mnl, rp, and the surrogate of rppanel),
%   so the information matrix is updated from the rows that change.
%   obj = sum_m w_m * criterion_m.

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
Rmax = 1;
for m = act
    X{m} = dd_X(C.model(m), D);
    Is{m} = dd_rowinfo(C.model(m), X{m}, DR(m));
    It{m} = sum(Is{m}, 3);
    fm(m) = dd_crit(It{m}, C.model(m), crit);
    Rmax = max(Rmax, DR(m).R);
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
st = struct('passes', 0, 'evaluations', 0, 'moves', 0, 'timedOut', false);
pageTarget = 6000;
Xn = cell(1, M);
In = cell(1, M);
dI = cell(1, M);
allRows = (1:S)';
% With many prior draws, the moves of a chunk are first ranked on a few draws
% and only the best are evaluated on all of them (a move is applied only if it
% improves the objective on all draws). When a pass finds nothing, a wider
% screen (four times the draws, twice the moves) is used before the search
% stops; with o.FullCheck a last pass evaluates every move on all draws.
canScreen = Rmax >= 4 * o.ScreenDraws;
level = 1;
if ~canScreen
    level = 3;
end

for pass = 1:o.MaxPasses
    improved = false;
    switch level
        case 1
            sDraws = o.ScreenDraws;
            sTop = o.ScreenTop;
        case 2
            sDraws = min(Rmax, 4 * o.ScreenDraws);
            sTop = 2 * o.ScreenTop;
        otherwise
            sDraws = Rmax;
            sTop = Inf;
    end
    screenOn = level < 3 && sDraws < Rmax;
    Reff = sDraws;
    for k = randperm(rs, K)
        Lk = C.L(k);
        chunk = max(1, min(S, floor(pageTarget / (Reff * (S * (1 - 1 / Lk) + Lk)))));
        rowOrder = randperm(rs, S);
        for c0 = 1:chunk:S
            rows = rowOrder(c0:min(S, c0 + chunk - 1))';
            col = D(:, k);
            % ---- rows with one level replaced that the moves of this chunk need
            need = false(S, Lk);
            need(rows, :) = true;
            need(:, unique(col(rows))) = true;
            need(allRows + (col - 1) * S) = false;
            [sI, lI] = find(need);
            n1 = numel(sI);
            if n1 == 0
                continue
            end
            Ds = D(sI, :);
            Ds(:, k) = lI;
            ok1 = dd_rowok(C, Ds);
            hN = Ds * hw;
            hM = cell(1, numel(mw));
            if C.unique
                ok1 = ok1 & ~ismember(hN, h);
                for i = 1:numel(mw)
                    hM{i} = Ds * mw{i};
                    ok1 = ok1 & ~ismember(hM{i}, h);
                end
            end
            pos = zeros(S, Lk);
            pos(sI + (lI - 1) * S) = (1:n1)';
            % ---- swap moves: (sa in chunk, sb any row with another level)
            inChunk = false(S, 1);
            inChunk(rows) = true;
            SA = repmat(rows, 1, S);
            SB = repmat(allRows', numel(rows), 1);
            keep = reshape(col(SA), size(SA)) ~= reshape(col(SB), size(SB)) & (~reshape(inChunk(SB), size(SB)) | SB > SA);
            sa = reshape(SA(keep), [], 1);
            sb = reshape(SB(keep), [], 1);
            pa = pos(sa + (col(sb) - 1) * S);
            pb = pos(sb + (col(sa) - 1) * S);
            v = ok1(pa) & ok1(pb);
            for i = 1:numel(mw)
                v = v & hM{i}(pa) ~= hN(pb) & hM{i}(pb) ~= hN(pa);
            end
            sa = sa(v); sb = sb(v); pa = pa(v); pb = pb(v);
            nS = numel(sa);
            % ---- level changes allowed by the frequency bounds
            ca = zeros(0, 1); cl = zeros(0, 1); pc = zeros(0, 1);
            canLeave = cnt{k}(col(rows)) - 1 >= C.lo{k}(col(rows));
            canEnter = (cnt{k} + 1 <= C.hi{k})';
            if any(canLeave) && any(canEnter)
                CA = repmat(rows, 1, Lk);
                CL = repmat(1:Lk, numel(rows), 1);
                keep = canLeave & canEnter & CL ~= reshape(col(CA), size(CA));
                ca = reshape(CA(keep), [], 1);
                cl = reshape(CL(keep), [], 1);
                pc = pos(ca + (cl - 1) * S);
                v = ok1(pc);
                ca = ca(v); cl = cl(v); pc = pc(v);
            end
            nC = numel(ca);
            nCand = nS + nC;
            if nCand == 0
                continue
            end
            % ---- information of the replaced rows
            for m = act
                Mm = C.model(m);
                tl = Mm.colTerms{k};
                if isempty(tl)
                    In{m} = [];
                    continue
                end
                Xn{m} = X{m}(:, :, sI) - dd_X(Mm, D(sI, :), tl) + dd_X(Mm, Ds, tl);
                In{m} = dd_rowinfo(Mm, Xn{m}, DR(m));
            end
            % ---- screening on a few prior draws
            if screenOn && nCand > sTop
                fs = zeros(nCand, 1);
                for m = act
                    if isempty(In{m}), continue; end
                    rsub = 1:min(sDraws, DR(m).R);
                    dIs = cat(3, In{m}(:, :, pa, rsub) + In{m}(:, :, pb, rsub) - Is{m}(:, :, sa, rsub) - Is{m}(:, :, sb, rsub), ...
                        In{m}(:, :, pc, rsub) - Is{m}(:, :, ca, rsub));
                    fs = fs + w(m) * dd_crit(It{m}(:, :, 1, rsub) + dIs, C.model(m), crit);
                end
                st.evaluations = st.evaluations + nCand;
                [~, os] = sort(fs);
                top = sort(os(1:sTop));
                selS = top(top <= nS);
                selC = top(top > nS) - nS;
                sa = sa(selS); sb = sb(selS); pa = pa(selS); pb = pb(selS);
                ca = ca(selC); cl = cl(selC); pc = pc(selC);
                nS = numel(sa);
                nC = numel(ca);
                nCand = nS + nC;
            end
            % ---- objective of every (remaining) move on all prior draws
            f = zeros(nCand, 1);
            for m = act
                if isempty(In{m})
                    dI{m} = [];
                    f = f + w(m) * fm(m);
                    continue
                end
                dI{m} = cat(3, In{m}(:, :, pa, :) + In{m}(:, :, pb, :) - Is{m}(:, :, sa, :) - Is{m}(:, :, sb, :), ...
                    In{m}(:, :, pc, :) - Is{m}(:, :, ca, :));
                f = f + w(m) * dd_crit(It{m} + dI{m}, C.model(m), crit);
            end
            st.evaluations = st.evaluations + nCand;
            % ---- apply the best move; re-evaluate the others while they improve
            alive = true(nCand, 1);
            rowA = [sa; ca];
            rowB = [sb; zeros(nC, 1)];
            posA = [pa; pc];
            posB = [pb; zeros(nC, 1)];
            while true
                [fb, ib] = min(f);
                if ~(fb < obj * (1 - o.Tol) || (~isfinite(obj) && isfinite(fb)))
                    break
                end
                s1 = rowA(ib);
                s2 = rowB(ib);
                p1 = posA(ib);
                changed = s1;
                if s2 > 0
                    p2 = posB(ib);
                    l1 = D(s2, k);
                    l2 = D(s1, k);
                    D(s1, k) = l1;
                    D(s2, k) = l2;
                    changed = [s1; s2];
                else
                    l1 = cl(ib - nS);
                    cnt{k}(D(s1, k)) = cnt{k}(D(s1, k)) - 1;
                    cnt{k}(l1) = cnt{k}(l1) + 1;
                    D(s1, k) = l1;
                end
                for m = act
                    if isempty(dI{m}), continue; end
                    X{m}(:, :, s1) = Xn{m}(:, :, p1);
                    Is{m}(:, :, s1, :) = In{m}(:, :, p1, :);
                    if s2 > 0
                        X{m}(:, :, s2) = Xn{m}(:, :, p2);
                        Is{m}(:, :, s2, :) = In{m}(:, :, p2, :);
                    end
                    It{m} = It{m} + dI{m}(:, :, ib, :);
                    fm(m) = dd_crit(It{m}, C.model(m), crit);
                end
                obj = sum(w(act) .* fm(act));
                h(changed) = D(changed, :) * hw;
                improved = true;
                st.moves = st.moves + 1;
                % moves that touch a changed row are no longer valid
                alive = alive & ~ismember(rowA, changed) & ~ismember(rowB, changed);
                if s2 == 0
                    alive(nS + 1:end) = false;
                end
                if C.unique
                    hc = h(changed);
                    dupNow = ismember(hN(posA), hc);
                    dupNow(1:nS) = dupNow(1:nS) | ismember(hN(posB(1:nS)), hc);
                    for i = 1:numel(mw)
                        dupNow = dupNow | ismember(hM{i}(posA), hc);
                        dupNow(1:nS) = dupNow(1:nS) | ismember(hM{i}(posB(1:nS)), hc);
                    end
                    alive = alive & ~dupNow;
                end
                idx = find(alive);
                f(:) = Inf;
                if isempty(idx)
                    break
                end
                fn = zeros(numel(idx), 1);
                for m = act
                    if isempty(dI{m})
                        fn = fn + w(m) * fm(m);
                    else
                        fn = fn + w(m) * dd_crit(It{m} + dI{m}(:, :, idx, :), C.model(m), crit);
                    end
                end
                f(idx) = fn;
                st.evaluations = st.evaluations + numel(idx);
            end
        end
    end
    st.passes = pass;
    % remove accumulated rounding from the running totals
    for m = act
        It{m} = sum(Is{m}, 3);
        fm(m) = dd_crit(It{m}, C.model(m), crit);
    end
    obj = sum(w(act) .* fm(act));
    if ~improved
        if level == 1 || (level == 2 && o.FullCheck)
            level = level + 1;
            continue
        end
        break
    elseif canScreen
        level = 1;
    end
    if posixtime(datetime('now')) > o.deadline
        st.timedOut = true;
        break
    end
end
end
