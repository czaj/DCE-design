function [D, st] = dd_panelpolish(C, D, y, o, scale)
%DD_PANELPOLISH Improve a blocked design on the objective with the panel model itself. Internal.
%   The search of dd_localsearch needs models that add up over choice tasks,
%   so a panel mixed logit enters it through a surrogate. This stage works on
%   the real objective: the panel information is computed per respondent, i.e.
%   per block (C.panelUnit 'block') or over all rows ('design').
%
%   Moves: exchange the levels of two tasks in one column, both tasks in the
%   same block. This keeps the level frequencies of every column in the whole
%   design and inside every block, so the blocking stays as balanced as it was.
%   A move is applied as soon as it lowers the objective; passes repeat until
%   none does, o.PanelPolish passes have run, or o.MaxTime is reached.
%
%   Panel models are evaluated at their prior means in this stage. With
%   scale, the model criteria use the same references as candidate selection.
%   Where a block has too many choice sequences to list, they are simulated:
%   moves are screened on o.PanelRespondents respondents and those that lower
%   the objective are confirmed on o.EvalPanelRespondents others before they
%   are applied (st.confirmed). st.before and st.after are the criteria of the
%   models on the sample that decides (the confirming one, if there is one).

M = C.M;
w = C.w;
act = find(w > 0)';
crit = C.criterion;
Cp = C;
for m = act
    if strcmp(C.model(m).type, 'rppanel')
        Cp.model(m).stat = 'fixed';
    end
end
DR = dd_draws(Cp, o.BayesDraws, o.RandomDraws, o.Seed, o.DrawType, o.PanelRespondents, ...
    o.PanelDraws * strcmp(o.PanelIntegration, 'importance'));
S = size(D, 1);
K = C.K;
blocks = unique(y(:))';
rowsOf = arrayfun(@(b) find(y == b), blocks, 'UniformOutput', false);
if strcmp(C.panelUnit, 'block')
    grpOf = zeros(S, 1);
    for g = 1:numel(blocks)
        grpOf(rowsOf{g}) = g;
    end
    grpRows = rowsOf;
else
    grpOf = ones(S, 1);
    grpRows = {(1:S)'};
end

isPanel = false(1, M);
X = cell(1, M);
Is = cell(1, M);
Ig = cell(1, M);
It = cell(1, M);
f = zeros(M, 1);
for m = act
    Mm = Cp.model(m);
    isPanel(m) = strcmp(Mm.type, 'rppanel');
    X{m} = dd_X(Mm, D);
    if isPanel(m)
        Ig{m} = zeros(Mm.Pt, Mm.Pt, numel(grpRows));
        for g = 1:numel(grpRows)
            Ig{m}(:, :, g) = dd_panelinfo(Mm, X{m}(:, :, grpRows{g}), DR(m), 1, o.PanelExactMax);
        end
        It{m} = sum(Ig{m}, 3);
    else
        Is{m} = dd_rowinfo(Mm, X{m}, DR(m));
        It{m} = sum(Is{m}, 3);
    end
    f(m) = dd_crit(It{m}, Mm, crit);
end
obj = sum(w(act) .* f(act) ./ scale(act));

% Where the choice sequences of a block are simulated, a move that lowers the
% objective on the o.PanelRespondents respondents of the search is confirmed
% on a larger sample (o.EvalPanelRespondents other respondents) before it is
% applied: with few respondents the search would otherwise fit their sample.
Jalt = size(X{act(1)}, 1);
confirm = any(isPanel) && o.EvalPanelRespondents > o.PanelRespondents && ...
    any(cellfun(@(g) Jalt^numel(g) > o.PanelExactMax, grpRows));
fc = f;
if confirm
    DRc = dd_draws(Cp, o.BayesDraws, o.RandomDraws, o.Seed + 4409, o.DrawType, o.EvalPanelRespondents, ...
        o.PanelDraws * strcmp(o.PanelIntegration, 'importance'));
    Igc = cell(1, M);
    Itc = cell(1, M);
    for m = act(isPanel(act))
        Mm = Cp.model(m);
        Igc{m} = zeros(Mm.Pt, Mm.Pt, numel(grpRows));
        for g = 1:numel(grpRows)
            Igc{m}(:, :, g) = dd_panelinfo(Mm, X{m}(:, :, grpRows{g}), DRc(m), 1, o.PanelExactMax);
        end
        Itc{m} = sum(Igc{m}, 3);
        fc(m) = dd_crit(Itc{m}, Mm, crit);
    end
end
objc = sum(w(act) .* fc(act) ./ scale(act));
st = struct('passes', 0, 'moves', 0, 'evaluations', 0, 'before', fc', 'after', fc', 'timedOut', false, 'confirmed', confirm);
if ~all(isfinite(f(act))) || ~all(isfinite(fc(act)))
    return
end
rs = RandStream('mrg32k3a', 'Seed', o.Seed);
rs.Substream = 9001;

for pass = 1:o.PanelPolish
    improved = false;
    for b = randperm(rs, numel(blocks))
        rb = rowsOf{b};
        nb = numel(rb);
        if nb < 2
            continue
        end
        [i1, i2] = find(triu(true(nb), 1));
        for k = randperm(rs, K)
            ord = randperm(rs, numel(i1));
            for q = ord
                s1 = rb(i1(q));
                s2 = rb(i2(q));
                if D(s1, k) == D(s2, k)
                    continue
                end
                Dn = D;
                Dn(s1, k) = D(s2, k);
                Dn(s2, k) = D(s1, k);
                if ~all(dd_rowok(C, Dn([s1 s2], :))) || any(dd_dups(C, Dn))
                    continue
                end
                fn = f;
                Xn = cell(1, M);
                new = cell(1, M);
                Itn = cell(1, M);
                for m = act
                    Mm = Cp.model(m);
                    if isempty(Mm.colTerms{k})
                        continue
                    end
                    Xn{m} = dd_X(Mm, Dn([s1 s2], :));
                    if isPanel(m)
                        g = grpOf(s1);
                        Xg = X{m};
                        Xg(:, :, [s1 s2]) = Xn{m};
                        new{m} = dd_panelinfo(Mm, Xg(:, :, grpRows{g}), DR(m), 1, o.PanelExactMax);
                        Itn{m} = It{m} - Ig{m}(:, :, g) + new{m};
                    else
                        new{m} = dd_rowinfo(Mm, Xn{m}, DR(m));
                        Itn{m} = It{m} - Is{m}(:, :, s1, :) - Is{m}(:, :, s2, :) + sum(new{m}, 3);
                    end
                    fn(m) = dd_crit(Itn{m}, Mm, crit);
                end
                st.evaluations = st.evaluations + 1;
                objn = sum(w(act) .* fn(act) ./ scale(act));
                if objn >= obj * (1 - 1e-9)
                    continue
                end
                fcn = fn;
                objcn = objn;
                if confirm
                    newc = cell(1, M);
                    Itcn = cell(1, M);
                    for m = act(isPanel(act))
                        fcn(m) = fc(m);
                        if isempty(Xn{m})
                            continue
                        end
                        g = grpOf(s1);
                        Xg = X{m};
                        Xg(:, :, [s1 s2]) = Xn{m};
                        newc{m} = dd_panelinfo(Cp.model(m), Xg(:, :, grpRows{g}), DRc(m), 1, o.PanelExactMax);
                        Itcn{m} = Itc{m} - Igc{m}(:, :, g) + newc{m};
                        fcn(m) = dd_crit(Itcn{m}, Cp.model(m), crit);
                    end
                    st.evaluations = st.evaluations + 1;
                    objcn = sum(w(act) .* fcn(act) ./ scale(act));
                    if objcn >= objc * (1 - 1e-9)
                        continue
                    end
                end
                D = Dn;
                for m = act
                    if isempty(Xn{m})
                        continue
                    end
                    X{m}(:, :, [s1 s2]) = Xn{m};
                    if isPanel(m)
                        Ig{m}(:, :, grpOf(s1)) = new{m};
                        if confirm
                            Igc{m}(:, :, grpOf(s1)) = newc{m};
                            Itc{m} = Itcn{m};
                        end
                    else
                        Is{m}(:, :, [s1 s2], :) = new{m};
                    end
                    It{m} = Itn{m};
                end
                f = fn;
                obj = objn;
                fc = fcn;
                objc = objcn;
                improved = true;
                st.moves = st.moves + 1;
            end
            if posixtime(datetime('now')) > o.deadline
                st.timedOut = true;
                break
            end
        end
        if st.timedOut
            break
        end
    end
    st.passes = pass;
    if ~improved || st.timedOut
        break
    end
end
st.after = fc';
end
