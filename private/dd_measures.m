function out = dd_measures(C, D, y, o)
%DD_MEASURES Full evaluation of a design. Internal.
%   D is S-by-K (level indices), y the block of each row ([] if none).
%   Uses the evaluation draw sets (o.EvalBayesDraws, o.EvalRandomDraws,
%   o.EvalPanelDraws, o.EvalPanelRespondents).

S = size(D, 1);
out = struct();

% ---- per model ----------------------------------------------------------------
Cf = C;
for m = 1:C.M
    Cf.model(m).stat = 'fixed';
end
nIS = o.EvalPanelDraws * strcmp(o.PanelIntegration, 'importance');
DRf = dd_draws(Cf, 1, o.EvalRandomDraws, o.Seed + 1, o.DrawType, o.EvalPanelRespondents, nIS);
% prior draws per model: mixed logit models are costly per draw
nB = repmat(o.EvalBayesDraws, 1, C.M);
for m = 1:C.M
    if strcmp(C.model(m).type, 'rppanel')
        nB(m) = min(nB(m), o.PanelBayesDraws);
    elseif strcmp(C.model(m).type, 'rp')
        nB(m) = min(nB(m), 4 * o.PanelBayesDraws);
    end
end
DRe = dd_draws(C, nB, o.EvalRandomDraws, o.Seed + 1, o.DrawType, o.EvalPanelRespondents, nIS);
obj = 0;
for m = 1:C.M
    Mm = C.model(m);
    r = struct();
    r.name = Mm.name;
    r.type = Mm.type;
    r.parameters = Mm.thNames(:);
    r.prior = [Mm.b0; Mm.s0(Mm.ridx)];
    r.inCriterion = Mm.eff(:);
    % at the prior means
    [If, how, pdg] = dd_designinfo(Cf, m, D, y, DRf(m), o);
    r.panelIntegration = pdg;
    Mf = Cf.model(m);
    r.information = If;
    r.method = how;
    r.DerrorFixed = dd_crit(If, Mf, 'd');
    r.AerrorFixed = dd_crit(If, Mf, 'a');
    ld = dd_logdet(If);
    if isfinite(ld)
        Om = inv(If);
        r.avc = Om;
        r.se = sqrt(max(diag(Om), 0));
        r.tRatio = abs(r.prior) ./ r.se;
        r.Sestimate = (1.96 * r.se ./ abs(r.prior)).^2;
        r.Sestimate(r.prior == 0) = NaN;
    else
        r.avc = nan(Mm.Pt);
        r.se = nan(Mm.Pt, 1);
        r.tRatio = nan(Mm.Pt, 1);
        r.Sestimate = nan(Mm.Pt, 1);
    end
    r.identified = isfinite(ld);
    % over the prior draws
    r.stat = Mm.stat;
    if strcmp(Mm.stat, 'fixed')
        r.Derror = r.DerrorFixed;
        r.Aerror = r.AerrorFixed;
        r.bayes = [];
    else
        rsub = 1:DRe(m).R;
        if strcmp(Mm.type, 'rppanel')
            rsub = 1:min(DRe(m).R, o.PanelBayesDraws);
        end
        Ie = dd_designinfo(C, m, D, y, DRe(m), o, rsub);
        [fd, frd] = dd_crit(Ie, Mm, 'd');
        [fa, ~] = dd_crit(Ie, Mm, 'a');
        r.Derror = fd;
        r.Aerror = fa;
        fin = frd(isfinite(frd));
        r.bayes = struct('draws', numel(rsub), 'mean', mean(frd), 'median', median(frd), ...
            'sd', std(fin), 'min', min(frd), 'max', max(frd), 'q05', quant(fin, 0.05), 'q95', quant(fin, 0.95), ...
            'shareNotIdentified', mean(~isfinite(frd)));
    end
    if strcmp(C.criterion, 'd')
        r.criterion = r.Derror;
    else
        r.criterion = r.Aerror;
    end
    if C.w(m) > 0
        obj = obj + C.w(m) * r.criterion;
    end
    if m == 1
        out.models = r;
    else
        out.models(m) = r;
    end
end
out.objective = obj;

% ---- choice probabilities at the prior means of model 1 ---------------------------
Mm = C.model(1);
X = dd_X(Mm, D);
if strcmp(Mm.type, 'mnl')
    Bm = Mm.b0;
else
    Bm = DRf(1).Beta{1};
end
V = pagemtimes(X, reshape(Bm, Mm.P, 1, 1, []));
V = V - max(V, [], 1);
p = exp(V);
p = mean(p ./ sum(p, 1), 4);
P = reshape(p, C.J, S)';
out.probabilities = P;
out.utilityBalance = 100 * mean(prod(P * C.J, 2));
out.maxProbability = max(P, [], 2);
out.meanShare = mean(P, 1);

% ---- design properties -----------------------------------------------------------------
Lmax = max(C.L);
cnt = nan(Lmax, C.K);
for k = 1:C.K
    cnt(1:C.L(k), k) = accumarray(D(:, k), 1, [C.L(k) 1]);
end
out.levelCounts = cnt;
Vl = zeros(S, C.K);
for k = 1:C.K
    Vl(:, k) = C.LV{k}(D(:, k));
end
out.values = Vl;
sd = std(Vl, 0, 1);
cr = nan(C.K);
if all(sd > 0)
    cr = corrcoef(Vl);
end
out.correlations = cr;

% ---- checks ------------------------------------------------------------------------------
[ok, why] = dd_rowok(C, D);
dup = dd_dups(C, D);
freqOK = true;
for k = 1:C.K
    c = accumarray(D(:, k), 1, [C.L(k) 1]);
    freqOK = freqOK && all(c >= C.lo{k}) && all(c <= C.hi{k});
end
rules = {'reject expression', 'require expression', 'dominated alternative', 'identical alternatives', 'maxProb'};
out.checks = struct('constraintsOK', all(ok), 'violations', {[find(~ok) why(~ok)]}, ...
    'violationRules', {rules}, 'duplicateTasks', find(dup)', 'frequenciesOK', freqOK, ...
    'identified', all([out.models.identified]));
end

function q = quant(x, p)
if isempty(x)
    q = NaN;
    return
end
x = sort(x(:));
n = numel(x);
h = (n - 1) * p + 1;
q = x(floor(h)) + (h - floor(h)) * (x(min(n, floor(h) + 1)) - x(floor(h)));
end
