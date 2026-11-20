function adv = dd_advice(C, result)
%DD_ADVICE Remarks on a specification and, when given, on its design. Internal.
%   adv is a struct array with fields
%     topic  'size', 'rows', 'blocks', 'priors', 'constraints', 'parameters',
%            'search' or 'efficiency'
%     level  'info' (a fact worth knowing), 'suggestion' (something that could
%            be done better) or 'warning' (something that is likely wrong)
%     text   one or two sentences, with the numbers they rest on
%   Everything is computed from the specification and the design; nothing is
%   a general rule of thumb without a number behind it.

if nargin < 2
    result = [];
end
adv = struct('topic', {}, 'level', {}, 'text', {});
S = C.S;
J = C.J;
K = C.K;
B = C.B;
M1 = C.model(1);
Pmax = max([C.model.Pt]);

% ---- size of the problem ------------------------------------------------------
nTasks = prod(C.L);
rs = RandStream('mt19937ar', 'Seed', 4711);
enumerable = nTasks <= 2e5;
if enumerable
    g = cell(1, K);
    for k = 1:K
        g{k} = 1:C.L(k);
    end
    [g{:}] = ndgrid(g{:});
    Dall = cell2mat(cellfun(@(v) v(:), g, 'UniformOutput', false));
    ok = dd_rowok(C, Dall);
    feasShare = mean(ok);
    Dfeas = Dall(ok, :);
    nFeas = size(Dfeas, 1);
    if C.unique
        nFeas = nFeas - nnz(dd_dups(C, Dfeas));
    end
else
    Dr = zeros(4000, K);
    for k = 1:K
        Dr(:, k) = randi(rs, C.L(k), 4000, 1);
    end
    feasShare = mean(dd_rowok(C, Dr));
    nFeas = feasShare * nTasks;
    Dfeas = [];
end
% number of designs: sets of S feasible tasks (an upper bound when level
% frequencies are restricted)
log10designs = (gammaln(nFeas + 1) - gammaln(S + 1) - gammaln(max(nFeas - S, 0) + 1)) / log(10);
rate = evaluationRate(C);
years = 10^(log10designs - log10(rate) - log10(3600 * 24 * 365));
if log10designs <= log10(2e6) && nFeas >= S
    adv = add(adv, 'size', 'suggestion', sprintf(['There are %s feasible choice tasks and about %s possible designs. That is few enough to evaluate ' ...
        'every one (about %s): dce_design(spec, ''Method'', ''exhaustive'') returns the global optimum.'], ...
        big(nFeas), big(10^log10designs), duration(10^log10designs / rate)));
else
    adv = add(adv, 'size', 'info', sprintf(['%s possible choice tasks, %s of them feasible (%.0f%%); about 10^%.0f possible designs. ' ...
        'Evaluating every design at about %s per second would take %s, so the search improves designs step by step from several random starts.'], ...
        big(nTasks), big(nFeas), 100 * feasShare, log10designs, big(rate), longDuration(years)));
end

% ---- rows ------------------------------------------------------------------------
obs = S * (J - 1);
if obs < 2 * Pmax
    adv = add(adv, 'rows', 'suggestion', sprintf(['%d rows with %d alternatives give %d independent choice observations for %d parameters ' ...
        '(ratio %.1f). With a ratio below 2 the design depends strongly on the priors; %d rows or more would give a ratio of 2.'], ...
        S, J, obs, Pmax, obs / Pmax, ceil(2 * Pmax / (J - 1))));
end
Lu = unique(C.L);
lcmAll = 1;
for l = Lu
    lcmAll = lcm(lcmAll, l);
end
strictCols = cellfun(@(a, b) all(b - a <= 1), C.lo, C.hi);
unbal = find(mod(S, C.L) ~= 0 & strictCols);
if ~isempty(unbal)
    lowS = floor(S / lcmAll) * lcmAll;
    highS = ceil(S / lcmAll) * lcmAll;
    opts = unique([lowS highS]);
    opts = opts(opts >= ceil(Pmax / (J - 1)) & opts > 0);
    adv = add(adv, 'rows', 'suggestion', sprintf(['%d rows cannot balance attributes with %s levels exactly (%s): their levels appear ' ...
        'unequally often. Exact balance of all attributes needs a multiple of %d rows%s.'], S, ...
        strjoin(arrayfun(@num2str, unique(C.L(unbal)), 'UniformOutput', false), ' and '), ...
        strjoin(unique(C.colAttr(unbal), 'stable'), ', '), lcmAll, ...
        ternary(isempty(opts), '', sprintf(' (nearest: %s)', strjoin(arrayfun(@num2str, opts, 'UniformOutput', false), ' or ')))));
end
few = find(S ./ C.L < 2);
if ~isempty(few)
    adv = add(adv, 'rows', 'warning', sprintf('Each level of %s appears fewer than twice in %d rows; its effects rest on single choice tasks.', ...
        strjoin(unique(C.colAttr(few), 'stable'), ', '), S));
end

% ---- blocks ----------------------------------------------------------------------
if B > 1
    T = S / B;
    if T ~= round(T)
        adv = add(adv, 'blocks', 'suggestion', sprintf('%d rows do not divide into %d equal blocks; respondents will answer %d or %d tasks.', ...
            S, B, floor(T), ceil(T)));
    end
    Tf = floor(T);
    tooMany = find(C.L > Tf);
    if ~isempty(tooMany)
        adv = add(adv, 'blocks', 'warning', sprintf(['A block has %d tasks, but %s has more levels than that (%s): no respondent can see all its ' ...
            'levels. Fewer blocks (%s) would allow it.'], Tf, strjoin(unique(C.colAttr(tooMany), 'stable'), ', '), ...
            strjoin(arrayfun(@num2str, unique(C.L(tooMany)), 'UniformOutput', false), ', '), blockOptions(S, max(C.L), B)));
    else
        odd = find(mod(Tf, C.L) ~= 0);
        if ~isempty(odd) && T == Tf
            good = goodBlockCounts(S, C.L);
            adv = add(adv, 'blocks', 'suggestion', sprintf(['Blocks of %d tasks cannot show the levels of %s equally often (%s levels). ' ...
                '%s'], Tf, strjoin(unique(C.colAttr(odd), 'stable'), ', '), ...
                strjoin(arrayfun(@num2str, unique(C.L(odd)), 'UniformOutput', false), ' and '), ...
                ternary(isempty(good), 'No number of blocks of this design allows perfect balance within blocks; a different number of rows would.', ...
                sprintf('With %s block(s) every block could be perfectly balanced.', strjoin(arrayfun(@num2str, good, 'UniformOutput', false), ' or ')))));
        end
    end
    if T > 16
        adv = add(adv, 'blocks', 'suggestion', sprintf('Each respondent answers %g choice tasks; more blocks would shorten the task (%d blocks: %g tasks each).', ...
            T, 2 * B, T / 2));
    elseif T < 4
        adv = add(adv, 'blocks', 'suggestion', sprintf(['Each respondent answers only %g choice tasks. Fewer blocks would give more choices per respondent, ' ...
            'so fewer respondents for the same precision, and more information on preference heterogeneity in panel models.'], T));
    end
end

% ---- priors ----------------------------------------------------------------------
for m = 1:C.M
    Mm = C.model(m);
    if all(Mm.b0 == 0)
        adv = add(adv, 'priors', 'warning', sprintf(['Model %s: all priors are zero. The design is then optimal only if no attribute matters ' ...
            '(all alternatives equally likely), and dominated alternatives cannot be detected. Even rough priors with the right signs help.'], Mm.name));
    end
    ib = find(Mm.bdist == 1);
    wrong = zeros(size(ib));
    for q = 1:numel(ib)
        mu = Mm.ba(ib(q));
        sd = Mm.bb(ib(q));
        if mu ~= 0
            wrong(q) = 0.5 * erfc(abs(mu) / sd / sqrt(2));
        end
    end
    hit = wrong > 0.10;
    if any(hit)
        txt = arrayfun(@(q) sprintf('%s (%.0f%%)', Mm.parNames{ib(q)}, 100 * wrong(q)), find(hit), 'UniformOutput', false);
        adv = add(adv, 'priors', 'info', sprintf(['Model %s: the Bayesian priors put more than 10%% probability on the opposite sign for %s. ' ...
            'The design then hedges between both directions for these coefficients, at some cost if the sign is in fact known.'], Mm.name, strjoin(txt, ', ')));
    end
end

% ---- constraints -----------------------------------------------------------------
if feasShare < 0.25
    adv = add(adv, 'constraints', 'info', sprintf(['The constraints exclude %.0f%% of all possible choice tasks. Exact level balance is then hard to ' ...
        'combine with efficiency; spec.frequency with a minimum and a maximum count gives the search room.'], 100 * (1 - feasShare)));
end
if nFeas < 3 * S
    adv = add(adv, 'constraints', 'warning', sprintf(['Only about %s distinct feasible choice tasks exist for a design of %d rows, so the design uses a ' ...
        'large part of them and has little to choose from.'], big(nFeas), S));
end

if isempty(result)
    return
end

% ---- parameters ------------------------------------------------------------------
for m = 1:numel(result.models)
    r = result.models(m);
    if ~r.identified
        adv = add(adv, 'parameters', 'warning', sprintf('Model %s is not identified by this design at the prior means.', r.name));
        continue
    end
    if ~isempty(r.bayes) && r.bayes.shareNotIdentified > 0
        adv = add(adv, 'parameters', 'warning', sprintf('Model %s is not identified for %.1f%% of the prior draws.', r.name, 100 * r.bayes.shareNotIdentified));
    end
    if isfield(r, 'panelIntegration') && ~isempty(r.panelIntegration)
        pg = r.panelIntegration;
        if strcmp(pg.integration, 'prior') && pg.tasks > 4
            adv = add(adv, 'parameters', 'warning', sprintf(['Model %s: the panel information was simulated with the same draws for all choice ' ...
                'sequences (PanelIntegration = ''prior''). With %d tasks per respondent this overstates the information, so the D-error and ' ...
                'the sample sizes are too low unless the number of draws is very large. The default, ''importance'', does not have this problem.'], r.name, pg.tasks));
        elseif strcmp(pg.integration, 'importance') && (pg.pooled || pg.noise > 0.05)
            adv = add(adv, 'parameters', 'suggestion', sprintf(['Model %s: the simulated panel information is noisy for this design ' ...
                '(simulation noise %.0f%% of the information%s). More EvalPanelDraws would make the reported D-error more exact.'], r.name, ...
                100 * max(pg.noise, 0), ternary(pg.pooled, '; the noise-corrected estimate was not positive definite', '')));
        end
    end
    sest = r.Sestimate;
    sest(~r.inCriterion) = NaN;
    if any(isfinite(sest))
        [worst, iw] = max(sest);
        med = median(sest(isfinite(sest)));
        nb = max(result.blocks);
        if worst > 5 * med && worst > 30
            adv = add(adv, 'parameters', 'info', sprintf(['Model %s: %s is the weakest parameter. A t-ratio of 1.96 needs about %.0f replications of the ' ...
                'design (%.0f respondents), against a median of %.0f for the other parameters. A wider range of its attribute, or a larger prior if it ' ...
                'was set cautiously, would change this.'], r.name, r.parameters{iw}, worst, worst * nb, med));
        else
            adv = add(adv, 'parameters', 'info', sprintf(['Model %s: a t-ratio of 1.96 for every parameter with a non-zero prior needs about %.0f replications ' ...
                'of the design, i.e. %.0f respondents%s (weakest: %s). This holds only if the priors are right.'], r.name, worst, worst * nb, ...
                ternary(nb > 1, sprintf(' answering %g tasks each', result.rows / nb), ''), r.parameters{iw}));
        end
    end
end
P = result.probabilities;
dom = mean(result.maxProbability > 0.9);
if dom > 0.1
    adv = add(adv, 'priors', 'suggestion', sprintf(['In %.0f%% of the tasks one alternative is chosen with probability above 0.9 at the prior means. ' ...
        'Such tasks tell little; spec.maxProb = 0.9 excludes them.'], 100 * dom));
end
[shareMax, jmax] = max(result.meanShare);
if shareMax > 0.7 && size(P, 2) > 2
    adv = add(adv, 'priors', 'warning', sprintf(['%s is predicted to be chosen in %.0f%% of all tasks. If that is not expected in the survey, ' ...
        'check its constant''s prior; if it is, most choices will carry little information on trade-offs.'], result.alternatives{jmax}, 100 * shareMax));
end

% ---- efficiency bound -------------------------------------------------------------
if enumerable && ~isempty(Dfeas) && size(Dfeas, 1) <= 5000 && strcmp(M1.type, 'mnl') && all(M1.eff) && strcmp(C.criterion, 'd') && result.models(1).identified
    bound = dBound(M1, Dfeas, S);
    eff = bound / result.models(1).DerrorFixed;
    if isfinite(eff) && eff > 0 && eff <= 1 + 1e-6
        adv = add(adv, 'efficiency', 'info', sprintf(['Model %s at the prior means: no design of %d feasible tasks can have a D-error below %.5g, even without ' ...
            'level balance or whole numbers of tasks. This design has %.5g, i.e. at least %.1f%% of the highest possible D-efficiency.'], ...
            M1.name, S, bound, result.models(1).DerrorFixed, 100 * min(eff, 1)));
    end
end

% ---- search -----------------------------------------------------------------------
if ~isempty(result.search) && isfield(result.search, 'startObjectives')
    s = result.search;
    best = min(s.startObjectives);
    near = nnz(s.startObjectives <= best * 1.01);
    if isfield(s, 'method')
        adv = add(adv, 'search', 'info', sprintf('Search method: %s (%s).', methodName(s.method), s.methodReason));
    end
    if s.starts >= 3 && near < 2
        adv = add(adv, 'search', 'suggestion', sprintf(['Only the best of %d starts came within 1%% of the best objective (range %.4g to %.4g). ' ...
            'More starts (''Restarts'') may still find a better design.'], s.starts, best, max(s.startObjectives)));
    elseif s.starts >= 3
        adv = add(adv, 'search', 'info', sprintf(['%d of %d starts ended within 1%% of the best objective (range %.4g to %.4g), so more starts are ' ...
            'unlikely to change the result much.'], near, s.starts, best, max(s.startObjectives)));
    end
    if s.timedOut
        adv = add(adv, 'search', 'warning', 'The time limit stopped the search; the design is the best found until then.');
    end
end
if ~isempty(result.blocking) && isfield(result.blocking, 'absentAvoidable') && result.blocking.absentAvoidable > 0
    adv = add(adv, 'blocks', 'warning', sprintf('%d attribute level(s) never appear in some block although the block is large enough.', result.blocking.absentAvoidable));
end
end

% =============================================================================
function adv = add(adv, topic, level, text)
adv(end + 1) = struct('topic', topic, 'level', level, 'text', text);
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end

function s = big(x)
if x < 1e4
    s = sprintf('%.0f', x);
elseif x < 1e6
    s = sprintf('%.0f thousand', x / 1e3);
elseif x < 1e9
    s = sprintf('%.1f million', x / 1e6);
elseif x < 1e12
    s = sprintf('%.1f billion', x / 1e9);
else
    s = sprintf('10^%.0f', log10(x));
end
end

function s = duration(sec)
if sec < 90
    s = sprintf('%.0f s', max(sec, 1));
elseif sec < 5400
    s = sprintf('%.0f min', sec / 60);
else
    s = sprintf('%.1f h', sec / 3600);
end
end

function s = longDuration(years)
if years < 1 / 365
    s = sprintf('about %.0f hours', max(1, years * 365 * 24));
elseif years < 1
    s = sprintf('about %.0f days', years * 365);
elseif years < 1000
    s = sprintf('about %.0f years', years);
else
    s = sprintf('about 10^%.0f years', log10(years));
end
end

function rate = evaluationRate(C)
% Designs evaluated per second for model 1 at its prior means (order of magnitude).
Mm = C.model(1);
if ~strcmp(Mm.type, 'mnl')
    rate = 2e4;
    return
end
n = 2000;
rs = RandStream('mt19937ar', 'Seed', 1);
Dr = zeros(n, C.K);
for k = 1:C.K
    Dr(:, k) = randi(rs, C.L(k), n, 1);
end
d = struct('R', 1, 'B', Mm.b0);
I = dd_rowinfo(Mm, dd_X(Mm, Dr), d);
A = sum(I(:, :, 1:min(n, 4 * Mm.P), :), 3) + I;
t = tic;
for i = 1:5
    f = dd_crit(A, Mm, 'd'); %#ok<NASGU>
end
el = toc(t) / 5;
rate = 10^round(log10(max(n / max(el, 1e-6), 1)));
end

function txt = blockOptions(S, Lmax, B)
opts = [];
for b = 1:B - 1
    if mod(S, b) == 0 && S / b >= Lmax
        opts(end + 1) = b; %#ok<AGROW>
    end
end
if isempty(opts)
    txt = 'none with this number of rows';
else
    txt = sprintf('%s', strjoin(arrayfun(@num2str, opts(max(1, end - 2):end), 'UniformOutput', false), ' or '));
end
end

function good = goodBlockCounts(S, L)
good = [];
for b = 2:S
    if mod(S, b) == 0 && all(mod(S / b, L) == 0)
        good(end + 1) = b; %#ok<AGROW>
    end
end
end

function n = methodName(m)
switch m
    case 'swap', n = 'exchanges within columns';
    case 'rowexchange', n = 'row exchange';
    case 'both', n = 'row exchange and exchanges within columns';
    case 'exhaustive', n = 'all designs evaluated';
    otherwise, n = m;
end
end

function bound = dBound(Mm, Dfeas, S)
% Lower bound on the D-error of any design of S tasks drawn from Dfeas: the
% D-optimal approximate design (weights on the candidate tasks, found by the
% multiplicative algorithm), with the standard correction for not having
% converged: log det M* <= log det M + max_i tr(M^-1 I_i) - P.
P = Mm.P;
I = dd_rowinfo(Mm, dd_X(Mm, Dfeas), struct('R', 1, 'B', Mm.b0));
n = size(I, 3);
w = ones(n, 1) / n;
I2 = reshape(I, P * P, n);
bound = NaN;
for it = 1:2000
    M = reshape(I2 * w, P, P);
    [R, flag] = chol((M + M') / 2);
    if flag ~= 0
        return
    end
    Minv = R \ (R' \ eye(P));
    d = (Minv(:)' * I2)';
    if max(d) <= P * (1 + 1e-5)
        break
    end
    w = w .* d / P;
    w = w / sum(w);
end
ld = 2 * sum(log(diag(R))) + (max(d) - P);
bound = exp(-ld / P) / S;
end
