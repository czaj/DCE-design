function result = dce_design(spec, varargin)
%DCE_DESIGN Generate an efficient design for a discrete choice experiment.
%
%   result = DCE_DESIGN(spec) searches for the design that minimises the
%   D-error (or A-error) of the model(s) in spec, subject to its constraints,
%   assigns blocks, evaluates the design and returns everything in one struct.
%   spec is a struct, a JSON file with the same fields, or an Ngene script
%   (.ngs file, translated by dce_ngene). SPEC.md documents every field; the
%   essentials are:
%
%     spec.alternatives = {'A','B','sq'};
%     spec.rows   = 24;
%     spec.blocks = 3;
%     spec.attributes.price   = [10 20 30 40];
%     spec.attributes.quality = [0 1 2];
%     spec.model.type = 'mnl';                    % 'mnl', 'rp' or 'rppanel'
%     spec.model.utility.A  = 'b_price*price + b_q.dummy*quality';
%     spec.model.utility.B  = 'b_price*price + b_q.dummy*quality';
%     spec.model.utility.sq = 'asc_sq';
%     spec.model.priors.b_price = -0.05;          % fixed prior
%     spec.model.priors.b_q     = {'N(0.3,0.1)', 'N(0.6,0.1)'};   % Bayesian priors
%     spec.model.priors.asc_sq  = 0.2;
%
%   spec.type selects another kind of design: 'full' (all feasible choice
%   tasks, or a random fraction of spec.rows of them; see dce_full),
%   'orthogonal' or 'orthogonal-sequential' (see dce_orthogonal). The default,
%   'efficient', is described here.
%
%   result = DCE_DESIGN(spec, Name, Value, ...) sets options (also accepted in
%   spec.search). The main ones:
%     'Seed'            1     one seed fixes the starts and all draws
%     'Restarts'        12    independent random starts (4 for costly problems)
%     'Kicks'           3     perturb-and-reoptimise rounds per start (1)
%     'Method'          'auto' 'swap' (exchanges within columns), 'rowexchange'
%                             (whole tasks from the feasible candidates), 'both',
%                             or 'exhaustive' (every design; small problems).
%                             'auto' uses 'both' where the level frequencies
%                             have room and 'swap' where they are fixed exactly.
%     'MaxTime'         Inf   seconds; the best design so far is then returned
%     'BayesDraws'      128   prior draws in the search (reported values use
%                             'EvalBayesDraws', 2048)
%     'Parallel'        'auto' use an open parallel pool; true opens one
%     'Start'           []    starting design (level values) or an earlier result
%     'PanelSurrogate'  'mnl' how panel mixed logit models enter the search
%     'Normalize'       []    several models: scale each by its own optimum
%     'Verbose'         1     0 silent, 1 summary, 2 progress
%   See dd_options (private) or README.md for the full list.
%
%   Main fields of result:
%     design          table: Task, Block and one column per attribute of each
%                     alternative, in level values, sorted by block
%     objective       weighted criterion; models(m).Derror etc. per model
%     models(m)       parameters, priors, D-error, A-error, information matrix,
%                     covariance, standard errors and S-estimates for one
%                     replication of the design, Bayesian summary
%     probabilities   predicted choice shares of every task at the prior means
%     checks          constraints, duplicate tasks, level frequencies
%     blocking        balance of the blocks (see dce_block)
%     search          options, method, starts, time, candidates considered
%     advice          remarks on the specification and the design (dce_advise)
%     notes           what the tool assumed or wants you to know
%
%   The same seed and specification always give the same design.
%
%   See also dce_check, dce_evaluate, dce_advise, dce_block, dce_full,
%   dce_orthogonal, dce_report, dce_export, dce_ngene, dce_cli.

t0 = tic;
spec = dd_readspec(spec);
if isfield(spec, 'type') && ~isempty(spec.type)
    switch lower(char(spec.type))
        case 'efficient'
        case 'full'
            args = varargin;
            if isfield(spec, 'rows') && ~isempty(spec.rows)
                args = [{'Rows', spec.rows} args];
            end
            result = dce_full(spec, args{:});
            return
        case {'orthogonal', 'orthogonal-simultaneous'}
            result = dce_orthogonal(spec, 'Type', 'simultaneous', varargin{:});
            return
        case 'orthogonal-sequential'
            result = dce_orthogonal(spec, 'Type', 'sequential', varargin{:});
            return
        otherwise
            dd_err('spec:type', 'spec.type must be ''efficient'', ''full'', ''orthogonal'' or ''orthogonal-sequential''.');
    end
end
C = dd_compile(spec);
o = dd_options(C.spec, varargin{:});
say = @(lvl, varargin) printIf(o.Verbose >= lvl, varargin{:});
dd_precheck(C);

% ---- starting design ---------------------------------------------------------
D0 = [];
start = o.Start;
if isempty(start) && isfield(spec, 'start') && ~isempty(spec.start)
    start = spec.start;
end
if ~isempty(start)
    if isstruct(start) && isfield(start, 'design')
        start = start.design;
    elseif ischar(start) || (isstring(start) && isscalar(start))
        [~, ~, ext] = fileparts(char(start));
        if strcmpi(ext, '.ngd')
            start = dce_readngd(char(start));
        else
            start = readtable(char(start), 'VariableNamingRule', 'preserve');
        end
    end
    D0 = dd_tolevels(C, start, 'The starting design');
    if size(D0, 1) ~= C.S
        dd_err('design:rows', 'The starting design has %d rows; spec.rows is %d.', size(D0, 1), C.S);
    end
    bad = ~dd_rowok(C, D0) | dd_dups(C, D0);
    okFreq = true;
    for k = 1:C.K
        c = accumarray(D0(:, k), 1, [C.L(k) 1]);
        okFreq = okFreq && all(c >= C.lo{k} & c <= C.hi{k});
    end
    if any(bad) || ~okFreq
        C.notes{end + 1} = sprintf(['The starting design breaks the constraints or level frequencies (%d tasks), ' ...
            'so it was not used as a start.'], nnz(bad));
        D0 = [];
    end
end

% ---- search -------------------------------------------------------------------
Cs = dd_searchmodel(C, o);
DR = dd_draws(Cs, o.BayesDraws, o.RandomDraws, o.Seed, o.DrawType, 0);
act = find(C.w > 0)';
% default effort: cheap problems get more starts
costly = any([DR(act).R] > 1) || any(strcmp({Cs.model(act).type}, 'rp')) || max([Cs.model(act).Pt]) >= 25 || C.S * C.K >= 600;
if isempty(o.Restarts)
    if costly, o.Restarts = 4; else, o.Restarts = 12; end
end
if isempty(o.Kicks)
    if costly, o.Kicks = 1; else, o.Kicks = 3; end
end
% search method: exchanges within columns always work; whole tasks can be
% replaced only where the level frequencies are not fixed exactly
freedom = any(cellfun(@(a, b) any(b > a), C.lo, C.hi));
method = o.Method;
if strcmp(method, 'exhaustive')
    methodReason = 'every design is evaluated, as requested';
elseif strcmp(method, 'auto')
    if freedom
        method = 'both';
        methodReason = ['level frequencies have room, so whole tasks are replaced from the feasible candidates (row exchange) ' ...
            'in turn with exchanges within columns'];
    else
        method = 'swap';
        methodReason = 'level frequencies are fixed exactly, so the design can only change by exchanges within columns';
    end
elseif ~freedom && ~strcmp(method, 'swap')
    methodReason = sprintf('Method ''%s'' was requested, but the level frequencies are fixed exactly, so only exchanges within columns can change the design', method);
    method = 'swap';
else
    methodReason = 'as requested';
end
o.method = method;
pool = [];
if ~strcmp(method, 'swap')
    pool = dd_candidates(Cs, DR, o);
    say(2, 'Row exchange: %d candidate tasks%s\n', size(pool.D, 1), ternaryText(pool.all, ' (all feasible tasks)', ' (a random sample of the feasible tasks)'));
end
normalize = o.Normalize;
if isempty(normalize)
    normalize = numel(act) > 1;
end
normalize = normalize && numel(act) > 1;
ref = ones(C.M, 1);
kappa = ones(C.M, 1);
panelIdx = act(strcmp({C.model(act).type}, 'rppanel'));
if normalize
    % each model is scaled by the best value found for it alone
    say(2, 'Reference optimum of each model:\n');
    for m = act
        Cm = Cs;
        Cm.w = zeros(C.M, 1);
        Cm.w(m) = 1;
        if strcmp(method, 'exhaustive')
            [~, ref(m)] = dd_exhaustive(Cm, DR, o, pool);
        else
            cm = dd_search(Cm, DR, o, D0, max(2, ceil(o.Restarts / 3)), pool);
            ref(m) = cm(1).f;
        end
        say(2, '  %-12s %.6g\n', C.model(m).name, ref(m));
    end
    Cs.w = C.w ./ ref;
elseif numel(act) > 1 && ~isempty(panelIdx) && ~strcmp(method, 'exhaustive')
    % raw weights refer to the panel criterion; rescale its surrogate so that
    % the weighted sum in the search has the same balance
    c0 = dd_search(Cs, DR, o, D0, max(1, ceil(o.Restarts / 4)), pool);
    [Dk, yk] = blockDesign(C, c0(1).D, o);
    mk = dd_measures(C, Dk, yk, o);
    for m = panelIdx
        fsur = dd_crit(dd_designinfo(Cs, m, Dk, [], DR(m), o), Cs.model(m), C.criterion);
        kappa(m) = mk.models(m).criterion / fsur;
    end
    Cs.w = C.w .* kappa;
    D0 = c0(1).D;
    say(2, 'Surrogate scale of the panel model(s): %s\n', mat2str(kappa(panelIdx)', 4));
end
exhaustive = [];
if strcmp(method, 'exhaustive')
    tx = tic;
    [Dx, fx, nTried, nValid] = dd_exhaustive(Cs, DR, o, pool);
    cand = struct('D', Dx, 'f', fx);
    slog = struct('passes', 0, 'evaluations', nValid, 'moves', 0, 'timedOut', false, 'rowMoves', 0, 'kicks', 0, ...
        'seconds', toc(tx), 'objective', fx);
    exhaustive = struct('sets', nTried, 'designs', nValid);
else
    [cand, slog] = dd_search(Cs, DR, o, D0, [], pool);
end
say(2, 'Search (%s): %d starts, %d distinct local optima, best search objective %.6g\n', method, numel(slog), numel(cand), cand(1).f);
nDistinct = numel(cand);
panel = any(strcmp({C.model.type}, 'rppanel'));
nc = min(o.Candidates, numel(cand));
cand = cand(1:nc);
candMethod = repmat({method}, 1, nc);
if panel && strcmp(method, 'both') && strcmp(o.Method, 'auto')
    % A panel model is searched through a surrogate, and in some problems the
    % surrogate gains more from the freedom of row exchange than the panel
    % criterion does. So the search is also run with exchanges within columns
    % only, and the panel criterion chooses among the designs of both.
    o2 = o;
    o2.method = 'swap';
    [cand2, slog2] = dd_search(Cs, DR, o2, D0, [], []);
    nDistinct = nDistinct + numel(cand2);
    n2 = min(o.Candidates, numel(cand2));
    cand = [cand cand2(1:n2)];
    candMethod = [candMethod repmat({'swap'}, 1, n2)];
    slog = [slog slog2];
    nc = numel(cand);
    methodReason = [methodReason '; because of the panel model the search was also run with exchanges within columns only, ' ...
        'and the designs of both were compared on the panel criterion'];
end

% ---- block and evaluate the best candidates -------------------------------------
% The blocking changes the information only in panel models (a respondent
% answers one block). With panel models every candidate is blocked before it
% is evaluated; otherwise only the chosen design is blocked.
%
% A panel criterion is simulated. Choosing the best of several designs on the
% sample that also gives the reported values would flatter those values, so
% with panel models the choices (among the candidates, and whether to keep
% the result of the second stage) are made on a separate validation sample.
res = cell(nc, 1);
score = zeros(nc, 1);
crit = zeros(nc, C.M);
oV = o;
if panel
    oV.Seed = o.Seed + 7331;
end
for i = 1:nc
    D = cand(i).D;
    y = [];
    bout = [];
    if panel
        [D, y, bout] = blockDesign(C, D, o);
    end
    meas = dd_measures(C, D, y, oV);
    res{i} = struct('D', D, 'y', y, 'bout', bout, 'meas', meas);
    crit(i, :) = [meas.models.criterion];
end
scale = ones(1, C.M);
if normalize
    scale = min(crit, [], 1);
end
for i = 1:nc
    score(i) = sum(C.w(act)' .* crit(i, act) ./ scale(act));
end
[~, best] = min(score);
D = res{best}.D;
y = res{best}.y;
bout = res{best}.bout;
meas = res{best}.meas;
measV = meas;
if panel
    meas = dd_measures(C, D, y, o);       % the reported values: not the sample the design was chosen on
end
polish = [];
if ~panel
    [D, y, bout] = blockDesign(C, D, o);
    meas = dd_measures(C, D, y, o);
elseif o.PanelPolish > 0
    % second stage: the objective with the panel model itself, by exchanges
    % within blocks (the first stage could only use a surrogate)
    [Dp, polish] = dd_panelpolish(C, D, y, o, normalize);
    say(2, 'Panel stage: %d moves in %d pass(es), %d evaluations\n', polish.moves, polish.passes, polish.evaluations);
    polish.accepted = false;
    polish.evalBefore = [meas.models.criterion];
    polish.evalAfter = polish.evalBefore;
    if polish.moves > 0
        % keep the result only if it is also better on the validation sample
        measP = dd_measures(C, Dp, y, oV);
        cb = [measV.models.criterion];
        ca = [measP.models.criterion];
        sc = ones(1, C.M);
        if normalize
            sc = cb;
        end
        if sum(C.w(act)' .* ca(act) ./ sc(act)) < sum(C.w(act)' .* cb(act) ./ sc(act))
            D = Dp;
            meas = dd_measures(C, Dp, y, o);
            polish.accepted = true;
            polish.evalAfter = [meas.models.criterion];
        end
    end
end
result = dd_result(C, D, y, o, meas);
result.blocking = bout;

% ---- search record ----------------------------------------------------------------
sr = struct();
sr.options = rmfield(o, {'deadline', 'Start'});
sr.seconds = toc(t0);
sr.starts = numel(slog);
sr.evaluations = sum([slog.evaluations]);
sr.startObjectives = [slog.objective];
sr.timedOut = any([slog.timedOut]);
sr.distinctOptima = nDistinct;
sr.method = method;
sr.methodReason = methodReason;
sr.moves = sum([slog.moves]);
sr.rowMoves = sum([slog.rowMoves]);
sr.candidatePool = 0;
if ~isempty(pool)
    sr.candidatePool = size(pool.D, 1);
end
sr.exhaustive = exhaustive;
if ~isempty(exhaustive)
    result.notes{end + 1} = sprintf(['All %d designs that satisfy the constraints and level frequencies were evaluated (%d sets of feasible tasks); ' ...
        'this design is the best of them%s.'], exhaustive.designs, exhaustive.sets, ...
        ternaryText(any([DR(act).R] > 1), ' on the prior draws of the search', ''));
end
sr.normalized = normalize;
sr.reference = ref';
sr.surrogateScale = kappa';
sr.candidates = table((1:nc)', candMethod(:), [cand(1:nc).f]', score, crit, (1:nc)' == best, ...
    'VariableNames', {'Candidate', 'Method', 'SearchObjective', 'Score', 'Criterion', 'Chosen'});
sr.chosenMethod = candMethod{best};
sr.panelSurrogate = '';
sr.panelStage = polish;
if panel
    sr.panelSurrogate = o.PanelSurrogate;
    unitTxt = 'one block';
    if strcmp(C.panelUnit, 'design') || C.B < 2
        unitTxt = 'every row';
    end
    if isempty(polish)
        result.notes{end + 1} = sprintf(['Panel mixed logit (a respondent answers %s): the search used the %s surrogate; the reported ' ...
            'values and the choice among the %d best designs use the panel model itself.'], unitTxt, upper(o.PanelSurrogate), nc);
    else
        base = sprintf(['Panel mixed logit (a respondent answers %s): stage 1 searched with the %s surrogate and the best of ' ...
            '%d designs was chosen on the panel criterion; stage 2 optimised the objective with the panel model itself ' ...
            'by exchanges within blocks'], unitTxt, upper(o.PanelSurrogate), nc);
        if polish.accepted
            chg = arrayfun(@(m) sprintf('%s %.4g -> %.4g', C.model(m).name, polish.evalBefore(m), polish.evalAfter(m)), act, 'UniformOutput', false);
            result.notes{end + 1} = sprintf('%s: %d exchange(s), %s.', base, polish.moves, strjoin(chg, ', '));
        elseif polish.moves > 0
            result.notes{end + 1} = sprintf(['%s, but its %d exchange(s) did not improve the objective on the validation sample, ' ...
                'so the stage-1 design was kept. More PanelRespondents or PanelDraws make stage 2 more reliable.'], base, polish.moves);
        else
            result.notes{end + 1} = sprintf('%s and found nothing to improve.', base);
        end
        if polish.timedOut
            sr.timedOut = true;
        end
    end
end
if sr.timedOut
    result.notes{end + 1} = sprintf('MaxTime (%g s) was reached; the design is the best found so far.', o.MaxTime);
end
if sr.normalized
    result.notes{end + 1} = 'Several models: each criterion is divided by the best value found for that model alone before weighting.';
end
result.search = sr;
if ~result.checks.constraintsOK || ~isempty(result.checks.duplicateTasks) || ~result.checks.frequenciesOK
    dd_err('internal:check', 'The final design fails its own checks. Please report this specification.');
end
result.advice = dd_advice(C, result);
if o.Verbose >= 1
    dce_report(result);
end
end

% =============================================================================
function [D, y, bout] = blockDesign(C, D, o)
% Assign blocks with dce_block and sort the rows by block.
S = size(D, 1);
y = ones(S, 1);
bout = [];
if C.B < 2 || ~o.Block
    return
end
V = zeros(S, C.K);
for k = 1:C.K
    V(:, k) = C.LV{k}(D(:, k));
end
Mm = C.model(1);
X = dd_X(Mm, D);
args = {'Names', C.colName, 'Verbose', false, 'Seed', o.Seed};
info = struct('X', permute(X, [3 1 2]), 'beta', Mm.b0);
% lowest and highest level of numeric attributes that enter linearly
lin = false(1, C.K);
for t = 1:Mm.T
    for f = 1:Mm.tNF(t)
        k = Mm.tCol(t, f);
        if isequal(Mm.tTab{t, f}, C.LV{k}) && C.L(k) >= 3
            lin(k) = true;
        end
    end
end
if any(lin)
    args = [args {'Extremes', find(lin)}];
end
[y, b] = dce_block(V, C.B, args{:}, 'Info', info, 'Restarts', o.BlockRestarts, 'Iterations', max(1000, 20 * S));
[y, ord] = sort(y, 'ascend');
D = D(ord, :);
bout = struct('Psi', b.Psi, 'PsiBound', b.PsiBound, 'absentLevels', height(b.absent), ...
    'absentAvoidable', nnz(b.absent.Avoidable), 'blockSizes', b.blockSizes, 'extremeColumns', {C.colName(lin)}, ...
    'extremesOK', isempty(b.extremes) || all(b.extremes.HasLowest & b.extremes.HasHighest), ...
    'informationShare', [], 'leaveOneBlockOut', []);
if ~isempty(b.info)
    bout.informationShare = b.info.share;
    bout.leaveOneBlockOut = b.info.LOBOratio;
end
end

function printIf(cond, varargin)
if cond
    fprintf(varargin{:});
end
end

function t = ternaryText(c, a, b)
if c, t = a; else, t = b; end
end
