function C = dd_compile(spec, evaluating)
%DD_COMPILE Validate a design specification and compile it for the engine. Internal.
%   See SPEC.md for the specification format. Every problem raises an error
%   with identifier dce:spec:* and a message that names the offending field.

if nargin < 2
    evaluating = false;
end
if ~isstruct(spec) || ~isscalar(spec)
    dd_err('spec:type', 'The specification must be a scalar struct (or a JSON object).');
end
known = {'name', 'type', 'alternatives', 'rows', 'blocks', 'attributes', 'altLevels', 'scenario', ...
    'fixed', 'base', 'labels', 'model', 'models', 'weights', 'criterion', 'constants', 'exclude', ...
    'excludeMode', 'reject', 'require', 'dominance', 'maxProb', 'balance', 'frequency', ...
    'unique', 'mirror', 'panelUnit', 'start', 'search', 'notes', 'source'};
checkFields(spec, known, 'spec');

C = struct();
C.version = dd_version();
C.name = char(getf(spec, 'name', 'design'));
C.notes = {};
if isfield(spec, 'notes') && ~isempty(spec.notes)
    C.notes = reshape(cellstr(spec.notes), 1, []);
end

% ---- alternatives ----------------------------------------------------------
if ~isfield(spec, 'alternatives')
    dd_err('spec:alternatives', 'spec.alternatives is required, e.g. {''A'',''B'',''sq''}.');
end
alts = toCellstr(spec.alternatives, 'spec.alternatives');
J = numel(alts);
if J < 2
    dd_err('spec:alternatives', 'At least two alternatives are required.');
end
for j = 1:J
    if ~isvarname(alts{j})
        dd_err('spec:alternatives', 'Alternative name ''%s'' is not a valid identifier (letters, digits, underscore).', alts{j});
    end
end
if numel(unique(alts)) < J
    dd_err('spec:alternatives', 'Alternative names must be unique.');
end
C.J = J;
C.altNames = alts;

% ---- rows and blocks -------------------------------------------------------
S = reqInt(spec, 'rows', 2);
B = 1;
if isfield(spec, 'blocks') && ~isempty(spec.blocks)
    B = reqInt(spec, 'blocks', 1);
end
if B > S
    dd_err('spec:blocks', 'spec.blocks (%d) exceeds spec.rows (%d).', B, S);
end
C.S = S;
C.B = B;

% ---- attributes ------------------------------------------------------------
if ~isfield(spec, 'attributes') || ~isstruct(spec.attributes) || ~isscalar(spec.attributes)
    dd_err('spec:attributes', 'spec.attributes must be a struct: attribute name -> vector of levels.');
end
attrNames = fieldnames(spec.attributes)';
G = struct();
for i = 1:numel(attrNames)
    a = attrNames{i};
    if any(strcmp(a, alts))
        dd_err('spec:attributes', 'Attribute ''%s'' has the same name as an alternative.', a);
    end
    G.(a) = levelVector(spec.attributes.(a), ['spec.attributes.' a], 1);
end

% alternative-specific level sets
AL = struct();
if isfield(spec, 'altLevels') && ~isempty(spec.altLevels)
    if ~isstruct(spec.altLevels)
        dd_err('spec:altLevels', 'spec.altLevels must be a struct: alternative -> attribute -> levels.');
    end
    fa = fieldnames(spec.altLevels)';
    for i = 1:numel(fa)
        mustBeAlt(fa{i}, alts, 'spec.altLevels');
        fb = fieldnames(spec.altLevels.(fa{i}))';
        for q = 1:numel(fb)
            mustBeAttr(fb{q}, attrNames, ['spec.altLevels.' fa{i}]);
            AL.(fa{i}).(fb{q}) = levelVector(spec.altLevels.(fa{i}).(fb{q}), ...
                sprintf('spec.altLevels.%s.%s', fa{i}, fb{q}), 2);
        end
    end
end

% scenario variables
scen = {};
if isfield(spec, 'scenario') && ~isempty(spec.scenario)
    scen = toCellstr(spec.scenario, 'spec.scenario');
    for i = 1:numel(scen)
        mustBeAttr(scen{i}, attrNames, 'spec.scenario');
    end
end

% fixed levels (status quo / opt-out)
FX = struct();
if isfield(spec, 'fixed') && ~isempty(spec.fixed)
    if ~isstruct(spec.fixed)
        dd_err('spec:fixed', 'spec.fixed must be a struct: alternative -> attribute -> value.');
    end
    fa = fieldnames(spec.fixed)';
    for i = 1:numel(fa)
        mustBeAlt(fa{i}, alts, 'spec.fixed');
        fb = fieldnames(spec.fixed.(fa{i}))';
        for q = 1:numel(fb)
            mustBeAttr(fb{q}, attrNames, ['spec.fixed.' fa{i}]);
            v = spec.fixed.(fa{i}).(fb{q});
            if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
                dd_err('spec:fixed', 'spec.fixed.%s.%s must be one finite number.', fa{i}, fb{q});
            end
            if any(strcmp(fb{q}, scen))
                dd_err('spec:fixed', 'Attribute ''%s'' is a scenario variable and cannot be fixed for alternative ''%s''.', fb{q}, fa{i});
            end
            FX.(fa{i}).(fb{q}) = double(v);
        end
    end
end
% an attribute with a single level is a constant wherever it is used
for i = 1:numel(attrNames)
    a = attrNames{i};
    if isscalar(G.(a))
        if any(strcmp(a, scen))
            dd_err('spec:scenario', 'Scenario variable ''%s'' needs at least two levels.', a);
        end
        for j = 1:J
            if ~(isfield(FX, alts{j}) && isfield(FX.(alts{j}), a))
                FX.(alts{j}).(a) = G.(a);
            end
        end
    end
end

% base levels of dummy / effects coded attributes
BASE = struct();
for i = 1:numel(attrNames)
    BASE.(attrNames{i}) = G.(attrNames{i})(1);
end
if isfield(spec, 'base') && ~isempty(spec.base)
    fb = fieldnames(spec.base)';
    for q = 1:numel(fb)
        mustBeAttr(fb{q}, attrNames, 'spec.base');
        v = double(spec.base.(fb{q}));
        if ~isscalar(v) || ~any(G.(fb{q}) == v)
            dd_err('spec:base', 'spec.base.%s must be one of the levels of %s (%s).', fb{q}, fb{q}, num2str(G.(fb{q})));
        end
        BASE.(fb{q}) = v;
    end
end

% ---- models: parse utilities -------------------------------------------------
if isfield(spec, 'model') && isfield(spec, 'models')
    dd_err('spec:models', 'Give either spec.model or spec.models, not both.');
elseif isfield(spec, 'model')
    ms = {spec.model};
elseif isfield(spec, 'models')
    ms = spec.models;
    if isstruct(ms)
        ms = num2cell(ms(:))';
    end
    if ~iscell(ms) || isempty(ms)
        dd_err('spec:models', 'spec.models must be a struct array or a cell array of model structs.');
    end
    ms = ms(:)';
else
    dd_err('spec:models', 'spec.model is required: a struct with fields type, utility, priors.');
end
M = numel(ms);
U = cell(M, J);
mtype = cell(1, M);
for m = 1:M
    mdl = ms{m};
    where = sprintf('model %d', m);
    if ~isstruct(mdl) || ~isscalar(mdl)
        dd_err('spec:models', '%s must be a struct.', where);
    end
    checkFields(mdl, {'name', 'type', 'utility', 'priors', 'random', 'stat'}, where);
    ty = lower(char(getf(mdl, 'type', 'mnl')));
    switch ty
        case 'mnl'
        case {'rp', 'mxl-cross', 'cross'}
            ty = 'rp';
        case {'rppanel', 'mxl', 'panel', 'mxl-panel'}
            ty = 'rppanel';
        otherwise
            dd_err('spec:models', '%s: unknown type ''%s''. Use ''mnl'', ''rp'' (cross-sectional mixed logit) or ''rppanel'' (panel mixed logit).', where, ty);
    end
    mtype{m} = ty;
    if ~isfield(mdl, 'utility') || ~isstruct(mdl.utility)
        dd_err('spec:utility', '%s: utility must be a struct: alternative -> utility string.', where);
    end
    fu = fieldnames(mdl.utility)';
    for q = 1:numel(fu)
        mustBeAlt(fu{q}, alts, [where '.utility']);
    end
    for j = 1:J
        if isfield(mdl.utility, alts{j})
            U{m, j} = parseUtility(mdl.utility.(alts{j}), sprintf('%s, utility of %s', where, alts{j}));
        else
            U{m, j} = parseUtility('', '');
        end
        for t = 1:numel(U{m, j})
            for f = 1:numel(U{m, j}(t).attrs)
                mustBeAttr(U{m, j}(t).attrs{f}, attrNames, sprintf('%s, utility of %s', where, alts{j}));
            end
        end
    end
end

% ---- design columns ----------------------------------------------------------
nA = numel(attrNames);
used = false(J, nA);
scenUsed = false(1, nA);
referenced = false(1, nA);
for m = 1:M
    for j = 1:J
        for t = 1:numel(U{m, j})
            for f = 1:numel(U{m, j}(t).attrs)
                ia = find(strcmp(U{m, j}(t).attrs{f}, attrNames));
                referenced(ia) = true;
                if any(strcmp(attrNames{ia}, scen))
                    scenUsed(ia) = true;
                elseif ~(isfield(FX, alts{j}) && isfield(FX.(alts{j}), attrNames{ia}))
                    used(j, ia) = true;
                end
            end
        end
    end
end
unusedAttr = attrNames(~referenced);
if ~isempty(unusedAttr)
    dd_err('spec:attributes', 'Attribute(s) not used in any utility function: %s. Remove them or add them to a utility.', strjoin(unusedAttr, ', '));
end
colName = {};
colAlt = [];
colAttr = {};
levels = {};
colIndex = zeros(J, nA);
scenCol = zeros(1, nA);
for ia = find(scenUsed)
    colName{end + 1} = attrNames{ia}; %#ok<AGROW>
    colAlt(end + 1) = 0; %#ok<AGROW>
    colAttr{end + 1} = attrNames{ia}; %#ok<AGROW>
    levels{end + 1} = G.(attrNames{ia}); %#ok<AGROW>
    scenCol(ia) = numel(colName);
end
for j = 1:J
    for ia = 1:nA
        if used(j, ia)
            a = attrNames{ia};
            colName{end + 1} = [alts{j} '.' a]; %#ok<AGROW>
            colAlt(end + 1) = j; %#ok<AGROW>
            colAttr{end + 1} = a; %#ok<AGROW>
            if isfield(AL, alts{j}) && isfield(AL.(alts{j}), a)
                levels{end + 1} = AL.(alts{j}).(a); %#ok<AGROW>
            else
                levels{end + 1} = G.(a); %#ok<AGROW>
            end
            colIndex(j, ia) = numel(colName);
        end
    end
end
K = numel(colName);
if K == 0
    dd_err('spec:attributes', 'No attribute varies in the design (all are fixed).');
end
C.K = K;
C.colName = colName;
C.colAlt = colAlt;
C.colAttr = colAttr;
C.levels = levels;
C.L = cellfun(@numel, levels);
C.LV = cellfun(@(v) v(:), levels, 'UniformOutput', false);
C.attrNames = attrNames;
C.colIndex = colIndex;
C.scenCol = scenCol;
C.fixed = FX;
C.base = BASE;

ctx = struct('alts', {alts}, 'attrNames', {attrNames}, 'G', G, 'BASE', BASE, 'FX', FX, ...
    'colIndex', colIndex, 'scenCol', scenCol, 'levels', {levels}, 'scen', {scen}, 'K', K, 'J', J);

% ---- compile models ------------------------------------------------------------
for m = 1:M
    Mm = compileModel(ms{m}, m, mtype{m}, U(m, :), ctx);
    if m == 1
        C.model = Mm;
    else
        C.model(m) = Mm;
    end
end
C.M = M;

% ---- criterion, weights, excluded parameters ------------------------------------
crit = lower(char(getf(spec, 'criterion', 'd')));
if ~any(strcmp(crit, {'d', 'a'}))
    dd_err('spec:criterion', 'spec.criterion must be ''d'' (D-error) or ''a'' (A-error).');
end
C.criterion = crit;
w = ones(M, 1);
if isfield(spec, 'weights') && ~isempty(spec.weights)
    w = double(spec.weights(:));
    if numel(w) ~= M || any(~isfinite(w)) || any(w < 0) || ~any(w > 0)
        dd_err('spec:weights', 'spec.weights must have one non-negative weight per model (%d).', M);
    end
end
C.w = w;
exMode = lower(char(getf(spec, 'excludeMode', 'nuisance')));
if ~any(strcmp(exMode, {'nuisance', 'known'}))
    dd_err('spec:excludeMode', 'spec.excludeMode must be ''nuisance'' or ''known''.');
end
ex = {};
if isfield(spec, 'exclude') && ~isempty(spec.exclude)
    ex = toCellstr(spec.exclude, 'spec.exclude');
end
% constants (parameters that multiply no attribute) are part of the criterion
% unless spec.constants is false
C.constantNames = {};
for m = 1:M
    for j = 1:J
        for t = 1:numel(U{m, j})
            if isempty(U{m, j}(t).attrs)
                C.constantNames = union(C.constantNames, {U{m, j}(t).par}, 'stable');
            end
        end
    end
end
cin = getf(spec, 'constants', true);
if ~(islogical(cin) || isnumeric(cin)) || ~isscalar(cin)
    dd_err('spec:constants', 'spec.constants must be true (constants are part of the criterion; default) or false.');
end
C.constantsIn = logical(cin);
if ~C.constantsIn && ~isempty(C.constantNames)
    ex = union(ex, C.constantNames, 'stable');
    C.notes{end + 1} = sprintf('spec.constants = false: the constants (%s) are estimated but left out of the criterion.', strjoin(C.constantNames, ', '));
end
for m = 1:M
    Mm = C.model(m);
    eff = true(Mm.Pt, 1);
    for q = 1:numel(ex)
        hit = strcmp(ex{q}, Mm.thBase);
        if ~any(hit) && m == 1 && ~any(cellfun(@(mm) any(strcmp(ex{q}, mm)), {C.model.thBase}))
            dd_err('spec:exclude', 'spec.exclude names ''%s'', which is not a parameter of any model. Parameters: %s.', ...
                ex{q}, strjoin(unique(Mm.thBase, 'stable'), ', '));
        end
        eff(hit) = false;
    end
    if ~any(eff)
        dd_err('spec:exclude', 'Model %d: all parameters are excluded from the criterion.', m);
    end
    C.model(m).eff = eff;
    C.model(m).exMode = exMode;
    C.model(m).Peff = nnz(eff);
end

% ---- identification count -------------------------------------------------------
for m = 1:M
    if S * (J - 1) < C.model(m).Pt
        dd_err('spec:rows', 'Model %d has %d parameters but %d rows with %d alternatives give at most %d independent observations.', ...
            m, C.model(m).Pt, S, J, S * (J - 1));
    end
end

% ---- level frequency bounds ------------------------------------------------------
bal = lower(char(getf(spec, 'balance', 'strict')));
if ~any(strcmp(bal, {'strict', 'none'}))
    dd_err('spec:balance', 'spec.balance must be ''strict'' (default) or ''none''. Use spec.frequency for bounds per attribute.');
end
C.balance = bal;
lo = cell(1, K);
hi = cell(1, K);
for k = 1:K
    if strcmp(bal, 'strict')
        lo{k} = floor(S / C.L(k)) * ones(C.L(k), 1);
        hi{k} = ceil(S / C.L(k)) * ones(C.L(k), 1);
    else
        lo{k} = zeros(C.L(k), 1);
        hi{k} = S * ones(C.L(k), 1);
    end
end
if isfield(spec, 'frequency') && ~isempty(spec.frequency)
    if ~isstruct(spec.frequency)
        dd_err('spec:frequency', 'spec.frequency must be a struct: attribute -> struct with exact, or min and max.');
    end
    fb = fieldnames(spec.frequency)';
    for q = 1:numel(fb)
        mustBeAttr(fb{q}, attrNames, 'spec.frequency');
        fq = spec.frequency.(fb{q});
        if ~isstruct(fq)
            dd_err('spec:frequency', 'spec.frequency.%s must be a struct with field exact, or fields min and max.', fb{q});
        end
        checkFields(fq, {'exact', 'min', 'max'}, ['spec.frequency.' fb{q}]);
        for k = find(strcmp(fb{q}, colAttr))
            Lk = C.L(k);
            if isfield(fq, 'exact')
                e = expandFreq(fq.exact, Lk, fb{q}, 'exact');
                lo{k} = e;
                hi{k} = e;
            else
                if isfield(fq, 'min'), lo{k} = expandFreq(fq.min, Lk, fb{q}, 'min'); else, lo{k} = zeros(Lk, 1); end
                if isfield(fq, 'max'), hi{k} = expandFreq(fq.max, Lk, fb{q}, 'max'); else, hi{k} = S * ones(Lk, 1); end
            end
        end
    end
end
for k = 1:K
    if ~evaluating && (any(lo{k} > hi{k}) || sum(lo{k}) > S || sum(hi{k}) < S)
        dd_err('spec:frequency', 'Level frequencies of %s cannot be met with %d rows (min %s, max %s).', ...
            colName{k}, S, mat2str(lo{k}'), mat2str(hi{k}'));
    end
end
C.lo = lo;
C.hi = hi;

% ---- row constraints ---------------------------------------------------------------
C.rejectExpr = {};
C.rejectFun = {};
C.requireExpr = {};
C.requireFun = {};
if isfield(spec, 'reject') && ~isempty(spec.reject)
    C.rejectExpr = toCellstr(spec.reject, 'spec.reject');
    for q = 1:numel(C.rejectExpr)
        C.rejectFun{q} = compileExpr(C.rejectExpr{q}, ctx, sprintf('spec.reject{%d}', q));
    end
end
if isfield(spec, 'require') && ~isempty(spec.require)
    C.requireExpr = toCellstr(spec.require, 'spec.require');
    for q = 1:numel(C.requireExpr)
        C.requireFun{q} = compileExpr(C.requireExpr{q}, ctx, sprintf('spec.require{%d}', q));
    end
end

% ---- dominance ------------------------------------------------------------------------
[C.dom, note] = compileDominance(getf(spec, 'dominance', 'auto'), U(1, :), C, ctx);
C.notes = [C.notes note];

mp = getf(spec, 'maxProb', []);
if ~isempty(mp) && (~isnumeric(mp) || ~isscalar(mp) || mp <= 1 / J || mp > 1)
    dd_err('spec:maxProb', 'spec.maxProb must be a number in (1/J, 1].');
end
C.maxProb = double(mp);
pu = lower(char(getf(spec, 'panelUnit', 'block')));
if ~any(strcmp(pu, {'block', 'design'}))
    dd_err('spec:panelUnit', 'spec.panelUnit must be ''block'' (a respondent answers one block; default) or ''design'' (a respondent answers every row, as Ngene assumes).');
end
C.panelUnit = pu;
C.unique = logical(getf(spec, 'unique', true));
C.mirror = logical(getf(spec, 'mirror', true));
C.mirrorPairs = mirrorPairs(C);

rsH = RandStream('mt19937ar', 'Seed', 12345);
C.hashW = randi(rsH, 2^20, K, 1);

% ---- labels ------------------------------------------------------------------------------
C.labels = struct();
if isfield(spec, 'labels') && ~isempty(spec.labels)
    fb = fieldnames(spec.labels)';
    for q = 1:numel(fb)
        mustBeAttr(fb{q}, attrNames, 'spec.labels');
        lab = toCellstr(spec.labels.(fb{q}), ['spec.labels.' fb{q}]);
        if numel(lab) ~= numel(G.(fb{q}))
            dd_err('spec:labels', 'spec.labels.%s needs one label per level (%d).', fb{q}, numel(G.(fb{q})));
        end
        C.labels.(fb{q}) = lab;
    end
end
C.globalLevels = G;

% ---- notes ---------------------------------------------------------------------------------
for m = 1:M
    z = C.model(m).zeroPrior;
    if ~isempty(z)
        C.notes{end + 1} = sprintf('Model %d: zero priors for %s (the design is then optimised as if these attributes did not matter).', ...
            m, strjoin(z, ', '));
    end
end
C.spec = spec;
C.hash = dd_spechash(spec);
end

% =============================================================================
function Mm = compileModel(mdl, m, ty, U, ctx)
J = ctx.J;
where = sprintf('model %d', m);
reg = struct('name', {}, 'coding', {}, 'attr', {}, 'idx', {}, 'nb', {}, 'inline', {});
parNames = {};
parBase = {};
tAlt = [];
tPar = [];
tConst = [];
tNF = [];
tCol = zeros(0, 3);
tTab = cell(0, 3);
for j = 1:J
    for t = 1:numel(U{j})
        term = U{j}(t);
        ri = find(strcmp(term.par, {reg.name}));
        coded = ~strcmp(term.coding, 'linear');
        if coded && isempty(term.attrs)
            dd_err('spec:utility', '%s: parameter %s is %s coded but multiplies no attribute.', where, term.par, term.coding);
        end
        if isempty(ri)
            r = struct('name', term.par, 'coding', term.coding, 'attr', '', 'idx', [], 'nb', [], 'inline', term.prior);
            if coded
                a = term.attrs{1};
                gl = ctx.G.(a);
                if numel(gl) < 2
                    dd_err('spec:utility', '%s: attribute %s has one level and cannot be %s coded.', where, a, term.coding);
                end
                r.attr = a;
                r.nb = gl(gl ~= ctx.BASE.(a));
                nm = arrayfun(@(v) sprintf('%s[%s]', term.par, num2str(v)), r.nb, 'UniformOutput', false);
            else
                nm = {term.par};
            end
            r.idx = numel(parNames) + (1:numel(nm));
            parNames = [parNames nm]; %#ok<AGROW>
            parBase = [parBase repmat({term.par}, 1, numel(nm))]; %#ok<AGROW>
            reg(end + 1) = r; %#ok<AGROW>
            ri = numel(reg);
        else
            if ~strcmp(reg(ri).coding, term.coding)
                dd_err('spec:utility', '%s: parameter %s is used with two codings (%s and %s). Write the coding suffix wherever the parameter appears.', ...
                    where, term.par, reg(ri).coding, term.coding);
            end
            if coded && ~strcmp(reg(ri).attr, term.attrs{1})
                dd_err('spec:utility', '%s: %s coded parameter %s is used with two attributes (%s and %s).', ...
                    where, term.coding, term.par, reg(ri).attr, term.attrs{1});
            end
            if ~isempty(term.prior)
                if isempty(reg(ri).inline)
                    reg(ri).inline = term.prior;
                elseif ~strcmp(strrep(reg(ri).inline, ' ', ''), strrep(term.prior, ' ', ''))
                    dd_err('spec:priors', '%s: parameter %s has two different inline priors.', where, term.par);
                end
            end
        end
        for e = 1:numel(reg(ri).idx)
            const = 1;
            nf = 0;
            cols = [0 0 0];
            tabs = cell(1, 3);
            for f = 1:numel(term.attrs)
                a = term.attrs{f};
                ia = find(strcmp(a, ctx.attrNames));
                isCoded = coded && f == 1;
                if ctx.scenCol(ia) > 0 && any(strcmp(a, ctx.scen))
                    k = ctx.scenCol(ia);
                    src = ctx.levels{k}(:);
                elseif isfield(ctx.FX, ctx.alts{j}) && isfield(ctx.FX.(ctx.alts{j}), a)
                    k = 0;
                    src = ctx.FX.(ctx.alts{j}).(a);
                    if isCoded && ~any(ctx.G.(a) == src)
                        dd_err('spec:fixed', ['%s: attribute %s is %s coded, but its fixed value for alternative %s (%g) is not one of its levels (%s). ' ...
                            'Add the value to spec.attributes.%s and restrict the designed alternatives with spec.altLevels.'], ...
                            where, a, term.coding, ctx.alts{j}, src, num2str(ctx.G.(a)), a);
                    end
                else
                    k = ctx.colIndex(j, ia);
                    src = ctx.levels{k}(:);
                end
                if isCoded
                    val = double(src == reg(ri).nb(e));
                    if strcmp(term.coding, 'effects')
                        val = val - double(src == ctx.BASE.(a));
                    end
                else
                    val = src;
                end
                if k == 0
                    const = const * val;
                else
                    nf = nf + 1;
                    if nf > 3
                        dd_err('spec:utility', '%s: a term of parameter %s has more than three varying factors.', where, term.par);
                    end
                    cols(nf) = k;
                    tabs{nf} = val;
                end
            end
            if const ~= 0
                tAlt(end + 1) = j; %#ok<AGROW>
                tPar(end + 1) = reg(ri).idx(e); %#ok<AGROW>
                tConst(end + 1) = const; %#ok<AGROW>
                tNF(end + 1) = nf; %#ok<AGROW>
                tCol(end + 1, :) = cols; %#ok<AGROW>
                tTab(end + 1, :) = tabs; %#ok<AGROW>
            end
        end
    end
end
P = numel(parNames);
if P == 0
    dd_err('spec:utility', '%s has no parameters.', where);
end
missing = setdiff(1:P, tPar);
if ~isempty(missing)
    dd_err('spec:utility', '%s: parameter(s) %s never enter a utility with a non-zero regressor (check fixed levels and base levels).', ...
        where, strjoin(parNames(missing), ', '));
end

% priors of the coefficients (means, for random coefficients)
pri = struct();
if isfield(mdl, 'priors') && ~isempty(mdl.priors)
    if ~isstruct(mdl.priors)
        dd_err('spec:priors', '%s: priors must be a struct: parameter -> value.', where);
    end
    pri = mdl.priors;
    fb = fieldnames(pri)';
    for q = 1:numel(fb)
        if ~any(strcmp(fb{q}, {reg.name}))
            dd_err('spec:priors', '%s: priors.%s does not match a parameter. Parameters: %s.', where, fb{q}, strjoin({reg.name}, ', '));
        end
    end
end
b0 = zeros(P, 1);
bdist = zeros(P, 1);
ba = zeros(P, 1);
bb = zeros(P, 1);
zeroPrior = {};
for ri = 1:numel(reg)
    n = numel(reg(ri).idx);
    hasInline = ~isempty(reg(ri).inline);
    hasField = isfield(pri, reg(ri).name);
    if hasInline && hasField
        dd_err('spec:priors', '%s: parameter %s has a prior both inline and in priors.', where, reg(ri).name);
    elseif hasInline
        q = dd_parseprior(reg(ri).inline, n, sprintf('%s, prior of %s', where, reg(ri).name));
    elseif hasField
        q = dd_parseprior(pri.(reg(ri).name), n, sprintf('%s, priors.%s', where, reg(ri).name));
    else
        q = dd_parseprior(zeros(1, n), n, '');
        zeroPrior{end + 1} = reg(ri).name; %#ok<AGROW>
    end
    b0(reg(ri).idx) = q.mean;
    bdist(reg(ri).idx) = q.dist;
    ba(reg(ri).idx) = q.a;
    bb(reg(ri).idx) = q.b;
end

% random coefficients
rdist = zeros(P, 1);
rsign = ones(P, 1);
s0 = zeros(P, 1);
sdist = zeros(P, 1);
sa = zeros(P, 1);
sb = zeros(P, 1);
hasRandom = isfield(mdl, 'random') && ~isempty(mdl.random);
if strcmp(ty, 'mnl')
    if hasRandom
        dd_err('spec:random', '%s is of type mnl but lists random parameters. Set type to ''rp'' or ''rppanel''.', where);
    end
else
    if ~hasRandom || ~isstruct(mdl.random)
        dd_err('spec:random', '%s is of type %s and needs random: a struct parameter -> distribution, e.g. ''normal(0.3)''.', where, ty);
    end
    fb = fieldnames(mdl.random)';
    for q = 1:numel(fb)
        ri = find(strcmp(fb{q}, {reg.name}));
        if isempty(ri)
            dd_err('spec:random', '%s: random.%s does not match a parameter. Parameters: %s.', where, fb{q}, strjoin({reg.name}, ', '));
        end
        n = numel(reg(ri).idx);
        r = dd_parseprior(mdl.random.(fb{q}), n, sprintf('%s, random.%s', where, fb{q}), true);
        rdist(reg(ri).idx) = r.rdist;
        rsign(reg(ri).idx) = r.rsign;
        s0(reg(ri).idx) = r.mean;
        sdist(reg(ri).idx) = r.dist;
        sa(reg(ri).idx) = r.a;
        sb(reg(ri).idx) = r.b;
    end
end
ridx = find(rdist > 0);

Mm = struct();
Mm.name = char(getf(mdl, 'name', sprintf('M%d', m)));
Mm.type = ty;
Mm.J = J;
Mm.P = P;
Mm.parNames = parNames;
Mm.parBase = parBase;
Mm.T = numel(tAlt);
Mm.tAlt = tAlt;
Mm.tPar = tPar;
Mm.tConst = tConst;
Mm.tNF = tNF;
Mm.tCol = tCol;
Mm.tTab = tTab;
Mm.colTerms = cell(1, ctx.K);
for k = 1:ctx.K
    Mm.colTerms{k} = find(any(tCol == k, 2))';
end
Mm.b0 = b0;
Mm.bdist = bdist;
Mm.ba = ba;
Mm.bb = bb;
Mm.rdist = rdist;
Mm.rsign = rsign;
Mm.s0 = s0;
Mm.sdist = sdist;
Mm.sa = sa;
Mm.sb = sb;
Mm.ridx = ridx(:)';
Mm.Pr = numel(ridx);
Mm.Pt = P + numel(ridx);
sdNames = cellfun(@(s) ['sd(' s ')'], parNames(ridx), 'UniformOutput', false);
Mm.thNames = [parNames sdNames];
Mm.thBase = [parBase parBase(ridx)];
Mm.isBayes = any(bdist > 0) || any(sdist(ridx) > 0);
st = lower(char(getf(mdl, 'stat', 'mean')));
if ~any(strcmp(st, {'mean', 'median', 'fixed'}))
    dd_err('spec:stat', '%s: stat must be ''mean'', ''median'' or ''fixed'' (how Bayesian prior draws are summarised).', where);
end
if ~Mm.isBayes
    st = 'fixed';
end
Mm.stat = st;
Mm.reg = reg;
Mm.zeroPrior = zeroPrior;
Mm.eff = true(Mm.Pt, 1);
Mm.exMode = 'nuisance';
Mm.Peff = Mm.Pt;
end

% =============================================================================
function terms = parseUtility(str, where)
terms = struct('par', {}, 'coding', {}, 'prior', {}, 'attrs', {});
str = strtrim(char(str));
if isempty(str) || strcmp(str, '0')
    return
end
parts = splitTop(str, '+');
for i = 1:numel(parts)
    part = strtrim(parts{i});
    if isempty(part)
        dd_err('spec:utility', '%s: empty term in ''%s''.', where, str);
    end
    fac = strtrim(splitTop(part, '*'));
    head = fac{1};
    t = struct('par', '', 'coding', 'linear', 'prior', '', 'attrs', {{}});
    ib = find(head == '[', 1);
    okHead = true;
    if ~isempty(ib)
        okHead = head(end) == ']';
        if okHead
            t.prior = strtrim(head(ib + 1:end - 1));
            head = strtrim(head(1:ib - 1));
        end
    end
    tk = regexp(head, '^([A-Za-z]\w*)\.(dummy|effects|linear|d|e)$', 'tokens', 'once');
    if ~isempty(tk)
        t.par = tk{1};
        switch tk{2}
            case {'dummy', 'd'}, t.coding = 'dummy';
            case {'effects', 'e'}, t.coding = 'effects';
        end
    elseif isvarname(head)
        t.par = head;
    else
        okHead = false;
    end
    if ~okHead
        dd_err('spec:utility', ['%s: cannot read ''%s''. A term is parameter[.dummy|.effects][[prior]] * attribute * ..., ' ...
            'terms are joined by +, and negative effects are written as negative priors.'], where, part);
    end
    for f = 2:numel(fac)
        if ~isvarname(fac{f})
            dd_err('spec:utility', '%s: ''%s'' in term ''%s'' is not an attribute name.', where, fac{f}, part);
        end
        t.attrs{end + 1} = fac{f};
    end
    terms(end + 1) = t; %#ok<AGROW>
end
end

function parts = splitTop(str, delim)
depth = 0;
parts = {};
st = 1;
for i = 1:numel(str)
    ch = str(i);
    if ch == '[' || ch == '('
        depth = depth + 1;
    elseif ch == ']' || ch == ')'
        depth = depth - 1;
    elseif ch == delim && depth == 0
        parts{end + 1} = str(st:i - 1); %#ok<AGROW>
        st = i + 1;
    end
end
parts{end + 1} = str(st:end);
end

% =============================================================================
function f = compileExpr(expr, ctx, where)
expr = char(expr);
funcs = {'abs', 'min', 'max', 'mod', 'floor', 'ceil', 'round', 'and', 'or', 'not', 'xor', 'true', 'false'};
[tok, st, en] = regexp(expr, '(\d+\.?\d*([eE][-+]?\d+)?)|(\.\d+([eE][-+]?\d+)?)|([A-Za-z_]\w*(\.[A-Za-z_]\w*)?)', 'match', 'start', 'end');
out = expr;
nCol = 0;
for i = numel(tok):-1:1
    name = tok{i};
    if ~isempty(regexp(name(1), '[\d\.]', 'once'))
        continue
    end
    dot = find(name == '.', 1);
    rep = '';
    if ~isempty(dot)
        aN = name(1:dot - 1);
        tN = name(dot + 1:end);
        ja = find(strcmp(aN, ctx.alts));
        ia = find(strcmp(tN, ctx.attrNames));
        if isempty(ja) || isempty(ia)
            dd_err('spec:constraint', '%s: ''%s'' is not alternative.attribute. Alternatives: %s. Attributes: %s.', ...
                where, name, strjoin(ctx.alts, ', '), strjoin(ctx.attrNames, ', '));
        end
        if isfield(ctx.FX, aN) && isfield(ctx.FX.(aN), tN)
            rep = sprintf('(%.17g)', ctx.FX.(aN).(tN));
        elseif ctx.scenCol(ia) > 0 && any(strcmp(tN, ctx.scen))
            rep = sprintf('LV{%d}(D(:,%d))', ctx.scenCol(ia), ctx.scenCol(ia));
            nCol = nCol + 1;
        elseif ctx.colIndex(ja, ia) > 0
            rep = sprintf('LV{%d}(D(:,%d))', ctx.colIndex(ja, ia), ctx.colIndex(ja, ia));
            nCol = nCol + 1;
        else
            dd_err('spec:constraint', '%s: alternative %s has no attribute %s in the design.', where, aN, tN);
        end
    elseif any(strcmp(name, ctx.scen)) && ctx.scenCol(strcmp(name, ctx.attrNames)) > 0
        k = ctx.scenCol(strcmp(name, ctx.attrNames));
        rep = sprintf('LV{%d}(D(:,%d))', k, k);
        nCol = nCol + 1;
    elseif any(strcmp(name, funcs))
        continue
    else
        dd_err('spec:constraint', '%s: unknown name ''%s''. Refer to attributes as alternative.attribute (scenario variables by name).', where, name);
    end
    out = [out(1:st(i) - 1) rep out(en(i) + 1:end)];
end
rest = regexprep(expr, '(\d+\.?\d*([eE][-+]?\d+)?)|(\.\d+([eE][-+]?\d+)?)|([A-Za-z_]\w*(\.[A-Za-z_]\w*)?)', '');
if ~isempty(regexp(rest, '[^\s\(\)&\|~=<>\+\-\*/,]', 'once'))
    dd_err('spec:constraint', '%s: unsupported characters in ''%s''. Use MATLAB operators: & | ~ == ~= < <= > >= + - * /.', where, expr);
end
if nCol == 0
    dd_err('spec:constraint', '%s: ''%s'' does not refer to any design attribute.', where, expr);
end
try
    f = str2func(['@(D,LV) ' out]);
    r = f(ones(2, ctx.K), cellfun(@(v) v(:), ctx.levels, 'UniformOutput', false));
catch err
    dd_err('spec:constraint', '%s: cannot evaluate ''%s'' (%s).', where, expr, err.message);
end
if ~isequal(size(r), [2 1])
    dd_err('spec:constraint', '%s: ''%s'' must give one true/false value per choice task.', where, expr);
end
end

% =============================================================================
function [dom, notes] = compileDominance(opt, U1, C, ctx)
dom = struct('a', {}, 'b', {}, 'ka', {}, 'kb', {}, 'ua', {}, 'ub', {}, 'comp', {});
notes = {};
J = ctx.J;
explicit = {};
if islogical(opt) || isnumeric(opt)
    if isscalar(opt) && opt, mode = 'on'; else, mode = 'off'; end
elseif ischar(opt) || (isstring(opt) && isscalar(opt))
    mode = lower(char(opt));
    if ~any(strcmp(mode, {'auto', 'on', 'off'}))
        dd_err('spec:dominance', 'spec.dominance must be ''auto'', true, false or a list of alternatives.');
    end
else
    explicit = toCellstr(opt, 'spec.dominance');
    for i = 1:numel(explicit)
        mustBeAlt(explicit{i}, ctx.alts, 'spec.dominance');
    end
    mode = 'list';
end
if strcmp(mode, 'off')
    return
end
% structural signature of each alternative in model 1
sig = cell(1, J);
hasCol = any(ctx.colIndex > 0, 2)';
for j = 1:J
    s = cell(1, numel(U1{j}));
    for t = 1:numel(U1{j})
        s{t} = sprintf('%s.%s*%s', U1{j}(t).par, U1{j}(t).coding, strjoin(U1{j}(t).attrs, '*'));
    end
    sig{j} = strjoin(sort(s), '+');
end
pairs = zeros(0, 2);
for a = 1:J - 1
    for b = a + 1:J
        if ~hasCol(a) || ~hasCol(b), continue; end
        if strcmp(mode, 'list')
            take = any(strcmp(ctx.alts{a}, explicit)) && any(strcmp(ctx.alts{b}, explicit));
        else
            take = strcmp(sig{a}, sig{b}) && isequal(ctx.colIndex(a, :) > 0, ctx.colIndex(b, :) > 0);
        end
        if take
            pairs(end + 1, :) = [a b]; %#ok<AGROW>
        end
    end
end
if isempty(pairs)
    if strcmp(mode, 'on')
        dd_err('spec:dominance', 'Dominance checks were requested, but no two designed alternatives share the same attributes.');
    elseif strcmp(mode, 'list')
        notes{end + 1} = sprintf('Dominance: not checked; among %s no two alternatives have attributes that vary in the design.', strjoin(explicit, ', '));
    else
        notes{end + 1} = 'Dominance: not checked (no two designed alternatives have the same utility structure).';
    end
    return
end
Mm = C.model(1);
interact = false;
unknown = {};
for i = 1:size(pairs, 1)
    a = pairs(i, 1);
    b = pairs(i, 2);
    ia = find(ctx.colIndex(a, :) > 0 & ctx.colIndex(b, :) > 0);
    d = struct('a', a, 'b', b, 'ka', ctx.colIndex(a, ia), 'kb', ctx.colIndex(b, ia), 'ua', {{}}, 'ub', {{}}, 'comp', []);
    for q = 1:numel(ia)
        [d.ua{q}, c1, i1] = partialUtility(Mm, a, d.ka(q), C.L(d.ka(q)));
        [d.ub{q}, c2, i2] = partialUtility(Mm, b, d.kb(q), C.L(d.kb(q)));
        d.comp(q) = c1 || c2;
        interact = interact || i1 || i2;
        if ~d.comp(q)
            unknown{end + 1} = ctx.attrNames{ia(q)}; %#ok<AGROW>
        end
    end
    d.comp = logical(d.comp);
    dom(end + 1) = d; %#ok<AGROW>
end
pn = arrayfun(@(i) sprintf('%s-%s', ctx.alts{pairs(i, 1)}, ctx.alts{pairs(i, 2)}), 1:size(pairs, 1), 'UniformOutput', false);
notes{end + 1} = sprintf('Dominance: tasks with a dominated or identical alternative are excluded for %s (ranking from the prior means of model 1).', strjoin(pn, ', '));
if ~isempty(unknown)
    notes{end + 1} = sprintf('Dominance: zero priors give no ranking for %s; tasks that differ on these attributes are never treated as dominated.', strjoin(unique(unknown, 'stable'), ', '));
end
if interact
    notes{end + 1} = 'Dominance: interaction terms are ignored in the ranking (main effects only).';
end
end

function [u, comp, interact] = partialUtility(Mm, j, k, L)
u = zeros(L, 1);
comp = false;
interact = false;
for t = find(Mm.tAlt == j & any(Mm.tCol == k, 2)')
    if Mm.tNF(t) == 1
        u = u + Mm.b0(Mm.tPar(t)) * Mm.tConst(t) * Mm.tTab{t, 1};
        comp = comp || Mm.b0(Mm.tPar(t)) ~= 0;
    else
        interact = true;
    end
end
end

function mp = mirrorPairs(C)
% Pairs of alternatives whose columns can be exchanged to detect mirror-image
% duplicate tasks (same attributes with the same level sets).
mp = struct('ka', {}, 'kb', {});
for a = 1:C.J - 1
    for b = a + 1:C.J
        ia = find(C.colIndex(a, :) > 0 | C.colIndex(b, :) > 0);
        if isempty(ia) || ~isequal(C.colIndex(a, ia) > 0, C.colIndex(b, ia) > 0), continue; end
        ka = C.colIndex(a, ia);
        kb = C.colIndex(b, ia);
        same = all(arrayfun(@(q) isequal(C.levels{ka(q)}, C.levels{kb(q)}), 1:numel(ka)));
        if same && any(arrayfun(@(d) d.a == a && d.b == b, C.dom))
            mp(end + 1) = struct('ka', ka, 'kb', kb); %#ok<AGROW>
        end
    end
end
end

% =============================================================================
function checkFields(s, known, where)
f = fieldnames(s);
for i = 1:numel(f)
    if ~any(strcmp(f{i}, known))
        dd_err('spec:unknownField', '%s has an unknown field ''%s''. Known fields: %s.', where, f{i}, strjoin(known, ', '));
    end
end
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function c = toCellstr(v, where)
if ischar(v)
    c = {v};
elseif isstring(v)
    c = cellstr(v(:))';
elseif iscellstr(v) %#ok<ISCLSTR>
    c = v(:)';
elseif iscell(v) && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), v(:)))
    c = cellfun(@char, v(:)', 'UniformOutput', false);
else
    dd_err('spec:type', '%s must be text or a list of text values.', where);
end
c = strtrim(c);
end

function n = reqInt(spec, name, minVal)
if ~isfield(spec, name) || isempty(spec.(name))
    dd_err(['spec:' name], 'spec.%s is required.', name);
end
n = double(spec.(name));
if ~isscalar(n) || ~isfinite(n) || n ~= round(n) || n < minVal
    dd_err(['spec:' name], 'spec.%s must be an integer >= %d.', name, minVal);
end
end

function v = levelVector(x, where, minN)
if ~isnumeric(x) && ~islogical(x)
    dd_err('spec:levels', '%s must be a numeric vector of levels (use spec.labels for text labels).', where);
end
v = double(x(:)');
if numel(v) < minN || any(~isfinite(v))
    dd_err('spec:levels', '%s needs at least %d finite levels.', where, minN);
end
if numel(unique(v)) < numel(v)
    dd_err('spec:levels', '%s has repeated levels.', where);
end
end

function mustBeAlt(name, alts, where)
if ~any(strcmp(name, alts))
    dd_err('spec:alternatives', '%s refers to ''%s'', which is not an alternative (%s).', where, name, strjoin(alts, ', '));
end
end

function mustBeAttr(name, attrNames, where)
if ~any(strcmp(name, attrNames))
    dd_err('spec:attributes', '%s refers to ''%s'', which is not an attribute (%s).', where, name, strjoin(attrNames, ', '));
end
end

function e = expandFreq(x, L, attr, what)
e = double(x(:));
if isscalar(e)
    e = e * ones(L, 1);
end
if numel(e) ~= L || any(e < 0) || any(e ~= round(e))
    dd_err('spec:frequency', 'spec.frequency.%s.%s must be one non-negative integer, or one per level (%d).', attr, what, L);
end
end
