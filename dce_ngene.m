function [spec, notes] = dce_ngene(src, varargin)
%DCE_NGENE Translate an Ngene script into a design specification.
%
%   spec = DCE_NGENE(src) reads an Ngene script (the text itself, an .ngs file,
%   or an .ngd file, which carries its script) and returns the specification
%   that dce_design and dce_evaluate use. notes (also stored in spec.notes)
%   lists every place where the translation had to decide something.
%
%   Supported: ;alts (with * for the dominance and repetition checks), ;rows,
%   ;block, ;eff with one or several weighted models of type mnl, rp or
%   rppanel and measure d or a (statistic fixed, mean or median), ;model with
%   generic and alternative-specific parameters, constants, fixed, Bayesian
%   ((n,..), (u,..)) and random (n,.. and u,..) priors, dummy and effects
%   coding, interactions, level lists, level frequencies such as (8,8,8) or
%   (2-4,2-4), scenario variables x[x], single-level (status quo) attributes,
%   ;cond, ;reject, ;require, ;con, ;bdraws, ;rdraws, ;rep, ;start and stop
%   rules of ;alg. Anything else raises an error that names it, so nothing is
%   dropped silently.
%
%   Where Ngene's conventions differ from this toolbox, the translation uses
%   the toolbox's and says so in the notes:
%     - constants are part of the efficiency measure (Ngene leaves them out
%       unless ;con is given);
%     - Bayesian priors are integrated over (Ngene evaluates them at their
%       means unless ;eff asks for mean or median);
%     - in panel models a respondent answers one block (Ngene assumes every
%       row);
%     - several models are each scaled by their own optimum and weighted
%       equally (in Ngene the weights of ;eff multiply raw D-errors and are
%       usually chosen by hand to equalise their scales).
%   spec = DCE_NGENE(src, 'Emulate', true) keeps Ngene's four conventions
%   instead, so that values can be compared with the ones Ngene reported.
%   In both cases the base level of a dummy or effects coded attribute is the
%   last level listed, as in Ngene.
%
%   See also dce_design, dce_evaluate, dce_readngd.

p = inputParser;
p.FunctionName = 'dce_ngene';
addParameter(p, 'Emulate', false, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
emulate = logical(p.Results.Emulate);

src = char(src);
folder = '';
name = 'ngene';
if size(src, 1) == 1 && numel(src) < 1000 && isfile(src)
    [folder, name, ext] = fileparts(src);
    if strcmpi(ext, '.ngd')
        [~, txt] = dce_readngd(src);
        if isempty(txt)
            error('dce:ngene:format', '%s contains no script.', src);
        end
    else
        txt = fileread(src);
    end
else
    txt = src;
end
notes = {};
txt = strrep(txt, char(13), '');
txt = strrep(txt, char(65279), '');
txt = regexprep(txt, '(^|\n)\s*//[^\n]*', '$1');
txt = regexprep(txt, '\?[^\n]*', '');
cut = find(txt == '$', 1);
if ~isempty(cut)
    txt = txt(1:cut - 1);
end
txt = regexprep(strtrim(txt), '^design\s*', '', 'once', 'ignorecase');

% ---- properties ---------------------------------------------------------------
P = struct('name', {}, 'model', {}, 'value', {});
parts = splitTop(txt, ';');
for i = 1:numel(parts)
    s = strtrim(parts{i});
    if isempty(s)
        continue
    end
    tk = regexp(s, '^[A-Za-z]+', 'match', 'once');
    if isempty(tk)
        error('dce:ngene:syntax', 'Cannot read the property starting with ''%s''.', s(1:min(30, end)));
    end
    rest = strtrim(s(numel(tk) + 1:end));
    mdl = '';
    if ~isempty(rest) && rest(1) == '('
        close = matchParen(rest, 1);
        mdl = strtrim(rest(2:close - 1));
        rest = strtrim(rest(close + 1:end));
    end
    if ~isempty(rest) && any(rest(1) == '=:')
        rest = strtrim(rest(2:end));
    end
    P(end + 1) = struct('name', lower(tk), 'model', mdl, 'value', rest); %#ok<AGROW>
end
names = {P.name};
unsupported = intersect(unique(names), {'fisher', 'foldover', 'gfim'});
if ~isempty(unsupported)
    error('dce:ngene:unsupported', 'Ngene properties not supported by this toolbox: %s.', strjoin(strcat(';', unsupported), ', '));
end
ignored = intersect(unique(names), {'store', 'prec', 'bseed', 'rseed', 'trimdist', 'wtp', 'eval'});
ignored = [ignored unique(names(strncmp(names, 'format', 6)))];
if ~isempty(ignored)
    notes{end + 1} = sprintf(['Ignored Ngene properties: %s (seeds are set with the Seed option; ' ...
        'an existing design is evaluated with dce_evaluate).'], strjoin(strcat(';', ignored), ', '));
end
known = {'alts', 'rows', 'block', 'eff', 'alg', 'model', 'cond', 'reject', 'require', 'con', 'bdraws', 'rdraws', 'rep', 'start', ...
    'fact', 'rand', 'orth'};
other = setdiff(unique(names), [known unsupported ignored]);
if ~isempty(other)
    error('dce:ngene:unsupported', 'Unknown Ngene properties: %s.', strjoin(strcat(';', other), ', '));
end

spec = struct();
spec.name = name;

% ---- alternatives ----------------------------------------------------------------
ia = find(strcmp(names, 'alts'));
if isempty(ia)
    error('dce:ngene:syntax', ';alts is missing.');
end
alts = {};
starred = {};
for i = ia
    v = P(i).value;
    if any(v == '(')
        error('dce:ngene:unsupported', 'Grouped alternatives in ;alts (Ngene Online syntax) are not supported; mark them with *.');
    end
    a = strtrim(strsplit(v, ','));
    st = cellfun(@(s) ~isempty(s) && s(end) == '*', a);
    a = strtrim(strrep(a, '*', ''));
    [~, first] = unique(lower(a), 'stable');
    if numel(first) < numel(a)
        notes{end + 1} = sprintf(';alts lists an alternative more than once (%s); the repeats are dropped.', strjoin(a, ', ')); %#ok<AGROW>
        a = a(first);
        st = st(first);
    end
    if isempty(alts)
        alts = a;
    elseif ~isequal(sort(alts), sort(a))
        error('dce:ngene:unsupported', 'The models list different alternatives in ;alts; this is not supported.');
    end
    starred = union(starred, a(st), 'stable');
end
spec.alternatives = alts;

% ---- rows, blocks ------------------------------------------------------------------
ir = find(strcmp(names, 'rows'), 1);
if isempty(ir)
    error('dce:ngene:syntax', ';rows is missing.');
end
rv = str2double(P(ir).value);
isFull = any(strcmp(names, 'fact')) || any(strcmp(names, 'rand'));
if isFull
    spec.type = 'full';
    if ~isnan(rv)
        spec.rows = rv;
        notes{end + 1} = sprintf('Full factorial with ;rows = %d: a random fraction of that size.', rv);
    elseif ~strcmpi(strtrim(P(ir).value), 'all')
        error('dce:ngene:syntax', 'Cannot read ;rows = %s.', P(ir).value);
    end
else
    if isnan(rv)
        error('dce:ngene:unsupported', ';rows = %s is not supported (give the number of choice situations).', P(ir).value);
    end
    spec.rows = rv;
end
io = find(strcmp(names, 'orth'), 1);
if ~isempty(io)
    ov = lower(strtrim(P(io).value));
    if isempty(ov) || strcmp(ov, 'sim')
        spec.type = 'orthogonal';
    elseif any(strcmp(ov, {'seq', 'seq2'}))
        spec.type = 'orthogonal-sequential';
    else
        error('dce:ngene:unsupported', ';orth = %s is not supported (sim and seq are).', P(io).value);
    end
end
ib = find(strcmp(names, 'block'), 1);
if ~isempty(ib)
    bv = sscanf(P(ib).value, '%f', 1);
    if ~isempty(bv) && bv > 1
        spec.blocks = bv;
    end
    if any(P(ib).value == ',')
        notes{end + 1} = 'Blocking options of ;block are replaced by the balanced blocking of dce_block.';
    end
end

% ---- efficiency: models, weights, statistics ---------------------------------------
ie = find(strcmp(names, 'eff'), 1);
E = struct('w', {}, 'name', {}, 'type', {}, 'measure', {}, 'stat', {});
if isempty(ie)
    E(1) = struct('w', 1, 'name', '', 'type', 'mnl', 'measure', 'd', 'stat', 'fixed');
else
    terms = splitTop(P(ie).value, '+');
    for i = 1:numel(terms)
        t = strtrim(terms{i});
        op = find(t == '(', 1);
        if isempty(op) || t(end) ~= ')'
            error('dce:ngene:syntax', 'Cannot read the ;eff term ''%s''.', t);
        end
        before = strtrim(t(1:op - 1));
        w = 1;
        nm = before;
        star = find(before == '*', 1);
        if ~isempty(star)
            w = str2double(before(1:star - 1));
            nm = strtrim(before(star + 1:end));
            if isnan(w)
                error('dce:ngene:syntax', 'Cannot read the weight in the ;eff term ''%s''.', t);
            end
        end
        arg = strtrim(strsplit(t(op + 1:end - 1), ','));
        ty = lower(arg{1});
        if strcmp(ty, 'imbalance')
            error('dce:ngene:unsupported', 'The ;eff term ''%s'' is not supported (use spec.frequency for level balance).', t);
        end
        if isscalar(arg)
            arg{2} = 'd';
        end
        if ~any(strcmp(ty, {'mnl', 'rp', 'rppanel'}))
            error('dce:ngene:unsupported', 'Model type ''%s'' in ;eff is not supported (mnl, rp and rppanel are).', ty);
        end
        me = lower(arg{2});
        if ~any(strcmp(me, {'d', 'a'}))
            error('dce:ngene:unsupported', 'Efficiency measure ''%s'' in ;eff is not supported (d and a are).', me);
        end
        stt = 'fixed';
        if numel(arg) >= 3
            stt = lower(arg{3});
            if ~any(strcmp(stt, {'fixed', 'mean', 'median'}))
                error('dce:ngene:unsupported', 'Statistic ''%s'' in ;eff is not supported (fixed, mean and median are).', stt);
            end
        end
        if numel(arg) >= 4 && ~strcmpi(arg{4}, 'rum')
            error('dce:ngene:unsupported', 'Decision rule ''%s'' in ;eff is not supported.', arg{4});
        end
        E(end + 1) = struct('w', w, 'name', nm, 'type', ty, 'measure', me, 'stat', stt); %#ok<AGROW>
    end
end
if numel(unique({E.measure})) > 1
    error('dce:ngene:unsupported', 'Mixing D- and A-error in one ;eff is not supported.');
end
spec.criterion = E(1).measure;

% ---- models --------------------------------------------------------------------------
im = find(strcmp(names, 'model'));
if isempty(im)
    error('dce:ngene:syntax', ';model is missing.');
end
A = struct('name', {}, 'defs', {}, 'defAlt', {}, 'freq', {}, 'scenario', {}, 'coded', {});
constants = {};
models = cell(1, numel(E));
for e = 1:numel(E)
    if isempty(E(e).name)
        hit = im(cellfun(@isempty, {P(im).model}));
        if isempty(hit) && isscalar(im)
            hit = im;
        elseif isempty(hit) && isempty(ie)
            % no ;eff (full factorial, orthogonal): the first model describes the attributes
            hit = im(1);
            if numel(im) > 1
                notes{end + 1} = sprintf('No ;eff: of the %d ;model blocks the first (%s) is used.', numel(im), P(im(1)).model); %#ok<AGROW>
            end
        end
    else
        hit = im(strcmpi({P(im).model}, E(e).name));
    end
    if numel(hit) ~= 1
        error('dce:ngene:syntax', 'Cannot match the ;eff term %d to exactly one ;model block.', e);
    end
    [models{e}, A, cst] = parseModel(P(hit).value, alts, A, E(e));
    constants = union(constants, cst, 'stable');
    if isempty(E(e).name)
        models{e}.name = upper(E(e).type);
    else
        models{e}.name = E(e).name;
    end
end

% ---- attributes ------------------------------------------------------------------------
spec.attributes = struct();
scen = {};
for i = 1:numel(A)
    a = A(i);
    multi = find(cellfun(@numel, a.defs) > 1);
    if isempty(multi)
        if isempty(a.defs)
            error('dce:ngene:syntax', 'Attribute %s has no levels.', a.name);
        end
        spec.attributes.(a.name) = a.defs{1};
        for q = 2:numel(a.defs)
            if ~isequal(a.defs{q}, a.defs{1})
                spec.fixed.(a.defAlt{q}).(a.name) = a.defs{q};
            end
        end
    else
        glob = a.defs{multi(1)};
        spec.attributes.(a.name) = glob;
        for q = 1:numel(a.defs)
            if isscalar(a.defs{q})
                spec.fixed.(a.defAlt{q}).(a.name) = a.defs{q};
            elseif ~isequal(a.defs{q}, glob)
                spec.altLevels.(a.defAlt{q}).(a.name) = a.defs{q};
            end
        end
        if a.coded
            spec.base.(a.name) = glob(end);
        end
        if ~isempty(a.freq)
            if isequal(a.freq.min, a.freq.max)
                spec.frequency.(a.name).exact = a.freq.min;
            else
                spec.frequency.(a.name).min = a.freq.min;
                spec.frequency.(a.name).max = a.freq.max;
            end
        end
    end
    if a.scenario
        scen{end + 1} = a.name; %#ok<AGROW>
    end
end
if ~isempty(scen)
    spec.scenario = scen;
end

% ---- assemble the models -----------------------------------------------------------------
for e = 1:numel(E)
    m = struct('name', models{e}.name, 'type', E(e).type, 'stat', E(e).stat, 'utility', models{e}.utility, ...
        'priors', models{e}.priors, 'random', models{e}.random);
    if isempty(fieldnames(m.random))
        m = rmfield(m, 'random');
        if ~strcmp(E(e).type, 'mnl')
            error('dce:ngene:syntax', 'Model %s is of type %s but has no random parameters.', m.name, E(e).type);
        end
    elseif strcmp(E(e).type, 'mnl')
        m = rmfield(m, 'random');
        notes{end + 1} = sprintf('Model %s is evaluated as mnl: random parameters enter with their means.', m.name); %#ok<AGROW>
    end
    if models{e}.bayes && strcmp(E(e).stat, 'fixed')
        if emulate
            notes{end + 1} = sprintf(['Model %s has Bayesian priors, but ;eff does not ask for mean or median: as in Ngene, ' ...
                'the design is evaluated at the prior means only.'], m.name); %#ok<AGROW>
        else
            m.stat = 'mean';
            notes{end + 1} = sprintf(['Model %s has Bayesian priors, but ;eff does not ask for mean or median. Ngene would then ' ...
                'use only the prior means; here the D-error is averaged over the prior distributions.'], m.name); %#ok<AGROW>
        end
    end
    models{e} = m;
end
if numel(E) == 1
    spec.model = models{1};
else
    spec.models = models;
    if emulate
        spec.weights = [E.w];
        spec.search.Normalize = false;
        notes{end + 1} = 'The weights of ;eff multiply the raw D-errors, as in Ngene.';
    else
        spec.weights = ones(1, numel(E));
        notes{end + 1} = sprintf(['The weights of ;eff (%s) multiply raw D-errors in Ngene. Here each model is scaled by its own ' ...
            'optimum and the models are weighted equally; set spec.weights for another balance.'], strtrim(sprintf('%g ', [E.w])));
    end
end

% ---- constants in the criterion ------------------------------------------------------------
if ~any(strcmp(names, 'con')) && ~isempty(constants)
    if emulate
        spec.constants = false;
        notes{end + 1} = sprintf('No ;con: as in Ngene, the constants (%s) are estimated but left out of the efficiency measure.', strjoin(constants, ', '));
    else
        notes{end + 1} = sprintf(['No ;con: Ngene would leave the constants (%s) out of the efficiency measure; here they are part of it ' ...
            '(set spec.constants = false to leave them out).'], strjoin(constants, ', '));
    end
end

% ---- constraints -----------------------------------------------------------------------------
rej = {};
ic = find(strcmp(names, 'cond'));
for i = ic
    items = splitTop(P(i).value, ',');
    for q = 1:numel(items)
        it = strtrim(regexprep(items{q}, 'fractional\s*=\s*[\d\.]+\s*%?', '', 'ignorecase'));
        if isempty(it)
            continue
        end
        tk = regexpi(it, '^if\s*\((.*)\)$', 'tokens', 'once');
        if isempty(tk)
            error('dce:ngene:syntax', 'Cannot read the ;cond item ''%s''.', it);
        end
        ab = splitTop(tk{1}, ',');
        if numel(ab) ~= 2
            error('dce:ngene:syntax', 'A ;cond item needs if(condition, consequence): ''%s''.', it);
        end
        rej{end + 1} = sprintf('(%s) & ~(%s)', convExpr(ab{1}), convExpr(ab{2})); %#ok<AGROW>
        % a consequence that restricts one attribute to some of its levels nests
        % that attribute in another; as in Ngene it then loses exact level balance
        nest = regexp(strtrim(ab{2}), '^\w+\.(\w+)\s*(?:=|<>|<=|>=|<|>)\s*(?:\[[^\]]*\]|[-+]?[\d\.]+)$', 'tokens', 'once');
        if ~isempty(nest)
            an = fieldnames(spec.attributes);
            hit = find(strcmpi(nest{1}, an), 1);
            if ~isempty(hit) && numel(spec.attributes.(an{hit})) > 1 && ~(isfield(spec, 'frequency') && isfield(spec.frequency, an{hit}))
                bal = spec.rows / numel(spec.attributes.(an{hit}));
                spec.frequency.(an{hit}) = struct('min', floor(0.5 * bal), 'max', ceil(1.5 * bal));
                notes{end + 1} = sprintf([';cond restricts %s to some of its levels depending on another attribute; ' ...
                    'its level counts may vary between %d and %d instead of being exactly balanced.'], ...
                    an{hit}, floor(0.5 * bal), ceil(1.5 * bal)); %#ok<AGROW>
            end
        end
    end
end
for i = find(strcmp(names, 'reject'))
    items = splitTop(P(i).value, ',');
    for q = 1:numel(items)
        if ~isempty(strtrim(items{q}))
            rej{end + 1} = convExpr(items{q}); %#ok<AGROW>
        end
    end
end
attrAll = fieldnames(spec.attributes)';
if ~isempty(rej)
    spec.reject = cellfun(@(e) canonNames(e, alts, attrAll), rej, 'UniformOutput', false);
end
req = {};
for i = find(strcmp(names, 'require'))
    items = splitTop(P(i).value, ',');
    for q = 1:numel(items)
        if ~isempty(strtrim(items{q}))
            req{end + 1} = convExpr(items{q}); %#ok<AGROW>
        end
    end
end
if ~isempty(req)
    spec.require = cellfun(@(e) canonNames(e, alts, attrAll), req, 'UniformOutput', false);
end
if isempty(starred)
    spec.dominance = false;
    if emulate
        spec.mirror = false;
    end
else
    spec.dominance = starred(:)';
end

% ---- algorithm and draws -------------------------------------------------------------------------
il = find(strcmp(names, 'alg'), 1);
if ~isempty(il)
    v = lower(P(il).value);
    if contains(v, 'neldermead') || ~isempty(regexp(v, '^\s*all\b', 'once')) || contains(v, 'eval(')
        error('dce:ngene:unsupported', ';alg = %s is not supported.', P(il).value);
    end
    if contains(v, 'mfed')
        spec.balance = 'none';
        notes{end + 1} = ';alg = mfederov does not keep level balance: spec.balance is set to ''none''. Add spec.frequency bounds if needed.';
    end
    tk = regexp(v, 'total\s*\(\s*([\d\.]+)\s*(secs?|mins?|iterations?)', 'tokens', 'once');
    if ~isempty(tk) && ~strncmp(tk{2}, 'iter', 4)
        tmax = str2double(tk{1});
        if strncmp(tk{2}, 'min', 3)
            tmax = tmax * 60;
        end
        spec.search.MaxTime = tmax;
    end
    notes{end + 1} = 'The search settings of ;alg are replaced by the multi-start swap search of this toolbox (it stops by itself).';
end
[spec, notes] = drawsOption(spec, notes, P, names, 'bdraws', 'BayesDraws');
[spec, notes] = drawsOption(spec, notes, P, names, 'rdraws', 'RandomDraws');
ip = find(strcmp(names, 'rep'), 1);
if ~isempty(ip)
    % fewer than 1000 simulated respondents make the panel criterion too noisy to optimise
    spec.search.EvalPanelRespondents = max(2000, str2double(P(ip).value));
end
if any(strcmp({E.type}, 'rppanel'))
    if emulate
        spec.panelUnit = 'design';
        notes{end + 1} = ['Panel mixed logit: as in Ngene, a respondent is assumed to answer every row. ' ...
            'Set spec.panelUnit = ''block'' to evaluate what respondents actually see.'];
    elseif isfield(spec, 'blocks')
        notes{end + 1} = 'Panel mixed logit: a respondent answers one block (Ngene assumes that a respondent answers every row).';
    end
end
is = find(strcmp(names, 'start'), 1);
if ~isempty(is)
    sf = strtrim(P(is).value);
    cand = {sf, fullfile(folder, sf)};
    hit = cand(cellfun(@isfile, cand));
    if isempty(hit)
        notes{end + 1} = sprintf('The starting design ''%s'' was not found and is not used.', sf);
    else
        spec.start = hit{1};
    end
end
spec.notes = notes(:);
end

% =============================================================================
function [mdl, A, constants] = parseModel(body, alts, A, E)
body = regexprep(body, '\s+', ' ');
mdl = struct('utility', struct(), 'priors', struct(), 'random', struct(), 'bayes', false);
constants = {};
seen = struct();
utils = splitTop(body, '/');
for u = 1:numel(utils)
    s = strtrim(utils{u});
    if isempty(s)
        continue
    end
    tk = regexp(s, '^U\s*\(\s*(\w+)\s*\)\s*=\s*(.*)$', 'tokens', 'once', 'ignorecase');
    if isempty(tk)
        error('dce:ngene:syntax', 'Cannot read the utility function ''%s''.', s(1:min(40, end)));
    end
    alt = tk{1};
    ja = find(strcmpi(alt, alts), 1);
    if isempty(ja)
        error('dce:ngene:syntax', 'U(%s): %s is not listed in ;alts.', alt, alt);
    end
    alt = alts{ja};
    terms = splitTop(tk{2}, '+');
    out = cell(1, numel(terms));
    for t = 1:numel(terms)
        fac = strtrim(splitTop(strtrim(terms{t}), '*'));
        % --- parameter
        head = fac{1};
        prior = '';
        % the coding suffix may also follow the prior: b[0.1|0.2].dummy
        tk0 = regexp(head, '^(\w+)\s*\[(.*)\]\s*\.(dummy|effects|d|e)$', 'tokens', 'once', 'ignorecase');
        if ~isempty(tk0)
            head = [tk0{1} '.' tk0{3}];
            prior = tk0{2};
        else
            ob = find(head == '[', 1);
            if ~isempty(ob)
                if head(end) ~= ']'
                    error('dce:ngene:syntax', 'Cannot read the parameter ''%s''.', head);
                end
                prior = head(ob + 1:end - 1);
                head = strtrim(head(1:ob - 1));
            end
        end
        coding = '';
        tk2 = regexp(head, '^(\w+)\.(dummy|effects|d|e)$', 'tokens', 'once', 'ignorecase');
        if ~isempty(tk2)
            par = tk2{1};
            if lower(tk2{2}(1)) == 'd', coding = 'dummy'; else, coding = 'effects'; end
        elseif ~isempty(regexp(head, '^\w+$', 'once'))
            par = head;
        else
            error('dce:ngene:unsupported', 'Parameter syntax ''%s'' is not supported.', fac{1});
        end
        % Ngene is not case sensitive: keep the first spelling of each name
        key = ['p_' lower(par)];
        if isfield(seen, key)
            par = seen.(key).name;
            if isempty(coding), coding = seen.(key).coding; end
        else
            if ~isvarname(par)
                error('dce:ngene:unsupported', 'Parameter name ''%s'' is not a valid identifier here (start it with a letter).', par);
            end
            seen.(key) = struct('name', par, 'coding', coding);
        end
        if ~isempty(prior)
            [pv, rv, isB] = convPrior(prior, par);
            mdl.priors.(par) = pv;
            if ~isempty(rv)
                mdl.random.(par) = rv;
            end
            mdl.bayes = mdl.bayes || isB;
        elseif ~isfield(mdl.priors, par) && numel(fac) == 1
            mdl.priors.(par) = 0;
        end
        % --- attributes
        txt = par;
        if ~isempty(coding)
            txt = [par '.' coding];
        end
        if numel(fac) == 1
            constants{end + 1} = par; %#ok<AGROW>
        end
        for f = 2:numel(fac)
            [an, lev, frq, isCopy] = parseAttr(fac{f});
            ix = find(strcmpi(an, {A.name}), 1);
            if isempty(ix)
                if ~isvarname(an)
                    error('dce:ngene:unsupported', 'Attribute name ''%s'' is not a valid identifier here (start it with a letter).', an);
                end
                A(end + 1) = struct('name', an, 'defs', {{}}, 'defAlt', {{}}, 'freq', [], 'scenario', false, 'coded', false); %#ok<AGROW>
                ix = numel(A);
            end
            an = A(ix).name;
            if isCopy
                A(ix).scenario = true;
            elseif ~isempty(lev)
                A(ix).defs{end + 1} = lev;
                A(ix).defAlt{end + 1} = alt;
                if ~isempty(frq) && isempty(A(ix).freq)
                    A(ix).freq = frq;
                end
            end
            if f == 2 && ~isempty(coding)
                A(ix).coded = true;
            end
            txt = [txt '*' an]; %#ok<AGROW>
        end
        out{t} = txt;
    end
    mdl.utility.(alt) = strjoin(out, ' + ');
end
% parameters that never received a prior default to zero in the compiler
if strcmp(E.type, 'mnl')
    mdl.random = struct();
end
end

function [an, lev, frq, isCopy] = parseAttr(s)
s = strtrim(s);
lev = [];
frq = [];
isCopy = false;
tk = regexp(s, '^(\w+)', 'tokens', 'once');
if isempty(tk)
    error('dce:ngene:syntax', 'Cannot read the attribute ''%s''.', s);
end
an = tk{1};
rest = strtrim(s(numel(an) + 1:end));
if isempty(rest)
    return
end
if rest(1) == '.'
    error('dce:ngene:unsupported', 'Attribute modifier in ''%s'' (.dummy[level], .ref, .piv, .covar) is not supported.', s);
end
if rest(1) ~= '['
    error('dce:ngene:syntax', 'Cannot read the attribute ''%s''.', s);
end
cb = find(rest == ']', 1);
inner = strtrim(rest(2:cb - 1));
rest = strtrim(rest(cb + 1:end));
if strcmpi(inner, an)
    isCopy = true;
elseif ~isempty(regexp(inner, 'fcn\s*\(', 'once'))
    error('dce:ngene:unsupported', 'Attribute functions (fcn) are not supported: ''%s''.', s);
elseif any(inner == ':')
    pr = str2double(strsplit(inner, ':'));
    if numel(pr) ~= 3 || any(isnan(pr))
        error('dce:ngene:unsupported', 'Continuous attribute levels ''%s'' are not supported.', s);
    end
    lev = pr(1):pr(3):pr(2) + 1e-9 * pr(3);
elseif ~isempty(regexp(inner, '^[A-Za-z_]\w*$', 'once'))
    error('dce:ngene:unsupported', 'Attribute ''%s'' copies another attribute; only the form x[x] (scenario variable) is supported.', s);
else
    lev = str2double(strtrim(strsplit(inner, ',')));
    if any(isnan(lev))
        error('dce:ngene:syntax', 'Cannot read the levels in ''%s''.', s);
    end
end
if ~isempty(rest)
    if rest(1) ~= '(' || rest(end) ~= ')'
        error('dce:ngene:syntax', 'Cannot read the level frequencies in ''%s''.', s);
    end
    items = strtrim(strsplit(rest(2:end - 1), ','));
    mn = zeros(1, numel(items));
    mx = mn;
    for i = 1:numel(items)
        r = regexp(items{i}, '^(\d+)\s*-\s*(\d+|inf)$', 'tokens', 'once', 'ignorecase');
        if isempty(r)
            mn(i) = str2double(items{i});
            mx(i) = mn(i);
        else
            mn(i) = str2double(r{1});
            if strcmpi(r{2}, 'inf'), mx(i) = 1e9; else, mx(i) = str2double(r{2}); end
        end
    end
    if any(isnan(mn)) || numel(mn) ~= numel(lev)
        error('dce:ngene:syntax', 'The level frequencies in ''%s'' must give one count or range per level.', s);
    end
    frq = struct('min', mn, 'max', mx);
end
end

function [pv, rv, isBayes] = convPrior(txt, par)
% Ngene prior text -> priors value and random specification of the toolbox.
ent = strtrim(splitTop(txt, '|'));
n = numel(ent);
pv = cell(1, n);
rv = cell(1, n);
isRand = false(1, n);
isBayes = false;
for i = 1:n
    e = ent{i};
    if ~isempty(regexp(e, '^ec\b', 'once', 'ignorecase'))
        error('dce:ngene:unsupported', 'Error components (parameter %s) are not supported.', par);
    end
    if any(strcmp(e, {'+', '-'}))
        error('dce:ngene:unsupported', 'Sign-only priors (parameter %s) are not supported; give a small value with that sign.', par);
    end
    tk = regexp(e, '^(n|u)\s*,(.*)$', 'tokens', 'once', 'ignorecase');
    if ~isempty(tk)
        a = strtrim(splitTop(tk{2}, ','));
        if numel(a) ~= 2
            error('dce:ngene:syntax', 'Cannot read the random parameter prior ''%s'' of %s.', e, par);
        end
        isRand(i) = true;
        if lower(tk{1}) == 'n'
            pv{i} = oneValue(a{1});
            rv{i} = sprintf('normal(%s)', a{2});
            if str2double(a{2}) == 0
                rv{i} = 'fixed';
                isRand(i) = false;
            end
            isBayes = isBayes || any(a{1} == '(') || any(a{2} == '(');
        else
            lo = str2double(a{1});
            hi = str2double(a{2});
            if isnan(lo) || isnan(hi)
                error('dce:ngene:unsupported', 'A uniform random parameter with Bayesian bounds (%s) is not supported.', par);
            end
            pv{i} = (lo + hi) / 2;
            rv{i} = sprintf('uniform(%.15g)', (hi - lo) / 2);
        end
    else
        pv{i} = oneValue(e);
        rv{i} = 'fixed';
        isBayes = isBayes || any(e == '(');
    end
end
if all(cellfun(@isnumeric, pv))
    pv = [pv{:}];
end
if any(isRand)
    if n == 1
        rv = rv{1};
    end
else
    rv = [];
end
end

function v = oneValue(e)
x = str2double(e);
if ~isnan(x)
    v = x;
else
    v = strtrim(e);
end
end

function out = convExpr(s)
% Ngene logical expression -> MATLAB expression understood by the toolbox.
out = strtrim(regexprep(s, '\s+', ' '));
% sets: x = [a,b] and x <> [a,b]
while true
    [tk, st, en] = regexp(out, '([\w\.]+)\s*(=|<>)\s*\[([^\]]*)\]', 'tokens', 'start', 'end', 'once');
    if isempty(tk)
        break
    end
    vals = strtrim(strsplit(tk{3}, ','));
    if strcmp(tk{2}, '=')
        rep = ['(' strjoin(cellfun(@(v) sprintf('%s == %s', tk{1}, v), vals, 'UniformOutput', false), ' | ') ')'];
    else
        rep = ['(' strjoin(cellfun(@(v) sprintf('%s ~= %s', tk{1}, v), vals, 'UniformOutput', false), ' & ') ')'];
    end
    out = [out(1:st - 1) rep out(en + 1:end)];
end
out = strrep(out, '<>', '~=');
out = regexprep(out, '(?<![<>=~])=(?!=)', '==');
out = regexprep(out, '\<and\>', '&', 'ignorecase');
out = regexprep(out, '\<or\>', '|', 'ignorecase');
end

function out = canonNames(expr, alts, attrs)
% Ngene is not case sensitive: write alternative.attribute as declared.
out = expr;
[tok, st, en] = regexp(expr, '([A-Za-z_]\w*)\.([A-Za-z_]\w*)', 'tokens', 'start', 'end');
for i = numel(tok):-1:1
    ja = find(strcmpi(tok{i}{1}, alts), 1);
    ia = find(strcmpi(tok{i}{2}, attrs), 1);
    if ~isempty(ja) && ~isempty(ia)
        out = [out(1:st(i) - 1) alts{ja} '.' attrs{ia} out(en(i) + 1:end)];
    end
end
end

function [spec, notes] = drawsOption(spec, notes, P, names, prop, opt)
i = find(strcmp(names, prop), 1);
if isempty(i)
    return
end
tk = regexp(lower(P(i).value), '^(\w+)\s*\(\s*([\d,\s]+)\)', 'tokens', 'once');
if isempty(tk)
    return
end
n = prod(str2double(strsplit(tk{2}, ',')));
% the search needs far fewer draws than the reported values, which use the
% evaluation draws of the toolbox
spec.search.(opt) = max(128, min(n, 256));
if any(strcmp(tk{1}, {'sobol', 'halton'}))
    spec.search.DrawType = tk{1};
else
    notes{end + 1} = sprintf(';%s = %s: scrambled Sobol draws are used instead.', prop, P(i).value);
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

function close = matchParen(s, open)
depth = 0;
close = numel(s);
for i = open:numel(s)
    if s(i) == '('
        depth = depth + 1;
    elseif s(i) == ')'
        depth = depth - 1;
        if depth == 0
            close = i;
            return
        end
    end
end
end
