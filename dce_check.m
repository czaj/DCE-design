function info = dce_check(spec, varargin)
%DCE_CHECK Validate a specification and show how the toolbox reads it.
%
%   DCE_CHECK(spec) checks a specification (struct, .json file or Ngene .ngs /
%   .ngd file) without searching, and prints:
%     - the design columns with their levels and allowed level frequencies,
%     - every coefficient of every model with the level it belongs to (dummy
%       and effects coding), its base level, its prior and, for random
%       coefficients, its distribution,
%     - the constraints, and the share of random choice tasks they allow,
%     - notes about what the toolbox assumed,
%     - advice: the size of the problem, rows and blocks against level balance,
%       priors, constraints (see dce_advise).
%   Any problem raises an error (identifier dce:spec:*) that names the field.
%   info = DCE_CHECK(spec) returns the same as a struct and prints nothing;
%   DCE_CHECK(spec, 'Verbose', true) does both.
%
%   Run it before dce_design: it takes a moment and shows whether priors are
%   attached to the intended levels.
%
%   See also dce_design, dce_evaluate, dce_ngene.

verbose = nargout == 0;
for i = 1:2:numel(varargin)
    if strcmpi(char(varargin{i}), 'Verbose')
        verbose = logical(varargin{i + 1});
    else
        error('dce:options:unknown', 'Unknown option ''%s''.', char(varargin{i}));
    end
end
spec = dd_readspec(spec);
% full factorial and orthogonal designs set their own number of rows
ownRows = isfield(spec, 'type') && ~isempty(spec.type) && ~strcmpi(char(spec.type), 'efficient');
if ownRows
    if ~any(strcmpi(spec.type, {'full','orthogonal','orthogonal-simultaneous','orthogonal-sequential'}))
        dd_err('spec:type', 'Unknown design type. Use efficient, full, orthogonal or orthogonal-sequential.');
    end
    specC = spec;
    if ~isfield(specC, 'rows') || isempty(specC.rows)
        specC.rows = 1e6;
    end
    specC.balance = 'none';
    for f = {'frequency', 'blocks'}
        if isfield(specC, f{1}), specC = rmfield(specC, f{1}); end
    end
    C = dd_compile(specC);
    C.notes{end + 1} = sprintf(['Design type ''%s'': the number of rows and the level frequencies are set by the construction ' ...
        '(see dce_full, dce_orthogonal); the counts below are not binding.'], char(spec.type));
else
    C = dd_compile(spec);
end
dd_options(C.spec);

info = struct();
info.name = C.name;
info.rows = C.S;
info.blocks = C.B;
info.alternatives = C.altNames;
info.columns = table(C.colName(:), cellfun(@(v) mat2str(v), C.levels(:), 'UniformOutput', false), ...
    cellfun(@(v) mat2str(v'), C.lo(:), 'UniformOutput', false), cellfun(@(v) mat2str(min(v, C.S)'), C.hi(:), 'UniformOutput', false), ...
    'VariableNames', {'Column', 'Levels', 'MinCount', 'MaxCount'});
fx = {};
fa = fieldnames(C.fixed);
for i = 1:numel(fa)
    fb = fieldnames(C.fixed.(fa{i}));
    for q = 1:numel(fb)
        fx{end + 1} = sprintf('%s.%s = %g', fa{i}, fb{q}, C.fixed.(fa{i}).(fb{q})); %#ok<AGROW>
    end
end
info.fixed = fx(:);
info.models = cell(1, C.M);
distName = {'fixed', 'normal', 'lognormal', 'uniform', 'triangular'};
for m = 1:C.M
    Mm = C.model(m);
    n = Mm.P;
    par = Mm.parNames(:);
    coding = repmat({'linear'}, n, 1);
    attr = repmat({''}, n, 1);
    level = nan(n, 1);
    base = nan(n, 1);
    for ri = 1:numel(Mm.reg)
        r = Mm.reg(ri);
        coding(r.idx) = {r.coding};
        if ~strcmp(r.coding, 'linear')
            attr(r.idx) = {r.attr};
            level(r.idx) = r.nb;
            base(r.idx) = C.base.(r.attr);
        end
    end
    prior = cell(n, 1);
    random = cell(n, 1);
    for p = 1:n
        switch Mm.bdist(p)
            case 0, prior{p} = sprintf('%g', Mm.b0(p));
            case 1, prior{p} = sprintf('N(%g, %g)', Mm.ba(p), Mm.bb(p));
            case 2, prior{p} = sprintf('U(%g, %g)', Mm.ba(p), Mm.bb(p));
        end
        if Mm.rdist(p) == 0
            random{p} = '';
        else
            sg = '';
            if Mm.rsign(p) < 0, sg = '-'; end
            switch Mm.sdist(p)
                case 0, sp = sprintf('%g', Mm.s0(p));
                case 1, sp = sprintf('N(%g, %g)', Mm.sa(p), Mm.sb(p));
                case 2, sp = sprintf('U(%g, %g)', Mm.sa(p), Mm.sb(p));
            end
            random{p} = sprintf('%s%s(%s)', sg, distName{Mm.rdist(p) + 1}, sp);
        end
    end
    inCrit = Mm.eff(1:n);
    info.models{m} = struct('name', Mm.name, 'type', Mm.type, 'stat', Mm.stat, 'weight', C.w(m), ...
        'coefficients', table(par, coding, attr, level, base, prior, random, inCrit, ...
        'VariableNames', {'Coefficient', 'Coding', 'Attribute', 'Level', 'BaseLevel', 'Prior', 'Random', 'InCriterion'}));
end
info.reject = C.rejectExpr(:);
info.require = C.requireExpr(:);
info.dominancePairs = arrayfun(@(d) sprintf('%s-%s', C.altNames{d.a}, C.altNames{d.b}), C.dom, 'UniformOutput', false);
rs = RandStream('mt19937ar', 'Seed', 4711);
Dr = zeros(4000, C.K);
for k = 1:C.K
    Dr(:, k) = randi(rs, C.L(k), 4000, 1);
end
info.feasibleShare = mean(dd_rowok(C, Dr));
info.notes = C.notes(:);
if ownRows
    info.advice = struct('topic', {}, 'level', {}, 'text', {});
else
    info.advice = dd_advice(C);
end
info.specHash = C.hash;

if verbose
    fprintf('\nSpecification "%s": %d rows, %d block(s), alternatives %s\n', C.name, C.S, C.B, strjoin(C.altNames, ', '));
    fprintf('\nDesign columns (%d):\n', C.K);
    disp(info.columns);
    if ~isempty(fx)
        fprintf('Fixed levels: %s\n\n', strjoin(fx, '; '));
    end
    for m = 1:C.M
        mm = info.models{m};
        fprintf('Model %s (%s, weight %g, prior draws summarised by: %s):\n', mm.name, mm.type, mm.weight, mm.stat);
        disp(mm.coefficients);
    end
    if ~isempty(info.reject), fprintf('Rejected when true:\n'); fprintf('  %s\n', info.reject{:}); end
    if ~isempty(info.require), fprintf('Required:\n'); fprintf('  %s\n', info.require{:}); end
    fprintf('Random choice tasks that satisfy all constraints: %.1f%%\n', 100 * info.feasibleShare);
    for i = 1:numel(info.notes)
        fprintf('Note: %s\n', info.notes{i});
    end
    for i = 1:numel(info.advice)
        fprintf('Advice (%s, %s): %s\n', info.advice(i).topic, info.advice(i).level, info.advice(i).text);
    end
    fprintf('\n');
end
if nargout == 0
    clear info
end
end
