function result = dce_full(spec, varargin)
%DCE_FULL Full factorial design, or a random fraction of it.
%
%   result = DCE_FULL(spec) returns every possible choice task of the
%   specification: all combinations of the levels of all design columns (the
%   attributes of every alternative, and the scenario variables). spec is the
%   specification of dce_design; spec.rows is not needed. Tasks that break
%   the constraints of spec (reject, require, dominance) and repeats of a task
%   with interchangeable alternatives exchanged are left out, so the result is
%   the set of all distinct feasible tasks. The design is evaluated like any
%   other (D-error under the priors of spec, level counts, correlations) and
%   blocked if spec.blocks is given.
%
%   result = DCE_FULL(spec, 'Rows', n) draws n of these tasks at random
%   (without replacement), the same ones for a given 'Seed'.
%
%   Options (name-value):
%     'Rows'         []     size of a random fraction; [] = all tasks
%     'Constraints'  true   leave out infeasible and repeated tasks; false
%                           returns the raw full factorial
%     'MaxRows'      20000  refuse a larger design (use 'Rows')
%     and the evaluation options of dce_design ('Seed', 'EvalBayesDraws', ...).
%
%   A full factorial is the reference for what a fraction gives up, and a
%   source of candidate tasks. As a design it is only practical when it is
%   small or split into many blocks.
%
%   See also dce_design, dce_orthogonal, dce_evaluate.

spec = dd_readspec(spec);
nRows = [];
useCons = true;
maxRows = 20000;
keep = true(1, numel(varargin));
for i = 1:2:numel(varargin)
    switch lower(char(varargin{i}))
        case 'rows', nRows = varargin{i + 1}; keep([i i + 1]) = false;
        case 'constraints', useCons = logical(varargin{i + 1}); keep([i i + 1]) = false;
        case 'maxrows', maxRows = varargin{i + 1}; keep([i i + 1]) = false;
    end
end
varargin = varargin(keep);

% compile once to learn the columns (rows and blocks are set afterwards)
s0 = baseSpec(spec, 1e6);
C0 = dd_compile(s0);
o = dd_options(C0.spec, varargin{:});
nAll = prod(C0.L);
if nAll > 5e6
    dd_err('full:size', 'The full factorial has about 10^%.0f choice tasks; too many to list. Use dce_design for an efficient fraction.', log10(nAll));
end
g = cell(1, C0.K);
for k = 1:C0.K
    g{k} = 1:C0.L(k);
end
[g{:}] = ndgrid(g{:});
D = cell2mat(cellfun(@(v) v(:), g, 'UniformOutput', false));
nFeas = nAll;
if useCons
    D = D(dd_rowok(C0, D), :);
    nFeas = size(D, 1);
    D = D(~dd_dups(C0, D), :);
end
nDistinct = size(D, 1);
if ~isempty(nRows)
    if ~isnumeric(nRows) || ~isscalar(nRows) || nRows < 2 || nRows ~= round(nRows)
        dd_err('full:rows', 'Rows must be an integer of at least 2.');
    end
    if nRows > nDistinct
        dd_err('full:rows', 'Rows = %d, but only %d distinct feasible choice tasks exist.', nRows, nDistinct);
    end
    rs = RandStream('mt19937ar', 'Seed', mod(o.Seed + 31, 2^32 - 1));
    D = D(sort(randperm(rs, nDistinct, nRows)), :);
end
S = size(D, 1);
if S > maxRows
    dd_err('full:size', ['The full factorial has %d choice tasks (limit MaxRows = %d). Give ''Rows'' for a random fraction, ' ...
        'raise ''MaxRows'', or use dce_design for an efficient fraction.'], S, maxRows);
end
if S < 2
    dd_err('full:size', 'Fewer than two feasible choice tasks exist.');
end

% evaluate as an existing design
s1 = baseSpec(spec, S);
if isfield(spec, 'blocks') && ~isempty(spec.blocks) && spec.blocks > 1
    s1.blocks = spec.blocks;
end
if ~useCons
    s1.unique = false;
end
C = dd_compile(s1);
y = [];
bout = [];
if C.B > 1
    V = zeros(S, C.K);
    for k = 1:C.K
        V(:, k) = C.LV{k}(D(:, k));
    end
    [y, b] = dce_block(V, C.B, 'Names', C.colName, 'Verbose', false, 'Seed', o.Seed, 'Restarts', min(o.BlockRestarts, 10), 'Iterations', max(1000, 5 * S));
    [y, ord] = sort(y);
    D = D(ord, :);
    bout = struct('Psi', b.Psi, 'PsiBound', b.PsiBound, 'absentLevels', height(b.absent), 'absentAvoidable', nnz(b.absent.Avoidable), ...
        'blockSizes', b.blockSizes, 'extremeColumns', {{}}, 'extremesOK', true, 'informationShare', [], 'leaveOneBlockOut', []);
end
meas = dd_measures(C, D, y, o);
result = dd_result(C, D, y, o, meas);
result.blocking = bout;
if isempty(nRows)
    txt = sprintf('full factorial: %d possible choice tasks', nAll);
    if useCons
        txt = sprintf('%s, %d feasible, %d distinct (all used)', txt, nFeas, nDistinct);
    end
else
    txt = sprintf('random fraction of the full factorial: %d of %d distinct feasible choice tasks (%d possible)', S, nDistinct, nAll);
end
result.search = struct('construction', txt, 'possibleTasks', nAll, 'feasibleTasks', nFeas, 'distinctTasks', nDistinct, 'seed', o.Seed);
result.advice = dd_advice(C, result);
if o.Verbose >= 1
    dce_report(result);
end
end

function s = baseSpec(spec, rows)
s = spec;
s.rows = rows;
s.balance = 'none';
for f = {'frequency', 'blocks', 'type', 'start'}
    if isfield(s, f{1}), s = rmfield(s, f{1}); end
end
end
