function result = dce_orthogonal(spec, varargin)
%DCE_ORTHOGONAL Orthogonal design for a choice experiment.
%
%   result = DCE_ORTHOGONAL(spec) builds a design in which the attribute
%   columns are uncorrelated and every level appears equally often, without
%   using the priors. spec is the specification of dce_design.
%
%   'Type', 'simultaneous' (default): one orthogonal array holds all design
%   columns, so attributes are orthogonal within and between alternatives.
%   'Type', 'sequential': an orthogonal array is built for the attributes of
%   one alternative (with the scenario variables); the other alternatives use
%   the same array with the rows in another order. Attributes are orthogonal
%   within each alternative, not between alternatives. This needs far fewer
%   rows. The row orders are chosen among 'Orders' random ones: first the one
%   with the fewest tasks that break the constraints of spec, then the lowest
%   D-error under the priors of spec.
%
%   Rows. A strength-2 orthogonal array needs a number of rows that is a
%   multiple of Li*Lj for every pair of columns with Li and Lj levels. If
%   spec.rows is not such a number, or no array is found with it, the next
%   suitable number is used and the result says so. If none is found up to
%   'MaxRows', the most nearly orthogonal array found is returned, with the
%   largest correlation in the notes.
%
%   Blocks. With spec.blocks, the block is first tried as one more column of
%   the array, which makes the blocks orthogonal to the attributes (every
%   level equally often in every block); otherwise dce_block is used.
%
%   Options (name-value):
%     'Type'     'simultaneous' or 'sequential'
%     'Tries'    20     attempts to find an orthogonal array for a given size
%     'Orders'   200    row orders tried for the other alternatives (sequential)
%     'MaxRows'  4 times the first suitable number of rows
%     'Seed'     1
%     and the evaluation options of dce_design.
%
%   Orthogonal designs cannot respect constraints on choice tasks in general;
%   result.checks lists the tasks that break them (dominated alternatives are
%   common in unlabelled orthogonal designs). The design is evaluated like any
%   other, so its D-error can be compared with an efficient design under the
%   same priors.
%
%   The arrays are found with the J2 algorithm of Xu (2002), "An algorithm for
%   constructing orthogonal and nearly-orthogonal arrays with mixed levels and
%   small runs", Technometrics 44, 356-368.
%
%   See also dce_design, dce_full, dce_evaluate.

spec = dd_readspec(spec);
type = 'simultaneous';
tries = 20;
orders = 200;
maxRows = [];
keep = true(1, numel(varargin));
for i = 1:2:numel(varargin)
    switch lower(char(varargin{i}))
        case 'type', type = lower(char(varargin{i + 1})); keep([i i + 1]) = false;
        case 'tries', tries = varargin{i + 1}; keep([i i + 1]) = false;
        case 'orders', orders = varargin{i + 1}; keep([i i + 1]) = false;
        case 'maxrows', maxRows = varargin{i + 1}; keep([i i + 1]) = false;
    end
end
varargin = varargin(keep);
if any(strcmp(type, {'sim', 'simultaneous'}))
    type = 'simultaneous';
elseif any(strcmp(type, {'seq', 'sequential'}))
    type = 'sequential';
else
    dd_err('orthogonal:type', 'Type must be ''simultaneous'' or ''sequential''.');
end
rowsWanted = 2;
if isfield(spec, 'rows') && ~isempty(spec.rows)
    rowsWanted = spec.rows;
end
C0 = dd_compile(baseSpec(spec, 1e6));
o = dd_options(C0.spec, varargin{:});
rs = RandStream('mt19937ar', 'Seed', mod(o.Seed + 53, 2^32 - 1));
B = 1;
if isfield(spec, 'blocks') && ~isempty(spec.blocks)
    B = spec.blocks;
end
K = C0.K;
Pmax = max([C0.model.Pt]);
minRows = max(rowsWanted, ceil(Pmax / (C0.J - 1)));

% ---- which columns form the array --------------------------------------------
if strcmp(type, 'simultaneous')
    groups = {1:K};
    arrayCols = 1:K;
else
    % alternatives with attributes; the first one carries the scenario variables
    altsWith = unique(C0.colAlt(C0.colAlt > 0), 'stable');
    if isempty(altsWith)
        dd_err('orthogonal:type', 'A sequential orthogonal design needs alternatives with attributes.');
    end
    groups = cell(1, numel(altsWith));
    for a = 1:numel(altsWith)
        groups{a} = find(C0.colAlt == altsWith(a));
    end
    groups{1} = [find(C0.colAlt == 0) groups{1}];
    arrayCols = groups{1};
    for a = 2:numel(groups)
        if ~isequal(C0.L(groups{a}), C0.L(groups{1}(C0.colAlt(groups{1}) > 0)))
            dd_err('orthogonal:type', ['A sequential orthogonal design needs the same numbers of levels in every alternative ' ...
                '(alternative %s differs). Use ''Type'', ''simultaneous''.'], C0.altNames{altsWith(a)});
        end
    end
end
Lv = C0.L(arrayCols);

% ---- number of rows ---------------------------------------------------------------
step = 1;
for a = 1:numel(Lv)
    step = lcm(step, Lv(a));
    for b = a + 1:numel(Lv)
        step = lcm(step, Lv(a) * Lv(b));
    end
end
S0 = ceil(minRows / step) * step;
if isempty(maxRows)
    maxRows = max(4 * step, S0);
end
found = false;
bestA = [];
bestDev = Inf;
blockInArray = false;
S = S0;
while S <= maxRows
    % with blocks: first try the block as one more column of the array
    if B > 1 && mod(S, B) == 0 && all(mod(S, B * Lv) == 0)
        [A, isOA] = dd_oa(S, [Lv B], tries, rs);
        if isOA
            found = true;
            blockInArray = true;
            break
        end
    end
    [A, isOA, dev] = dd_oa(S, Lv, tries, rs);
    if isOA
        found = true;
        break
    end
    if dev < bestDev
        bestDev = dev;
        bestA = A;
    end
    S = S + step;
end
if ~found
    A = bestA;
    S = size(A, 1);
end
y = [];
if blockInArray
    y = A(:, end);
    A = A(:, 1:end - 1);
end

% ---- assemble the design -------------------------------------------------------------
D = zeros(S, K);
D(:, arrayCols) = A;
C = dd_compile(evalSpec(spec, S, B));
nOrders = 0;
if strcmp(type, 'sequential')
    own = groups{1}(C0.colAlt(groups{1}) > 0);
    Mm = C.model(1);
    bestScore = [Inf Inf];
    Dbest = D;
    for t = 1:orders
        Dt = D;
        for a = 2:numel(groups)
            Dt(:, groups{a}) = A(randperm(rs, S), ismember(arrayCols, own));
        end
        bad = nnz(~dd_rowok(C, Dt)) + nnz(dd_dups(C, Dt));
        It = sum(dd_rowinfo(Mm, dd_X(Mm, Dt), struct('R', 1, 'B', Mm.b0)), 3);
        Mfix = Mm;
        Mfix.stat = 'fixed';
        if strcmp(Mm.type, 'mnl')
            f = dd_crit(It, Mfix, 'd');
        else
            f = 0;
        end
        if bad < bestScore(1) || (bad == bestScore(1) && f < bestScore(2))
            bestScore = [bad f];
            Dbest = Dt;
        end
    end
    D = Dbest;
    nOrders = orders;
end

% ---- blocks, evaluation ---------------------------------------------------------------
bout = [];
if B > 1
    V = zeros(S, K);
    for k = 1:K
        V(:, k) = C.LV{k}(D(:, k));
    end
    if isempty(y)
        y = dce_block(V, B, 'Names', C.colName, 'Verbose', false, 'Seed', o.Seed, 'Restarts', o.BlockRestarts, 'Iterations', max(1000, 20 * S));
    end
    [y, ord] = sort(y);
    D = D(ord, :);
    V = V(ord, :);
    [~, b] = dce_block(V, B, 'Evaluate', y, 'Names', C.colName, 'Verbose', false);
    bout = struct('Psi', b.Psi, 'PsiBound', b.PsiBound, 'absentLevels', height(b.absent), 'absentAvoidable', nnz(b.absent.Avoidable), ...
        'blockSizes', b.blockSizes, 'extremeColumns', {{}}, 'extremesOK', true, 'informationShare', [], 'leaveOneBlockOut', []);
end
meas = dd_measures(C, D, y, o);
result = dd_result(C, D, y, o, meas);
spec.type = 'orthogonal';
if strcmp(type, 'sequential')
    spec.type = 'orthogonal-sequential';
end
result.spec = spec;
result.specHash = dd_spechash(spec);
result.blocking = bout;
cr = meas.correlations;
within = 0;
for g = 1:numel(groups)
    if strcmp(type, 'simultaneous')
        cols = 1:K;
    else
        cols = groups{g};
    end
    sub = abs(cr(cols, cols));
    sub(1:numel(cols) + 1:end) = 0;
    within = max([within; sub(:)]);
end
txt = sprintf('%s orthogonal design, %d rows', type, S);
if found
    txt = sprintf('%s; orthogonal array of strength 2 for %d columns (levels %s)', txt, numel(Lv), strtrim(sprintf('%d ', Lv)));
else
    txt = sprintf('%s; NOT orthogonal: no orthogonal array found up to %d rows, nearest array used', txt, maxRows);
end
result.search = struct('construction', txt, 'type', type, 'orthogonal', found, 'rows', S, 'rowsRequested', rowsWanted, ...
    'maxCorrelation', within, 'blockInArray', blockInArray, 'ordersTried', nOrders, 'seed', o.Seed);
if S ~= rowsWanted && isfield(spec, 'rows') && ~isempty(spec.rows)
    result.notes{end + 1} = sprintf(['%d rows were asked for; an orthogonal array for these columns needs a multiple of %d rows, ' ...
        'and %d is the smallest for which one was found.'], rowsWanted, step, S);
end
if ~found
    result.notes{end + 1} = sprintf(['No orthogonal array was found for %d to %d rows; the design is nearly orthogonal ' ...
        '(largest correlation between its columns %.3f). More ''Tries'' or a larger ''MaxRows'' may find one.'], S0, maxRows, within);
end
if strcmp(type, 'sequential')
    result.notes{end + 1} = sprintf(['Sequential: attributes are orthogonal within each alternative (largest correlation %.3f), ' ...
        'not between alternatives. The row order of the other alternatives is the best of %d random orders.'], within, nOrders);
end
if blockInArray
    result.notes{end + 1} = 'The blocks are a column of the orthogonal array: every level appears equally often in every block.';
end
if ~result.checks.constraintsOK
    result.notes{end + 1} = sprintf(['%d of %d tasks break the constraints of the specification (see checks.violations); ' ...
        'an orthogonal design cannot avoid them. dce_design can.'], size(result.checks.violations, 1), S);
end
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

function s = evalSpec(spec, rows, B)
s = baseSpec(spec, rows);
if B > 1
    s.blocks = B;
end
end
