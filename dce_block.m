function [y, out] = dce_block(X, B, varargin)
%DCE_BLOCK Balanced blocking of a choice-experiment design.
%
%   [y, out] = DCE_BLOCK(X, B) assigns the S rows (choice tasks) of design X
%   to B blocks of (nearly) equal size so that each block shows the levels of
%   every design column as evenly as possible. X is S-by-K in wide format: one
%   column per attribute per alternative, plus any scenario variables. Entries
%   can be design codes, actual levels or (in a table) categorical values.
%   y is S-by-1 with block numbers 1..B (block 1 contains row 1); out holds
%   diagnostics.
%
%   Criterion (Bliemer, Hess, Collins & Amaris 2026, "On the blocking of
%   non-orthogonal experimental designs for choice experiments", ITLS-WP-26-22):
%   the sum over blocks and column groups of the squared coefficient of
%   variation of relative level frequencies, h^2. Balance is reported as
%   Psi = (1 - h/hmax)*100, where 100 means that every block is a miniature
%   replicate of a level-balanced design. This implementation adds:
%     1. Absent levels have lexicographic priority. The search first minimises
%        the number of (block, group, level) cells with zero count among cells
%        that the block is large enough to fill, and only then h^2.
%     2. Lowest and highest levels of ordered numeric attributes (cost, time)
%        can be required in every block ('Extremes').
%     3. Optional MNL information diagnostics and information-based choice
%        among near-best blockings ('Info').
%     4. Swaps that leave the objective unchanged are accepted to cross
%        plateaus, and the best restart solutions are polished by a full
%        pairwise sweep.
%
%   Name-value options
%     'BlockSizes'  1-by-B block sizes (default: as equal as possible).
%     'Groups'      Column groups whose counts are pooled (across-alternative
%                   balance for generic attributes): a cell array of column
%                   index vectors, or a 1-by-K vector of group ids. Columns of
%                   a group should share the same level set. Default: every
%                   column on its own (within-alternative balance). Pooling can
%                   certify perfect balance while one alternative never shows a
%                   level in a block, so out.PsiWithin is always reported too.
%     'Names'       1-by-K column names (default: table names or X1..XK).
%     'RequireAll'  true (default): minimise absent levels first.
%     'Extremes'    Indices of numeric columns whose lowest and highest level
%                   must appear in every block (lexicographic priority).
%     'Restarts'    Random starts (default 100 for B <= 4, 20 for B <= 20,
%                   5 otherwise: with many blocks, longer runs beat more runs).
%     'Iterations'  Random swap proposals per start (default keeps a total
%                   budget of 100*max(2000,10*S) proposals).
%     'AcceptEqual' Accept swaps that leave the objective unchanged (default
%                   true); the best solution found is kept.
%     'Polish'      Sweep all pairs of tasks until no swap improves (default
%                   true when S <= 300). Applied to all restart solutions when
%                   S <= 60, otherwise to the 10 best.
%     'Info'        Struct for MNL information diagnostics:
%                     .X      S-by-J-by-P coded attributes of each task and
%                             alternative, exactly as they enter utility
%                             (dummies, interactions, ASCs),
%                     .beta   P-by-1 priors,
%                     .avail  optional S-by-J logical availability.
%     'Select'      'balance' (default without Info) or 'info' (default with
%                   Info). With 'info', among distinct restart solutions with
%                   the fewest absent levels and Psi within 'Tolerance' of the
%                   best, pick the one whose worst leave-one-block-out D-error
%                   is smallest, i.e. the design that loses least if one block
%                   is under-filled in fieldwork.
%     'Tolerance'   Percentage points of Psi that may be traded for
%                   information with 'Select','info' (default 1).
%     'Evaluate'    S-by-1 existing blocking (e.g. from Ngene) to evaluate
%                   without searching. Labels are kept as given.
%     'Seed'        Seed of a private random stream (default 1); the global
%                   random number generator is not touched.
%     'Verbose'     Print a summary (default true).
%
%   Main fields of out
%     Psi, h, hmax    balance under 'Groups'
%     PsiWithin       balance with every column on its own
%     PsiBound        upper bound on attainable Psi from per-group relaxations
%                     (integer counts and fixed design totals). Rows tie the
%                     columns together, so the attainable maximum can be lower.
%     absentRequired  number of absent cells the search tried to avoid
%     absent          table of (block, column, level) never shown in a block
%     counts          level counts by block and level index (one column per
%                     design column; levels are listed in out.levels)
%     extremes        lowest/highest level present per block ('Extremes')
%     CramerV         Cramer's V of the block-by-level table of each group
%     corrNgene       correlation criterion of Ngene (eq. 5 of the paper),
%                     for comparison only; it depends on arbitrary labels
%     designBalanced  true if every column is level balanced overall; Psi is
%                     comparable only between blockings of the same design
%     info            MNL diagnostics (with 'Info'): Derror of the full
%                     design; share(b) = trace(I\I_b)/P, block b's share of
%                     the information (equal blocks: 1/B); DerrorLOBO(b),
%                     the D-error without block b, and its ratio to the full
%                     D-error; rank and D-error of each block on its own (Inf
%                     when S_b*(J-1) < P, which is common, hence the shares)
%     candidates      distinct restart solutions and their criteria
%     search          settings, seed and timing
%
%   See test_dce_block.m for a worked example (the paper's hotel design) and
%   README.md.
%
%   Version 1.0, 2026-11-14.

t0 = tic;
p = inputParser;
p.FunctionName = 'dce_block';
addRequired(p, 'X');
addRequired(p, 'B', @(v) isnumeric(v) && isscalar(v) && v >= 2 && v == round(v));
addParameter(p, 'BlockSizes', []);
addParameter(p, 'Groups', []);
addParameter(p, 'Names', {});
addParameter(p, 'RequireAll', true);
addParameter(p, 'Extremes', []);
addParameter(p, 'Restarts', []);
addParameter(p, 'Iterations', []);
addParameter(p, 'AcceptEqual', true);
addParameter(p, 'Polish', []);
addParameter(p, 'Info', []);
addParameter(p, 'Select', '');
addParameter(p, 'Tolerance', 1);
addParameter(p, 'Evaluate', []);
addParameter(p, 'Seed', 1);
addParameter(p, 'Verbose', true);
parse(p, X, B, varargin{:});
o = p.Results;

% ---- design values -------------------------------------------------------
[V, names, isNum] = toValues(X, o.Names);
[S, K] = size(V);
if any(isnan(V(:)))
    error('dce_block:nan', ['X contains NaN. Use wide format (one column per ' ...
        'attribute per alternative) and drop columns that do not apply.']);
end
if B > S
    error('dce_block:B', 'B (%d) exceeds the number of tasks (%d).', B, S);
end

% ---- block sizes ---------------------------------------------------------
evaluateOnly = ~isempty(o.Evaluate);
if evaluateOnly
    y0 = o.Evaluate(:);
    if numel(y0) ~= S || any(y0 < 1 | y0 ~= round(y0)) || max(y0) ~= B
        error('dce_block:evaluate', 'Evaluate must be S-by-1 with labels 1..B.');
    end
    Sb = accumarray(y0, 1, [B 1])';
    if any(Sb == 0)
        error('dce_block:evaluate', 'Some block in Evaluate is empty.');
    end
elseif isempty(o.BlockSizes)
    Sb = floor(S / B) * ones(1, B);
    r = S - sum(Sb);
    Sb(1:r) = Sb(1:r) + 1;
else
    Sb = o.BlockSizes(:)';
    if numel(Sb) ~= B || sum(Sb) ~= S || any(Sb < 1 | Sb ~= round(Sb))
        error('dce_block:sizes', 'BlockSizes must be B positive integers summing to S.');
    end
end

% ---- groupings -----------------------------------------------------------
gid = parseGroups(o.Groups, K);
G = makeGrouping(V, gid);
Gw = makeGrouping(V, 1:K);
if G.M == 0
    error('dce_block:constant', 'All design columns are constant.');
end
grouped = any(G.gsize > 1);

% ---- cells that must not be empty ---------------------------------------
ext = unique(o.Extremes(:))';
if any(ext < 1 | ext > K | ext ~= round(ext))
    error('dce_block:extremes', 'Extremes must be column indices.');
end
if any(~isNum(ext))
    error('dce_block:extremes', 'Extremes columns must be numeric (ordered levels).');
end
extGroup = false(1, G.M);
extGroup(G.colGroup(ext(G.colGroup(ext) > 0))) = true;
req = false(B, G.M, G.Lmax);
for m = 1:G.M
    for b = 1:B
        n = Sb(b) * G.gsize(m);
        if o.RequireAll && n >= G.L(m)
            req(b, m, 1:G.L(m)) = true;
        end
        if extGroup(m) && n >= 2
            req(b, m, [1 G.L(m)]) = true;
        end
    end
end
w = G.L ./ (Sb(:).^2 * G.gsize.^2);          % B-by-M weights
Wcell = repmat(w, [1 1 G.Lmax]);
tol = 1e-9;

% ---- MNL information (optional) -------------------------------------------
haveInfo = ~isempty(o.Info);
if haveInfo
    [Is, Ifull, P] = mnlInfo(o.Info, S);
    if logdetPD(Ifull) == -Inf
        warning('dce_block:info', ['The full design does not identify all ' ...
            'parameters at these priors; information diagnostics skipped.']);
        haveInfo = false;
    end
end
sel = lower(char(o.Select));
if isempty(sel)
    if haveInfo, sel = 'info'; else, sel = 'balance'; end
end
if strcmp(sel, 'info') && ~haveInfo
    sel = 'balance';
end

% ---- search ----------------------------------------------------------------
if isempty(o.Restarts)
    if B <= 4, R = 100; elseif B <= 20, R = 20; else, R = 5; end
else
    R = o.Restarts;
end
if isempty(o.Iterations), I = round(100 * max(2000, 10 * S) / R); else, I = o.Iterations; end
if isempty(o.Polish), doPolish = S <= 300; else, doPolish = logical(o.Polish); end
cand = table();
restartPsi = [];

if evaluateOnly
    y = y0;
else
    rs = RandStream('mt19937ar', 'Seed', o.Seed);
    labels = repelem(1:B, Sb)';
    finals = zeros(S, R);
    fmiss = zeros(1, R);
    fh2 = zeros(1, R);
    dup = grouped;
    for r = 1:R
        y = zeros(S, 1);
        y(randperm(rs, S)) = labels;
        C = countCells(G, y, B);
        h2 = sum(Wcell(:) .* C(:).^2) - B * G.M;
        miss = nnz(req & C == 0);
        by = y; bh2 = h2; bmiss = miss;
        for it = 1:I
            s1 = randi(rs, S);
            b1 = y(s1);
            s2 = randi(rs, S);
            while y(s2) == b1
                s2 = randi(rs, S);
            end
            b2 = y(s2);
            [dm, dh, u, du] = evalSwap(s1, s2, b1, b2, G, C, Wcell, req, B, dup);
            if isempty(u)
                continue
            end
            if dm < 0 || (dm == 0 && dh < -tol) || (o.AcceptEqual && dm == 0 && abs(dh) <= tol)
                C(u) = C(u) + du;
                y(s1) = b2;
                y(s2) = b1;
                h2 = h2 + dh;
                miss = miss + dm;
                if miss < bmiss || (miss == bmiss && h2 < bh2 - tol)
                    by = y; bh2 = h2; bmiss = miss;
                end
            end
        end
        finals(:, r) = by;
        fmiss(r) = bmiss;
        fh2(r) = bh2;
    end

    if doPolish
        if S <= 60
            toPolish = 1:R;
        else
            [~, ord] = sortrows([fmiss(:) fh2(:)]);
            toPolish = ord(1:min(10, R))';
        end
        for r = toPolish
            finals(:, r) = polishSwaps(finals(:, r), G, B, Wcell, req, tol);
        end
    end

    for r = 1:R
        C = countCells(G, finals(:, r), B);
        fh2(r) = max(sum(Wcell(:) .* C(:).^2) - B * G.M, 0);
        fmiss(r) = nnz(req & C == 0);
        finals(:, r) = canonLabels(finals(:, r), Sb);
    end
    hmaxG = sqrt(B * sum(G.L - 1));
    restartPsi = 100 * (1 - sqrt(fh2) / hmaxG);

    [~, ia, ic] = unique(finals', 'rows', 'stable');
    nc = numel(ia);
    cMiss = fmiss(ia)';
    cPsi = restartPsi(ia)';
    cCount = accumarray(ic(:), 1, [nc 1]);
    cWorst = nan(nc, 1);
    cShareDev = nan(nc, 1);
    cPsiW = zeros(nc, 1);
    for j = 1:nc
        stj = balanceStats(Gw, finals(:, ia(j)), Sb, B);
        cPsiW(j) = stj.Psi;
    end
    bestMiss = min(cMiss);
    inBest = cMiss == bestMiss;
    bestPsi = max(cPsi(inBest));
    if strcmp(sel, 'info')
        pool = find(inBest & cPsi >= bestPsi - o.Tolerance - 1e-9);
    else
        pool = find(inBest & cPsi >= bestPsi - 1e-9);
    end
    if haveInfo
        for j = 1:nc
            st = infoStats(Is, Ifull, finals(:, ia(j)), B, P);
            cWorst(j) = max(st.LOBOratio);
            cShareDev(j) = max(abs(st.share - Sb / S));
        end
    end
    % ties are broken by within-alternative balance (matters for pooled groups)
    if strcmp(sel, 'info')
        [~, k] = sortrows([cWorst(pool) -cPsi(pool) -cPsiW(pool)]);
    else
        [~, k] = sortrows(-cPsiW(pool));
    end
    pick = pool(k(1));
    y = finals(:, ia(pick));
    cand = table((1:nc)', cPsi, cPsiW, cMiss, cWorst, cShareDev, cCount, (1:nc)' == pick, ...
        'VariableNames', {'Solution', 'Psi', 'PsiWithin', 'AbsentRequired', 'WorstLOBOratio', ...
        'MaxShareDeviation', 'Restarts', 'Chosen'});
    cand = sortrows(cand, {'AbsentRequired', 'Psi'}, {'ascend', 'descend'});
end

% ---- diagnostics -----------------------------------------------------------
stG = balanceStats(G, y, Sb, B);
stW = balanceStats(Gw, y, Sb, B);

out = struct();
out.y = y;
out.blockSizes = Sb;
out.Psi = stG.Psi;
out.h = stG.h;
out.hmax = stG.hmax;
out.PsiWithin = stW.Psi;
out.PsiBound = psiBound(G, Sb, B);
out.PsiWithinBound = psiBound(Gw, Sb, B);
out.absentRequired = nnz(req & stG.C == 0);
out.grouping = struct('groups', {G.cols}, 'levels', {G.levels}, 'dropped', G.dropped);
out.names = names;
out.levels = cell(1, K);
for k = Gw.active
    out.levels{k} = Gw.levels{Gw.colGroup(k)};
end
out.CramerV = stG.V;

% counts by block and level (within-alternative)
Lw = Gw.Lmax;
cnt = nan(B * Lw, numel(Gw.active));
for j = 1:numel(Gw.active)
    k = Gw.active(j);
    m = Gw.colGroup(k);
    for b = 1:B
        cnt((b - 1) * Lw + (1:Gw.L(m)), j) = reshape(stW.C(b, m, 1:Gw.L(m)), [], 1);
    end
end
vn = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names(Gw.active)), {'Block', 'Level'});
out.counts = array2table([repelem((1:B)', Lw) repmat((1:Lw)', B, 1) cnt], ...
    'VariableNames', [{'Block', 'Level'}, vn]);

% absent levels (within-alternative)
ab = zeros(0, 4);
for k = Gw.active
    m = Gw.colGroup(k);
    for b = 1:B
        for l = 1:Gw.L(m)
            if stW.C(b, m, l) == 0
                ab(end + 1, :) = [b k Gw.levels{m}(l) Sb(b) >= Gw.L(m)]; %#ok<AGROW>
            end
        end
    end
end
colNames = names(ab(:, 2));
out.absent = table(ab(:, 1), colNames(:), ab(:, 3), logical(ab(:, 4)), ...
    'VariableNames', {'Block', 'Column', 'Level', 'Avoidable'});

% extremes
ex = zeros(0, 4);
for k = ext
    m = Gw.colGroup(k);
    if m == 0, continue; end
    for b = 1:B
        z = Gw.Z(y == b, k);
        ex(end + 1, :) = [b k any(z == 1) any(z == Gw.L(m))]; %#ok<AGROW>
    end
end
exNames = names(ex(:, 2));
out.extremes = table(ex(:, 1), exNames(:), logical(ex(:, 3)), logical(ex(:, 4)), ...
    'VariableNames', {'Block', 'Column', 'HasLowest', 'HasHighest'});

% Ngene's correlation criterion (eq. 5), for comparison
r2 = [];
for k = Gw.active
    c = corrcoef(Gw.Z(:, k), y);
    r2(end + 1) = c(1, 2)^2; %#ok<AGROW>
end
out.corrNgene = sqrt(mean(r2));

% overall level balance of the design
bal = true;
for k = Gw.active
    t = accumarray(Gw.Z(:, k), 1, [Gw.L(Gw.colGroup(k)) 1]);
    bal = bal && all(t == t(1));
end
out.designBalanced = bal;

if haveInfo
    out.info = infoStats(Is, Ifull, y, B, P);
else
    out.info = [];
end
out.candidates = cand;
out.search = struct('evaluateOnly', evaluateOnly, 'restarts', R, 'iterations', I, ...
    'polish', doPolish, 'acceptEqual', o.AcceptEqual, 'select', sel, ...
    'tolerance', o.Tolerance, 'seed', o.Seed, 'restartPsi', restartPsi, ...
    'groupWarnings', {G.warn}, 'seconds', toc(t0));

if o.Verbose
    printSummary(out, S, K, B, grouped, G);
end
end

% ===========================================================================
function [V, names, isNum] = toValues(X, names)
if istable(X)
    K = width(X);
    V = zeros(height(X), K);
    isNum = false(1, K);
    if isempty(names), names = X.Properties.VariableNames; end
    for k = 1:K
        x = X{:, k};
        if isnumeric(x) || islogical(x)
            V(:, k) = double(x);
            isNum(k) = true;
        else
            [~, ~, V(:, k)] = unique(x);
        end
    end
else
    V = double(X);
    K = size(V, 2);
    isNum = true(1, K);
    if isempty(names)
        names = arrayfun(@(k) sprintf('X%d', k), 1:K, 'UniformOutput', false);
    end
end
names = reshape(cellstr(names), 1, []);
if numel(names) ~= K
    error('dce_block:names', 'Names must have one entry per column.');
end
end

function gid = parseGroups(Gr, K)
if isempty(Gr)
    gid = 1:K;
    return
end
if iscell(Gr)
    gid = zeros(1, K);
    for m = 1:numel(Gr)
        c = Gr{m}(:)';
        if any(c < 1 | c > K | c ~= round(c))
            error('dce_block:groups', 'Groups must contain column indices.');
        end
        if any(gid(c) > 0)
            error('dce_block:groups', 'A column appears in more than one group.');
        end
        gid(c) = m;
    end
    un = find(gid == 0);
    gid(un) = numel(Gr) + (1:numel(un));
else
    if numel(Gr) ~= K
        error('dce_block:groups', 'A group-id vector must have one entry per column.');
    end
    [~, ~, gid] = unique(Gr(:)', 'stable');
    gid = gid(:)';
end
end

function G = makeGrouping(V, gid)
[S, K] = size(V);
G.Z = zeros(S, K);
G.colGroup = zeros(1, K);
G.cols = {};
G.L = [];
G.gsize = [];
G.levels = {};
G.dropped = [];
G.warn = {};
m = 0;
for g = unique(gid, 'stable')
    cols = find(gid == g);
    vals = V(:, cols);
    lev = unique(vals(:))';
    if numel(lev) < 2
        G.dropped = [G.dropped cols];
        continue
    end
    m = m + 1;
    for k = cols
        [~, G.Z(:, k)] = ismember(V(:, k), lev);
        nk = numel(unique(V(:, k)));
        if numel(cols) > 1 && nk < numel(lev)
            G.warn{end + 1} = sprintf('column %d uses %d of the %d levels of its group', ...
                k, nk, numel(lev));
        end
    end
    G.colGroup(cols) = m;
    G.cols{m} = cols;
    G.L(m) = numel(lev);
    G.gsize(m) = numel(cols);
    G.levels{m} = lev;
end
G.M = m;
G.Lmax = max([G.L 1]);
G.active = find(G.colGroup > 0);
end

function C = countCells(G, y, B)
n = B * G.M * G.Lmax;
c = zeros(n, 1);
for k = G.active
    idx = y(:) + B * (G.colGroup(k) - 1) + B * G.M * (G.Z(:, k) - 1);
    c = c + accumarray(idx, 1, [n 1]);
end
C = reshape(c, B, G.M, G.Lmax);
end

function [dmiss, dh2, u, du] = evalSwap(s1, s2, b1, b2, G, C, Wcell, req, B, dup)
% Change in (absent cells, h^2) if tasks s1 (block b1) and s2 (block b2) swap.
nBM = B * G.M;
idx = zeros(4 * numel(G.active), 1);
d = idx;
n = 0;
for k = G.active
    l1 = G.Z(s1, k);
    l2 = G.Z(s2, k);
    if l1 ~= l2
        mo = B * (G.colGroup(k) - 1);
        o1 = nBM * (l1 - 1) + mo;
        o2 = nBM * (l2 - 1) + mo;
        idx(n + 1:n + 4) = [b1 + o1; b1 + o2; b2 + o1; b2 + o2];
        d(n + 1:n + 4) = [-1; 1; 1; -1];
        n = n + 4;
    end
end
if n == 0
    dmiss = 0; dh2 = 0; u = []; du = [];
    return
end
idx = idx(1:n);
d = d(1:n);
if dup
    [u, ~, j] = unique(idx);
    du = accumarray(j, d);
    keep = du ~= 0;
    u = u(keep);
    du = du(keep);
    if isempty(u)
        dmiss = 0; dh2 = 0;
        return
    end
else
    u = idx;
    du = d;
end
c = C(u);
cn = c + du;
dh2 = sum(Wcell(u) .* (cn.^2 - c.^2));
rq = req(u);
dmiss = nnz(rq & cn == 0) - nnz(rq & c == 0);
end

function y = polishSwaps(y, G, B, Wcell, req, tol)
C = countCells(G, y, B);
S = numel(y);
dup = any(G.gsize > 1);
improved = true;
sweeps = 0;
while improved && sweeps < 200
    improved = false;
    sweeps = sweeps + 1;
    for s1 = 1:S - 1
        for s2 = s1 + 1:S
            b1 = y(s1);
            b2 = y(s2);
            if b1 == b2, continue; end
            [dm, dh, u, du] = evalSwap(s1, s2, b1, b2, G, C, Wcell, req, B, dup);
            if ~isempty(u) && (dm < 0 || (dm == 0 && dh < -tol))
                C(u) = C(u) + du;
                y(s1) = b2;
                y(s2) = b1;
                improved = true;
            end
        end
    end
end
end

function yc = canonLabels(y, Sb)
% Relabel blocks by first appearance, within classes of equal block size.
B = numel(Sb);
map = zeros(1, B);
used = false(1, B);
yc = zeros(size(y));
for s = 1:numel(y)
    b = y(s);
    if map(b) == 0
        c = find(Sb == Sb(b) & ~used, 1, 'first');
        map(b) = c;
        used(c) = true;
    end
    yc(s) = map(b);
end
end

function st = balanceStats(G, y, Sb, B)
C = countCells(G, y, B);
w = G.L ./ (Sb(:).^2 * G.gsize.^2);
Q = sum(C.^2, 3);
st.h2 = max(sum(w(:) .* Q(:)) - B * G.M, 0);
st.h = sqrt(st.h2);
st.hmax = sqrt(B * sum(G.L - 1));
st.Psi = 100 * (1 - st.h / st.hmax);
st.C = C;
st.V = zeros(1, G.M);
for m = 1:G.M
    T = reshape(C(:, m, 1:G.L(m)), B, G.L(m));
    N = sum(T(:));
    E = sum(T, 2) * sum(T, 1) / N;
    ok = E > 0;
    chi2 = sum((T(ok) - E(ok)).^2 ./ E(ok));
    st.V(m) = sqrt(chi2 / (N * (min(B, G.L(m)) - 1)));
end
end

function pb = psiBound(G, Sb, B)
% Upper bound on Psi: per group, the larger of (a) the most even integer
% split of each block's level slots and (b) the proportional split of the
% design's level totals (continuous relaxation). Both relax the problem.
LB = 0;
for m = 1:G.M
    L = G.L(m);
    g = G.gsize(m);
    wm = L ./ (Sb.^2 * g^2);
    n = Sb * g;
    q = floor(n / L);
    r = n - q * L;
    lbInt = sum(wm .* (r .* (q + 1).^2 + (L - r) .* q.^2)) - B;
    T = zeros(L, 1);
    for k = G.cols{m}
        T = T + accumarray(G.Z(:, k), 1, [L 1]);
    end
    lbCont = B * L * sum((T / sum(T)).^2) - B;
    LB = LB + max([lbInt, lbCont, 0]);
end
pb = 100 * (1 - sqrt(LB) / sqrt(B * sum(G.L - 1)));
end

function [Is, Ifull, P] = mnlInfo(Info, S)
if ~isfield(Info, 'X') || ~isfield(Info, 'beta')
    error('dce_block:info', 'Info needs fields X (S-by-J-by-P) and beta (P-by-1).');
end
Xi = double(Info.X);
[S2, J, P] = size(Xi);
if S2 ~= S
    error('dce_block:info', 'Info.X must have one row per task (%d).', S);
end
beta = double(Info.beta(:));
if numel(beta) ~= P
    error('dce_block:info', 'Info.beta must have %d elements.', P);
end
if isfield(Info, 'avail') && ~isempty(Info.avail)
    avail = logical(Info.avail);
else
    avail = true(S, J);
end
Is = zeros(P, P, S);
for s = 1:S
    Xs = reshape(Xi(s, :, :), J, P);
    a = avail(s, :)';
    v = Xs * beta;
    e = zeros(J, 1);
    e(a) = exp(v(a) - max(v(a)));
    pr = e / sum(e);
    Is(:, :, s) = Xs' * (diag(pr) - pr * pr') * Xs;
end
Ifull = sum(Is, 3);
end

function st = infoStats(Is, Ifull, y, B, P)
st.I = Ifull;
st.Derror = exp(-logdetPD(Ifull) / P);
st.share = zeros(1, B);
st.DerrorLOBO = zeros(1, B);
st.blockRank = zeros(1, B);
st.DerrorBlock = zeros(1, B);
for b = 1:B
    Ib = sum(Is(:, :, y == b), 3);
    st.share(b) = trace(Ifull \ Ib) / P;
    st.DerrorLOBO(b) = exp(-logdetPD(Ifull - Ib) / P);
    st.blockRank(b) = rank((Ib + Ib') / 2);
    st.DerrorBlock(b) = exp(-logdetPD(Ib) / P);
end
st.LOBOratio = st.DerrorLOBO / st.Derror;
end

function ld = logdetPD(A)
A = (A + A') / 2;
[Rc, flag] = chol(A);
if flag ~= 0 || rank(A) < size(A, 1)
    ld = -Inf;
else
    ld = 2 * sum(log(diag(Rc)));
end
end

function printSummary(out, S, K, B, grouped, G)
fprintf('\ndce_block: %d tasks, %d columns, %d blocks (sizes %s)\n', ...
    S, K, B, strtrim(sprintf('%d ', out.blockSizes)));
if out.search.evaluateOnly
    fprintf('  Evaluating a given blocking (no search)\n');
else
    fprintf('  Search: %d restarts x %d swaps, polish %d, seed %d, %.1f s\n', ...
        out.search.restarts, out.search.iterations, out.search.polish, ...
        out.search.seed, out.search.seconds);
end
if grouped
    fprintf('  Block level balance Psi = %.2f%% (pooled groups); within alternatives %.2f%%\n', ...
        out.Psi, out.PsiWithin);
    fprintf('  Relaxation bounds: pooled %.2f%%, within %.2f%%\n', out.PsiBound, out.PsiWithinBound);
else
    fprintf('  Block level balance Psi = %.2f%% (within alternatives); relaxation bound %.2f%%\n', ...
        out.Psi, out.PsiBound);
end
if ~out.designBalanced
    fprintf('  Note: the design itself is not level balanced, so 100%% is unattainable.\n');
end
nAv = nnz(out.absent.Avoidable);
fprintf('  Absent levels (block x column): %d, of which %d in blocks large enough to show all levels\n', ...
    height(out.absent), nAv);
if ~isempty(out.extremes) && height(out.extremes) > 0
    fprintf('  Extremes present in every block: %d of %d checks\n', ...
        nnz(out.extremes.HasLowest & out.extremes.HasHighest), height(out.extremes));
end
fprintf('  Ngene correlation criterion (for comparison): %.3f\n', out.corrNgene);
if ~isempty(out.info)
    fprintf('  MNL D-error (full design) %.4g; information shares by block: %s\n', ...
        out.info.Derror, strtrim(sprintf('%.3f ', out.info.share)));
    fprintf('  D-error without one block, relative to full: %s\n', ...
        strtrim(sprintf('%.3f ', out.info.LOBOratio)));
end
if ~out.search.evaluateOnly
    c = out.candidates;
    fprintf('  Distinct restart solutions: %d; at best balance: %d; selected by %s\n', ...
        height(c), nnz(c.AbsentRequired == min(c.AbsentRequired) & ...
        c.Psi >= max(c.Psi(c.AbsentRequired == min(c.AbsentRequired))) - 1e-9), out.search.select);
end
for j = 1:numel(G.warn)
    fprintf('  Warning: %s\n', G.warn{j});
end
end
