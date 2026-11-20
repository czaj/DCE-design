function [cand, log] = dd_search(C, DR, o, D0, nRestarts, pool)
%DD_SEARCH Multi-start search. Internal.
%   Each start: a random feasible design (start 1 uses D0 when given), a local
%   optimisation, then o.Kicks rounds of perturbation and re-optimisation.
%   Start r draws from substream r of one generator, so the result does not
%   depend on whether the starts run in parallel.
%
%   The local optimisation follows o.method:
%     'swap'         exchanges within columns and level changes (dd_localsearch)
%     'rowexchange'  replacement of whole tasks from the candidate pool
%     'both'         the two in turn until neither finds an improvement
%   cand: distinct designs sorted by objective (fields D, f). log: per start.

if nargin < 5 || isempty(nRestarts)
    nRestarts = o.Restarts;
end
if nargin < 6
    pool = [];
end
nW = 0;
if isequal(o.Parallel, true) || (ischar(o.Parallel) && strcmpi(o.Parallel, 'on'))
    p = gcp();
    nW = p.NumWorkers;
elseif ischar(o.Parallel) && strcmpi(o.Parallel, 'auto')
    p = gcp('nocreate');
    if ~isempty(p)
        nW = p.NumWorkers;
    end
end
Dall = cell(nRestarts, 1);
fall = inf(nRestarts, 1);
st = cell(nRestarts, 1);
ck = o.Checkpoint;
parfor (r = 1:nRestarts, nW)
    rs = RandStream('mrg32k3a', 'Seed', o.Seed); %#ok<PFBNS>
    rs.Substream = r;
    tStart = tic;
    if r == 1 && ~isempty(D0)
        D = D0;
    else
        D = dd_initial(C, rs);
    end
    [D, f, s] = optimise(C, D, DR, o, rs, pool);
    kicks = 0;
    for q = 1:o.Kicks
        if posixtime(datetime('now')) > o.deadline
            break
        end
        D2 = dd_kick(C, D, rs, o.KickSize);
        [D2, f2, s2] = optimise(C, D2, DR, o, rs, pool);
        s.evaluations = s.evaluations + s2.evaluations;
        s.passes = s.passes + s2.passes;
        s.moves = s.moves + s2.moves;
        s.rowMoves = s.rowMoves + s2.rowMoves;
        s.timedOut = s.timedOut || s2.timedOut;
        kicks = kicks + 1;
        if f2 < f * (1 - o.Tol)
            D = D2;
            f = f2;
        end
    end
    s.kicks = kicks;
    s.seconds = toc(tStart);
    s.objective = f;
    Dall{r} = D;
    fall(r) = f;
    st{r} = s;
    if ~isempty(ck)
        saveCheckpoint(ck, r, D, f);
    end
end
log = [st{:}];
[fs, ord] = sort(fall);
cand = struct('D', {}, 'f', {});
seen = {};
for i = 1:nRestarts
    if ~isfinite(fs(i))
        continue
    end
    Dsorted = sortrows(Dall{ord(i)});
    if any(cellfun(@(q) isequal(q, Dsorted), seen))
        continue
    end
    seen{end + 1} = Dsorted; %#ok<AGROW>
    cand(end + 1) = struct('D', Dall{ord(i)}, 'f', fs(i)); %#ok<AGROW>
end
if isempty(cand)
    dd_err('search:notIdentified', ['No start produced a design that identifies all parameters. ' ...
        'Check for parameters without variation across alternatives (for example a main effect of a scenario variable), ' ...
        'for attributes that are perfectly related by the constraints, or increase the number of rows.']);
end
end

function [D, f, s] = optimise(C, D, DR, o, rs, pool)
s = struct('passes', 0, 'evaluations', 0, 'moves', 0, 'timedOut', false, 'rowMoves', 0);
useSwap = ~strcmp(o.method, 'rowexchange');
useRow = ~isempty(pool) && ~strcmp(o.method, 'swap');
f = Inf;
if useSwap
    [D, f, a] = dd_localsearch(C, D, DR, o, rs);
    s = accumulate(s, a, false);
end
if useRow
    for round = 1:25
        [D, f, b] = dd_rowexchange(C, D, DR, o, rs, pool);
        s = accumulate(s, b, true);
        if b.moves == 0 || ~useSwap || s.timedOut
            break
        end
        [D, f, a] = dd_localsearch(C, D, DR, o, rs);
        s = accumulate(s, a, false);
        if a.moves == 0 || s.timedOut
            break
        end
    end
end
end

function s = accumulate(s, a, isRow)
s.passes = s.passes + a.passes;
s.evaluations = s.evaluations + a.evaluations;
s.moves = s.moves + a.moves;
s.timedOut = s.timedOut || a.timedOut;
if isRow
    s.rowMoves = s.rowMoves + a.moves;
end
end

function saveCheckpoint(file, r, D, f)
[p, n, e] = fileparts(file);
if isempty(e), e = '.mat'; end
fn = fullfile(p, sprintf('%s_start%03d%s', n, r, e));
save(fn, 'D', 'f', 'r');
end
