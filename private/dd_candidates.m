function cand = dd_candidates(C, DR, o)
%DD_CANDIDATES Candidate choice tasks for the row-exchange moves. Internal.
%   All feasible tasks when they can be enumerated (at most o.CandidateMax),
%   otherwise a random sample of feasible tasks, the same for a given seed.
%   For every active model the coded regressors of the candidates are stored,
%   and their information matrices: for all prior draws when that fits in
%   memory, otherwise for the first draws only (used to rank candidates; the
%   best are then evaluated on all draws).
%
%   cand.D        nCand-by-K level indices
%   cand.all      true if every feasible task is a candidate
%   cand.X{m}     J-by-P-by-nCand
%   cand.I{m}     Pt-by-Pt-by-nCand-by-R or -by-numel(rsub{m})
%   cand.rsub{m}  prior draws that cand.I{m} covers
%   cand.full(m)  cand.I{m} covers all draws

K = C.K;
nFull = prod(C.L);
rs = RandStream('mt19937ar', 'Seed', mod(o.Seed + 7001, 2^32 - 1));
if nFull <= 20 * o.CandidateMax
    g = cell(1, K);
    for k = 1:K
        g{k} = 1:C.L(k);
    end
    [g{:}] = ndgrid(g{:});
    D = cell2mat(cellfun(@(v) v(:), g, 'UniformOutput', false));
    D = D(dd_rowok(C, D), :);
    cand.all = size(D, 1) <= o.CandidateMax;
    if ~cand.all
        D = D(sort(randperm(rs, size(D, 1), o.CandidateMax)), :);
    end
else
    D = zeros(0, K);
    tries = 0;
    while size(D, 1) < o.CandidateMax && tries < 200
        tries = tries + 1;
        Dr = zeros(4 * o.CandidateMax, K);
        for k = 1:K
            Dr(:, k) = randi(rs, C.L(k), size(Dr, 1), 1);
        end
        D = unique([D; Dr(dd_rowok(C, Dr), :)], 'rows', 'stable');
    end
    D = D(1:min(end, o.CandidateMax), :);
    cand.all = false;
end
cand.D = D;
n = size(D, 1);
cand.hash = D * C.hashW;
cand.X = cell(1, C.M);
cand.I = cell(1, C.M);
cand.rsub = cell(1, C.M);
cand.full = true(1, C.M);
for m = find(C.w > 0)'
    Mm = C.model(m);
    cand.X{m} = dd_X(Mm, D);
    R = DR(m).R;
    if n * R * Mm.Pt^2 <= 2.5e7
        cand.rsub{m} = 1:R;
        cand.I{m} = dd_rowinfo(Mm, cand.X{m}, DR(m));
    else
        cand.full(m) = false;
        rsub = 1:min(R, max(o.ScreenDraws, floor(2.5e7 / (n * Mm.Pt^2))));
        cand.rsub{m} = rsub;
        cand.I{m} = dd_rowinfo(Mm, cand.X{m}, subDraws(DR(m), rsub));
    end
end
end

function d = subDraws(d, rsub)
d.R = numel(rsub);
d.B = d.B(:, rsub);
if ~isempty(d.Sg), d.Sg = d.Sg(:, rsub); end
if ~isempty(d.Beta)
    d.Beta = d.Beta(rsub);
    d.dmu = d.dmu(rsub);
    d.dsig = d.dsig(rsub);
end
end
