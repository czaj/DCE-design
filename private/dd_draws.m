function DR = dd_draws(C, nBayes, nRand, seed, method, nSim, nIS)
%DD_DRAWS Prior and random-coefficient draws for every model. Internal.
%   DR(m).R       number of Bayesian (outer) draws; 1 when the priors are fixed
%                 or the model's stat is 'fixed'
%   DR(m).B       P-by-R coefficient draws (means, for random coefficients)
%   DR(m).Sg      Pr-by-R spreads of the random coefficients
%   DR(m).Beta{r}, .dmu{r}, .dsig{r}   individual-level coefficients for outer
%                 draw r and their derivatives (rp and rppanel models)
%   DR(m).etaSim, .uSim   standardised coefficients and uniforms of nSim
%                 simulated respondents (rppanel; common random numbers)
%   DR(m).zIS, .nIS   with nIS > 0 (rppanel): two pools of standard normal
%                 points for the importance sampling of dd_panelis, which
%                 uses nIS/2 points per choice sequence from each pool
%   Draws are scrambled Sobol points ('sobol', default) or randomised Halton
%   points ('halton'), fully determined by seed.

if nargin < 5 || isempty(method)
    method = 'sobol';
end
if nargin < 6
    nSim = 0;
end
if nargin < 7
    nIS = 0;
end
if isscalar(nBayes)
    nBayes = repmat(nBayes, 1, C.M);
end
DR = struct('R', {}, 'B', {}, 'Sg', {}, 'Beta', {}, 'dmu', {}, 'dsig', {}, 'eta', {}, 'etaSim', {}, 'uSim', {}, 'zIS', {}, 'nIS', {});
for m = 1:C.M
    Mm = C.model(m);
    ib = find(Mm.bdist > 0);
    is = Mm.ridx(Mm.sdist(Mm.ridx) > 0);
    nd = numel(ib) + numel(is);
    d = struct('R', 1, 'B', Mm.b0, 'Sg', Mm.s0(Mm.ridx), 'Beta', {{}}, 'dmu', {{}}, 'dsig', {{}}, ...
        'eta', [], 'etaSim', [], 'uSim', [], 'zIS', [], 'nIS', 0);
    if nd > 0 && ~strcmp(Mm.stat, 'fixed')
        R = nBayes(m);
        u = uniformPoints(R, nd, seed + 7919 * m, method);
        B = repmat(Mm.b0, 1, R);
        Sg = repmat(Mm.s0, 1, R);
        for i = 1:numel(ib)
            B(ib(i), :) = fromPrior(u(:, i)', Mm.bdist(ib(i)), Mm.ba(ib(i)), Mm.bb(ib(i)));
        end
        for i = 1:numel(is)
            Sg(is(i), :) = max(fromPrior(u(:, numel(ib) + i)', Mm.sdist(is(i)), Mm.sa(is(i)), Mm.sb(is(i))), 1e-8);
        end
        d.R = R;
        d.B = B;
        d.Sg = Sg(Mm.ridx, :);
    end
    if ~strcmp(Mm.type, 'mnl')
        d.eta = standardise(uniformPoints(nRand, Mm.Pr, seed + 104729 * m + 1, method), Mm);
        d.Beta = cell(1, d.R);
        d.dmu = cell(1, d.R);
        d.dsig = cell(1, d.R);
        for r = 1:d.R
            [d.Beta{r}, d.dmu{r}, d.dsig{r}] = dd_betas(Mm, d.B(:, r), d.Sg(:, r), d.eta);
        end
        if strcmp(Mm.type, 'rppanel') && nSim > 0
            d.etaSim = standardise(uniformPoints(nSim, Mm.Pr, seed + 15485863 * m + 2, method), Mm);
            rs = RandStream('mt19937ar', 'Seed', mod(seed + 32452843 * m + 3, 2^32 - 1));
            d.uSim = rand(rs, nSim, C.S);
        end
        if strcmp(Mm.type, 'rppanel') && nIS > 0
            rz = RandStream('mt19937ar', 'Seed', mod(seed + 49979687 * m + 4, 2^32 - 1));
            d.zIS = randn(rz, Mm.Pr, 2^17, 2);
            d.nIS = max(8, 2 * floor(nIS / 4));
        end
    end
    DR(m) = d;
end
end

function eta = standardise(u, Mm)
n = size(u, 1);
eta = zeros(Mm.Pr, n);
for i = 1:Mm.Pr
    ui = u(:, i)';
    switch Mm.rdist(Mm.ridx(i))
        case {1, 2}
            eta(i, :) = -sqrt(2) * erfcinv(2 * ui);
        case 3
            eta(i, :) = 2 * ui - 1;
        case 4
            lowHalf = ui < 0.5;
            t = zeros(size(ui));
            t(lowHalf) = sqrt(2 * ui(lowHalf)) - 1;
            t(~lowHalf) = 1 - sqrt(2 * (1 - ui(~lowHalf)));
            eta(i, :) = t;
    end
end
end

function x = fromPrior(u, dist, a, b)
if dist == 1
    x = a - b * sqrt(2) * erfcinv(2 * u);
else
    x = a + (b - a) * u;
end
end

function u = uniformPoints(n, d, seed, method)
% n-by-d points in (0,1), deterministic given seed.
rs = RandStream('mt19937ar', 'Seed', mod(seed, 2^32 - 1));
if d == 0
    u = zeros(n, 0);
    return
end
if strcmpi(method, 'sobol') && exist('sobolset', 'file') == 2 && license('test','Statistics_Toolbox') && d <= 1111
    old = RandStream.setGlobalStream(rs);
    restore = onCleanup(@() RandStream.setGlobalStream(old));
    p = scramble(sobolset(d), 'MatousekAffineOwen');
    u = net(p, n);
    clear restore
else
    pr = primes(max(20, ceil(d * (log(d + 1) + 3) * 1.3)));
    pr = pr(1:d);
    u = zeros(n, d);
    for j = 1:d
        u(:, j) = mod(radicalInverse((1:n)', pr(j)) + rand(rs), 1);
    end
end
u = min(max(u, 1e-10), 1 - 1e-10);
end

function h = radicalInverse(i, base)
h = zeros(size(i));
f = 1 / base;
while any(i > 0)
    h = h + f * mod(i, base);
    i = floor(i / base);
    f = f / base;
end
end
