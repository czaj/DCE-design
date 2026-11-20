function [I, method, dg] = dd_panelinfo(Mm, X, d, r, exactMax)
%DD_PANELINFO Expected Fisher information of one respondent in a panel mixed logit. Internal.
%   X is J-by-P-by-T: the T choice tasks this respondent answers (a block).
%   d are the model's draws (dd_draws), r the index of the Bayesian draw.
%   I = E_y[ s(y) s(y)' ], s(y) the score of the probability of the whole
%   choice sequence y:
%       L(y) = E_beta[ prod_t p_t,y_t(beta) ].
%   The expectation over sequences is exact (all J^T sequences) when
%   J^T <= exactMax, otherwise it is a Monte Carlo average over the simulated
%   respondents in d.etaSim / d.uSim (the same ones for every design).
%
%   The expectation over the coefficients is taken
%   - by importance sampling around the mode of each sequence, when d.zIS
%     holds the points for it (dd_panelis; the default), or
%   - with the same draws d.Beta for all sequences ('prior'). That is the
%     plain simulated likelihood; it overstates the information of sequences
%     of more than a few choices unless the number of draws is very large.
%
%   dg (importance sampling only): diagnostics of dd_panelis.

[J, P, T] = size(X);
dg = [];
exact = J^T <= exactMax;
if exact
    method = 'exact';
    n = J^T;
    Y = zeros(n, T);
    idx = (0:n - 1)';
    for t = 1:T
        Y(:, t) = mod(floor(idx / J^(t - 1)), J) + 1;
    end
else
    method = 'simulated';
    if isempty(d.etaSim)
        dd_err('panel:draws', 'Simulated respondents are required for this panel evaluation (set PanelRespondents > 0).');
    end
    n = size(d.etaSim, 2);
    BetaS = dd_betas(Mm, d.B(:, r), d.Sg(:, r), d.etaSim);
    Vs = pagemtimes(X, reshape(BetaS, P, 1, 1, n));
    Vs = Vs - max(Vs, [], 1);
    ps = exp(Vs);
    ps = cumsum(ps ./ sum(ps, 1), 1);
    u = reshape(d.uSim(:, 1:T)', 1, 1, T, n);
    Y = reshape(1 + sum(u > ps(1:J - 1, :, :, :), 1), T, n)';
end

if isfield(d, 'zIS') && ~isempty(d.zIS)
    [I, dg] = dd_panelis(Mm, X, d.B(:, r), d.Sg(:, r), Y, exact, d.zIS, d.nIS, 0.8);
    return
end

% ---- the same draws for all sequences ------------------------------------------------
Beta = d.Beta{r};
dmu = d.dmu{r};
dsig = d.dsig{r};
Rr = size(Beta, 2);
V = pagemtimes(X, reshape(Beta, P, 1, 1, Rr));
V = V - max(V, [], 1);
p = exp(V);
p = p ./ sum(p, 1);
lp = reshape(log(max(p, realmin)), J, T, Rr);
Bsum = reshape(sum(sum(p .* X, 1), 3), P, Rr);
Xt = permute(X, [1 3 2]);

% in chunks of sequences, to bound memory
I = zeros(Mm.Pt);
chunk = max(1, floor(4e6 / Rr));
for c0 = 1:chunk:n
    ii = c0:min(n, c0 + chunk - 1);
    nc = numel(ii);
    LP = zeros(nc, Rr);
    A = zeros(nc, P);
    for t = 1:T
        LP = LP + reshape(lp(Y(ii, t), t, :), nc, Rr);
        A = A + reshape(Xt(Y(ii, t), t, :), nc, P);
    end
    W = exp(LP);
    L = mean(W, 2);
    dLm = (A .* (W * dmu') - W * (dmu .* Bsum)') / Rr;
    dLs = (A(:, Mm.ridx) .* (W * dsig') - W * (dsig .* Bsum(Mm.ridx, :))') / Rr;
    dL = [dLm dLs];
    if exact
        I = I + dL' * (dL ./ L);
    else
        sc = dL ./ L;
        I = I + (sc' * sc) / n;
    end
end
end
