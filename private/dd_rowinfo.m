function I = dd_rowinfo(Mm, X, d)
%DD_ROWINFO Expected Fisher information of single choice tasks. Internal.
%   X is J-by-P-by-n (coded regressors of n tasks), d the draws of the model
%   (see dd_draws). I is Pt-by-Pt-by-n-by-R: the information one respondent
%   contributes by answering each task, for each of the R prior draws.
%
%   mnl:  I = sum_j p_j (x_j - xbar)(x_j - xbar)'
%   rp (cross-sectional mixed logit; also the additive surrogate of rppanel):
%         I = sum_j g_j g_j' / Pbar_j,  g_j = d Pbar_j / d theta,
%         Pbar_j = mean over random draws of the logit probability.

[J, P, n] = size(X);
if strcmp(Mm.type, 'mnl')
    R = d.R;
    V = pagemtimes(X, reshape(d.B, P, 1, 1, R));
    V = V - max(V, [], 1);
    p = exp(V);
    p = p ./ sum(p, 1);
    Z = sqrt(p) .* (X - sum(p .* X, 1));
    I = pagemtimes(Z, 'transpose', Z, 'none');
else
    I = zeros(Mm.Pt, Mm.Pt, n, d.R);
    for r = 1:d.R
        Rr = size(d.Beta{r}, 2);
        V = pagemtimes(X, reshape(d.Beta{r}, P, 1, 1, Rr));
        V = V - max(V, [], 1);
        p = exp(V);
        p = p ./ sum(p, 1);
        G0 = p .* (X - sum(p .* X, 1));
        Gm = mean(G0 .* reshape(d.dmu{r}, 1, P, 1, Rr), 4);
        Gs = mean(G0(:, Mm.ridx, :, :) .* reshape(d.dsig{r}, 1, Mm.Pr, 1, Rr), 4);
        Z = [Gm, Gs] ./ sqrt(mean(p, 4));
        I(:, :, :, r) = pagemtimes(Z, 'transpose', Z, 'none');
    end
end
if J < 2
    I(:) = 0;
end
end
