function [Beta, dmu, dsig] = dd_betas(Mm, mu, sg, eta)
%DD_BETAS Individual-level coefficients of a mixed logit model. Internal.
%   mu P-by-1 means, sg Pr-by-1 spreads, eta Pr-by-R standardised draws.
%   Beta P-by-R; dmu P-by-R = dBeta/dmean; dsig Pr-by-R = dBeta/dspread.
%   normal, uniform, triangular: beta = mean + spread*eta
%   lognormal: beta = sign*exp(mean + spread*eta), mean and spread being those
%   of the underlying normal.

R = size(eta, 2);
Beta = repmat(mu, 1, R);
dmu = ones(Mm.P, R);
dsig = zeros(Mm.Pr, R);
for i = 1:Mm.Pr
    p = Mm.ridx(i);
    if Mm.rdist(p) == 2
        e = Mm.rsign(p) * exp(mu(p) + sg(i) * eta(i, :));
        Beta(p, :) = e;
        dmu(p, :) = e;
        dsig(i, :) = e .* eta(i, :);
    else
        Beta(p, :) = mu(p) + sg(i) * eta(i, :);
        dsig(i, :) = eta(i, :);
    end
end
end
