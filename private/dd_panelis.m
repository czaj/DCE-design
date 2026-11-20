function [I, dg, S] = dd_panelis(Mm, X, mu, sg, Y, allSequences, Z, Rh, share)
%DD_PANELIS Panel mixed logit information by importance sampling. Internal.
%   X   J-by-P-by-T: the choice tasks one respondent answers.
%   mu  P-by-1 means, sg Pr-by-1 spreads of the random coefficients.
%   Y   n-by-T choice sequences. allSequences true: Y lists every possible
%       sequence and each is weighted by its probability; false: Y are the
%       sequences of n simulated respondents, equally weighted.
%   Z   Pr-by-M-by-2 standard normal points: two independent pools. Sequence
%       number i takes Rh consecutive points of each pool, starting at
%       (i-1)*Rh+1 (wrapping around), so that sequences have their own points
%       and their simulation errors average out over sequences.
%   Rh  draws per sequence and set.
%   share  share of each set drawn from the proposal around the mode (the
%       rest is drawn from the distribution of the coefficients).
%
%   The information is I = E_y[ s(y) s(y)' ], with s(y) the score of the
%   probability of the whole sequence, L(y) = E_beta[ prod_t p_t,y_t(beta) ].
%   s(y) is an expectation over the distribution of the coefficients given y.
%   After several choices that distribution is far narrower than the
%   distribution of the coefficients in the population, so draws from the
%   latter leave few effective draws per sequence, the simulated score is
%   noisy, and its outer product overstates the information (the D-error is
%   understated, the more the longer the sequence).
%
%   Here every sequence gets its own draws.
%   1. The coefficients are written as functions of standard normal z. The
%      mode of the density of z given y is found by Newton steps, and z is
%      drawn from a mixture of the normal approximation at the mode (Laplace)
%      and the standard normal itself, which keeps the weights bounded.
%   2. With u = d log prod_t p_t / dz, any smooth v(z) satisfies
%      E[ div v + v'(u - z) | y ] = 0 (integration by parts). The score is
%      written with the v that makes its integrand constant when the density
%      of z given y is exactly normal: for a normal or lognormal coefficient
%          d/d mean   = E[ z_i + c_i | y ] / spread_i
%          d/d spread = E[ z_i^2 + V_ii + (z_i + m_i) c_i - 1 | y ] / spread_i
%      with c = V (u - z), m the mode and V the Laplace covariance; for the
%      other coefficients g_k + (dg_k/dz at the mode) c, g = d log prod p / d beta.
%      What is left to simulate is the departure from normality.
%   3. The draws around the mode come in antithetic pairs, and the bias of the
%      weighted mean (a ratio of two means) is removed to first order by a
%      jackknife over two halves of the draws.
%   4. The score is computed twice, from two independent sets of draws, and
%      the information is the average of the cross products s_A s_B', whose
%      expectation does not contain the variance of the simulation noise.
%
%   I is symmetric but, being a difference, need not be positive definite: the
%   information of one block alone is often singular, and the simulation noise
%   can then make an eigenvalue negative. Callers add up the blocks first.
%   dg: .ess  mean effective share of the draws,
%       .noise  trace of the removed noise term relative to the trace of I,
%       .pooled  the estimate from the mean of the two scores: positive
%                semidefinite, but it contains the variance of the simulation
%                noise. For the caller to fall back on if the sum of I over
%                blocks is not positive definite.
%   S (Pt-by-n-by-2) are the scores of the sequences from the two sets of draws.

[J, P, T] = size(X);
n = size(Y, 1);
Pr = Mm.Pr;
Pt = Mm.Pt;
ridx = Mm.ridx(:)';
JT = J * T;
X2 = reshape(permute(X, [1 3 2]), JT, P);      % rows (j, t), j fastest
Xr = X2(:, ridx);
fix = setdiff(1:P, ridx);
v0 = X2(:, fix) * mu(fix);
rd = Mm.rdist(ridx);
rsg = Mm.rsign(ridx);
mur = mu(ridx);
sg = max(sg(:), 1e-8);
M = size(Z, 2);
Rh = max(8, 2 * floor(Rh / 2));
R1 = max(4, 4 * floor(min(Rh - 2, round(share * Rh)) / 4));    % a multiple of 4: antithetic pairs, split in two halves
R0 = Rh - R1;
h1 = [1:R1 / 4, R1 / 2 + (1:R1 / 4), R1 + (1:R0 / 2)];          % half of the draws, pairs kept together
h2 = setdiff(1:Rh, h1);
a1 = R1 / Rh;
a0 = 1 - a1;
eyeP = eye(Pr);

S = zeros(Pt, n, 2);
Lh = zeros(2, n);
essSum = 0;
chunk = max(1, floor(3e6 / (JT * Rh)));
for c0 = 1:chunk:n
    ii = c0:min(n, c0 + chunk - 1);
    nc = numel(ii);
    Ych = zeros(JT, nc);
    lin = ((0:T - 1) * J + Y(ii, :)) + JT * (0:nc - 1).';
    Ych(lin(:)) = 1;

    % ---- mode of the density of z given the sequence: damped Newton -------------
    z = zeros(Pr, nc);
    [b, dbz, ~, eta] = zmap(z, mur, sg, rd, rsg);
    [p, lw] = probs(v0, Xr, b, J, T, Ych);
    f = lw - 0.5 * sum(z.^2, 1);
    for it = 1:15
        g = dbz .* (Xr.' * (Ych - p)) - z;
        H = gnhess(p, dbz, Xr, J, T, eyeP);
        step = reshape(pagemldivide(H, reshape(g, Pr, 1, nc)), Pr, nc);
        t = ones(1, nc);
        for h = 1:8
            zn = z + t .* step;
            [bn, dbzn, ~, etan] = zmap(zn, mur, sg, rd, rsg);
            [pn, lwn] = probs(v0, Xr, bn, J, T, Ych);
            fn = lwn - 0.5 * sum(zn.^2, 1);
            ok = fn >= f - 1e-10;
            if all(ok)
                break
            end
            t(~ok) = t(~ok) / 2;
        end
        moved = max(abs(t .* step) .* ok, [], 'all');
        z(:, ok) = zn(:, ok);
        dbz(:, ok) = dbzn(:, ok);
        eta(:, ok) = etan(:, ok);
        p(:, ok) = pn(:, ok);
        f(ok) = fn(ok);
        if moved < 1e-5
            break
        end
    end

    % ---- at the mode --------------------------------------------------------------
    g0 = X2.' * (Ych - p);                                    % d log prod p / d beta, P-by-nc
    p4 = reshape(p, J, T, 1, nc);
    X4 = reshape(X2, J, T, P);
    Zf = reshape(sqrt(p4) .* (X4 - sum(p4 .* X4, 1)), JT, P, nc);
    Hpr = pagemtimes(Zf, 'transpose', Zf(:, ridx, :), 'none');     % minus d g / d beta_random, P-by-Pr-by-nc
    dz = reshape(dbz, 1, Pr, nc);
    H = Hpr(ridx, :, :) .* dz .* pagetranspose(dz) + eyeP;
    [Lc, ld] = pchol(H);
    Lt = pagetranspose(Lc);
    Ci = pagemldivide(Lt, repmat(eyeP, 1, 1, nc));            % Ci*Ci.' = inv(H)
    Vz = pagemtimes(Ci, 'none', Ci, 'transpose');             % Laplace covariance of z
    BV = pagemtimes(-Hpr .* dz, Vz);                          % (d g / d z) * Vz, P-by-Pr-by-nc
    zc = reshape(z, Pr, 1, nc);
    Yr = reshape(repmat(reshape(Ych, JT, 1, nc), 1, Rh, 1), JT, Rh * nc);

    % ---- two independent sets of draws ------------------------------------------
    cols = reshape((mod((ii(:) - 1) * Rh + (0:Rh - 1), M) + 1).', 1, []);
    for k = 1:2
        E = reshape(Z(:, cols, k), Pr, Rh, nc);               % the points of each sequence
        EL = E(:, 1:R1, :);
        EL(:, R1 / 2 + 1:R1, :) = -EL(:, 1:R1 / 2, :);         % antithetic pairs
        zL = zc + pagemtimes(Ci, EL);                         % around the mode
        z0 = E(:, R1 + 1:Rh, :);                              % from the standard normal
        q0 = sum(pagemtimes(Lt, z0 - zc).^2, 1);
        zz = cat(2, zL, z0);                                  % Pr-by-Rh-by-nc
        quad = cat(2, sum(EL.^2, 1), q0);
        tq = 0.5 * reshape(ld, 1, 1, nc) - 0.5 * quad + 0.5 * sum(zz.^2, 1);   % log(q_mode / phi)
        [b, dbzd, ~, etad] = zmap(reshape(zz, Pr, Rh * nc), mur, sg, rd, rsg);
        [p, lw] = probs(v0, Xr, b, J, T, Yr);
        G = X2.' * (Yr - p);                                  % P-by-(Rh*nc)
        dU = reshape(dbzd .* G(ridx, :), Pr, Rh, nc) - zz;    % u - z
        Cc = pagemtimes(Vz, dU);                              % Pr-by-Rh-by-nc
        G3 = reshape(G, P, Rh, nc);
        A = zeros(Pt, Rh, nc);
        A(1:P, :, :) = G3 + pagemtimes(BV, dU);
        for i = 1:Pr
            zi = zz(i, :, :);
            ci = Cc(i, :, :);
            if rd(i) <= 2
                A(ridx(i), :, :) = (zi + ci) / sg(i);
                A(P + i, :, :) = (zi.^2 + Vz(i, i, :) + (zi + zc(i, 1, :)) .* ci - 1) / sg(i);
            else
                gi = G3(ridx(i), :, :);
                A(P + i, :, :) = reshape(etad(i, :), 1, Rh, nc) .* gi ...
                    + reshape(dbz(i, :) / sg(i) .* g0(ridx(i), :), 1, 1, nc) .* ci ...
                    + reshape(eta(i, :), 1, 1, nc) .* (A(ridx(i), :, :) - gi);
            end
        end
        tm = max(tq, 0);
        lwt = reshape(lw, 1, Rh, nc) - tm - log(a1 * exp(tq - tm) + a0 * exp(-tm));
        mx = max(lwt, [], 2);
        w = exp(lwt - mx);
        sw = sum(w, 2);
        % the ratio of two means is biased in 1/draws; the jackknife over the two halves removes that term
        sw1 = sum(w(:, h1, :), 2);
        g1 = sum(A(:, h1, :) .* w(:, h1, :), 2);
        g2 = sum(A(:, h2, :) .* w(:, h2, :), 2);
        sk = 2 * (g1 + g2) ./ sw - 0.5 * (g1 ./ sw1 + g2 ./ (sw - sw1));
        S(:, ii, k) = reshape(sk, Pt, nc);
        Lh(k, ii) = reshape(exp(mx) .* sw / Rh, 1, nc);
        essSum = essSum + sum(sw.^2 ./ sum(w.^2, 2), 'all') / Rh;
    end
end

if allSequences
    wy = mean(Lh, 1);
    wy = wy / sum(wy);
else
    wy = ones(1, n) / n;
end
SA = S(:, :, 1);
SB = S(:, :, 2);
I = (SA .* wy) * SB.';
I = (I + I.') / 2;
Sm = (SA + SB) / 2;
Ip = (Sm .* wy) * Sm.';
dg = struct('ess', essSum / (2 * n), 'noise', (trace(Ip) - trace(I)) / max(trace(I), realmin), 'pooled', Ip);
end

function [b, dbz, dm, eta] = zmap(z, mur, sg, rd, rsg)
% Coefficients as functions of standard normal z: b, db/dz, db/dmean, and
% eta with db/dspread = dm .* eta. One row per random coefficient.
b = zeros(size(z));
dbz = b;
dm = ones(size(z));
eta = z;
for i = 1:numel(rd)
    zi = z(i, :);
    switch rd(i)
        case 1          % normal
            b(i, :) = mur(i) + sg(i) * zi;
            dbz(i, :) = sg(i);
        case 2          % lognormal
            e = rsg(i) * exp(mur(i) + sg(i) * zi);
            b(i, :) = e;
            dbz(i, :) = e * sg(i);
            dm(i, :) = e;
        case 3          % uniform on mean -/+ spread
            eta(i, :) = erf(zi / sqrt(2));
            b(i, :) = mur(i) + sg(i) * eta(i, :);
            dbz(i, :) = sg(i) * 2 * exp(-0.5 * zi.^2) / sqrt(2 * pi);
        case 4          % triangular on mean -/+ spread
            u = 0.5 * erfc(-zi / sqrt(2));
            uc = 0.5 * erfc(zi / sqrt(2));
            lo = zi < 0;
            e = zeros(size(zi));
            de = e;
            e(lo) = sqrt(2 * u(lo)) - 1;
            de(lo) = 1 ./ sqrt(2 * max(u(lo), realmin));
            e(~lo) = 1 - sqrt(2 * uc(~lo));
            de(~lo) = 1 ./ sqrt(2 * max(uc(~lo), realmin));
            eta(i, :) = e;
            b(i, :) = mur(i) + sg(i) * e;
            dbz(i, :) = sg(i) * de .* exp(-0.5 * zi.^2) / sqrt(2 * pi);
    end
end
end

function [p, lw] = probs(v0, Xr, b, J, T, Ych)
% Choice probabilities of every alternative of every task (JT-by-m) and the
% log-probability of the sequence, for m coefficient vectors.
m = size(b, 2);
V = reshape(v0 + Xr * b, J, T * m);
V = V - max(V, [], 1);
p = exp(V);
p = reshape(p ./ sum(p, 1), J * T, m);
lw = sum(Ych .* log(max(p, realmin)), 1);
end

function H = gnhess(p, dbz, Xr, J, T, eyeP)
% Minus the Hessian of the log density of z given the sequence, without the
% term in the second derivative of the coefficients (exact for normal
% coefficients; positive definite always).
m = size(p, 2);
Pr = size(Xr, 2);
p4 = reshape(p, J, T, 1, m);
X4 = reshape(Xr, J, T, Pr);
Zm = reshape(sqrt(p4) .* (X4 - sum(p4 .* X4, 1)), J * T, Pr, m);
d = reshape(dbz, Pr, 1, m);
H = pagemtimes(Zm, 'transpose', Zm, 'none') .* d .* pagetranspose(d) + eyeP;
end

function [L, ld] = pchol(H)
% Cholesky factors (lower, H = L*L.') and log-determinants of a stack of
% positive definite matrices.
Pr = size(H, 1);
m = size(H, 3);
L = zeros(Pr, Pr, m);
ld = zeros(1, m);
for k = 1:Pr
    if k == 1
        s = H(1, 1, :);
    else
        s = H(k, k, :) - sum(L(k, 1:k - 1, :).^2, 2);
    end
    s = sqrt(max(s, 1e-12));
    L(k, k, :) = s;
    ld = ld + 2 * log(reshape(s, 1, m));
    if k < Pr
        if k == 1
            L(2:Pr, 1, :) = H(2:Pr, 1, :) ./ s;
        else
            L(k + 1:Pr, k, :) = (H(k + 1:Pr, k, :) - pagemtimes(L(k + 1:Pr, 1:k - 1, :), pagetranspose(L(k, 1:k - 1, :)))) ./ s;
        end
    end
end
end
