function report = dce_selftest(varargin)
%DCE_SELFTEST Check that the DCE-design toolbox computes what it should.
%
%   DCE_SELFTEST runs the test suite and prints one PASS/FAIL line per check;
%   it raises an error if any check fails. report = DCE_SELFTEST returns the
%   checks as a table. DCE_SELFTEST('quick') skips the slower checks.
%
%   What is checked:
%     1. Published designs with published D-errors (Bliemer, Rose & Chorus
%        2017; Hole's dcreate examples; Zwerina, Huber & Kuhfeld), including
%        the shares of tasks with a dominated alternative.
%     2. Information matrices of the multinomial logit, the cross-sectional
%        and the panel mixed logit against numerical second derivatives of
%        the expected log-likelihood; the simulated panel information against
%        the exact one; the panel information by importance sampling against
%        plain simulation with very many draws, also for a long panel.
%     3. The predicted standard errors against a Monte Carlo study that
%        simulates choices and re-estimates the model.
%     4. The search: global optimum of a small problem found by enumeration,
%        constraints, level frequencies, duplicates, reproducibility, the
%        running objective against a fresh evaluation.
%     5. The Ngene translator, blocking, export and the error messages.

quick = any(strcmpi(varargin, 'quick'));
T = struct('group', {}, 'check', {}, 'pass', {}, 'detail', {});
t0 = tic;

% =========================== 1. published benchmarks =========================
% Bliemer, Rose & Chorus (2017), Tables 3-4: time and cost linear,
% beta = (-0.2, -1.2), 8 tasks, two alternatives; pairs are (time, cost).
sp = struct();
sp.alternatives = {'a', 'b'};
sp.rows = 8;
sp.attributes.time = [10 15 20 25];
sp.attributes.cost = [1 2 3 4];
sp.model.utility.a = 'b_t*time + b_c*cost';
sp.model.utility.b = 'b_t*time + b_c*cost';
sp.model.priors.b_t = -0.2;
sp.model.priors.b_c = -1.2;
brc = { ...
    'O1', [15 3 15 2; 25 1 25 4; 20 4 20 1; 10 2 15 3; 10 3 10 2; 20 1 25 1; 15 2 20 4; 25 4 10 3], 0.304, 8
    'O2', [20 2 10 4; 10 4 25 1; 25 1 15 3; 10 1 20 2; 20 3 15 2; 15 2 10 1; 25 4 20 3; 15 3 25 4], 0.076, 5
    'E1', [10 4 25 1; 25 1 10 4; 10 3 25 2; 25 2 10 3; 20 1 15 4; 15 4 20 1; 20 3 15 2; 15 2 20 3], 0.057, 2
    'E2', [10 4 20 1; 25 2 10 4; 15 4 20 2; 20 1 15 3; 15 3 25 2; 25 2 10 3; 10 3 25 1; 20 1 15 4], 0.064, 0};
for i = 1:size(brc, 1)
    r = dce_evaluate(sp, brc{i, 2}, 'Verbose', 0);
    T = add(T, 'benchmark', sprintf('Bliemer-Rose-Chorus %s: D-error %.3f', brc{i, 1}, brc{i, 3}), ...
        abs(r.models.Derror - brc{i, 3}) < 6e-4, sprintf('%.5f', r.models.Derror));
    T = add(T, 'benchmark', sprintf('Bliemer-Rose-Chorus %s: %d of 8 tasks with a dominated alternative', brc{i, 1}, brc{i, 4}), ...
        size(r.checks.violations, 1) == brc{i, 4}, sprintf('%d', size(r.checks.violations, 1)));
end

% Hole, dcreate slides: wait {0,1,2,5}, flex {0,1}, thoro {0,1}; 8 sets of 2;
% published D-efficiency = 1 / D-error. Rows are (wait, flex, thoro) twice.
sp = struct();
sp.alternatives = {'a', 'b'};
sp.rows = 8;
sp.attributes.wait = [0 1 2 5];
sp.attributes.flex = [0 1];
sp.attributes.thoro = [0 1];
sp.dominance = false;
sp.unique = false;
u = 'b_w*wait + b_f*flex + b_t*thoro';
sp.model.utility.a = u;
sp.model.utility.b = u;
sp.model.priors = struct('b_w', -0.13, 'b_f', 0.2, 'b_t', 1);
da = [2 0 0 5 0 0; 2 0 1 0 1 0; 5 0 0 5 0 1; 1 1 0 0 0 0; 5 1 0 1 0 0; 1 0 1 2 0 1; 2 1 1 1 1 1; 0 1 0 1 0 1];
db = [5 0 1 0 1 0; 0 0 0 5 1 1; 0 0 1 5 1 0; 0 1 0 5 0 1; 5 0 0 0 1 1; 5 1 1 0 0 0; 5 0 1 0 1 0; 5 1 0 0 0 1];
dc = [1 1 0 0 0 1; 0 1 0 5 0 1; 1 0 1 5 1 1; 0 1 0 2 0 1; 0 0 1 2 1 0; 5 1 1 1 0 0; 2 0 0 1 1 1; 5 0 0 2 1 1];
r = dce_evaluate(sp, da, 'Verbose', 0);
T = add(T, 'benchmark', 'Hole dcreate, design a: D-efficiency 1.3214859932', abs(1 / r.models.Derror - 1.3214859932) < 1e-8, sprintf('%.10f', 1 / r.models.Derror));
r = dce_evaluate(sp, db, 'Verbose', 0);
T = add(T, 'benchmark', 'Hole dcreate, design b: D-efficiency 4.2422926121', abs(1 / r.models.Derror - 4.2422926121) < 1e-8, sprintf('%.10f', 1 / r.models.Derror));
sp.model.utility.a = 'b_w.dummy*wait + b_f*flex + b_t*thoro';
sp.model.utility.b = sp.model.utility.a;
sp.model.priors = struct('b_w', [-0.13 -0.26 -0.65], 'b_f', 0.2, 'b_t', 1);
r = dce_evaluate(sp, dc, 'Verbose', 0);
T = add(T, 'benchmark', 'Hole dcreate, design c (dummy coded): D-efficiency 0.8952664336', abs(1 / r.models.Derror - 0.8952664336) < 1e-8, sprintf('%.10f', 1 / r.models.Derror));

% Zwerina, Huber & Kuhfeld: 3^3, 9 sets of 3, effects coding, beta = 0;
% D-error = 27^(-1/2).
sp = struct();
sp.alternatives = {'a', 'b', 'c'};
sp.rows = 9;
sp.attributes.x1 = [1 2 3];
sp.attributes.x2 = [1 2 3];
sp.attributes.x3 = [1 2 3];
sp.dominance = false;
u = 'b1.effects*x1 + b2.effects*x2 + b3.effects*x3';
sp.model.utility = struct('a', u, 'b', u, 'c', u);
z = [111 222 333; 122 233 311; 133 211 322; 213 321 132; 221 332 113; 232 313 121; 312 123 231; 323 131 212; 331 112 223];
dz = zeros(9, 9);
for i = 1:9
    for j = 1:3
        dz(i, (j - 1) * 3 + (1:3)) = [floor(z(i, j) / 100), mod(floor(z(i, j) / 10), 10), mod(z(i, j), 10)];
    end
end
r = dce_evaluate(sp, dz, 'Verbose', 0);
T = add(T, 'benchmark', 'Zwerina-Huber-Kuhfeld, effects coding, zero priors: D-error 27^(-1/2)', abs(r.models.Derror - 27^-0.5) < 1e-10, sprintf('%.8f', r.models.Derror));
T = add(T, 'benchmark', 'Zwerina-Huber-Kuhfeld: all variances 0.222', all(abs(diag(r.models.avc) - 2 / 9) < 1e-10), sprintf('%.4f', max(diag(r.models.avc))));

% =========================== 2. information matrices ==========================
[specH, VH] = hotel();
C = dd_compile(setfield(specH, 'balance', 'none')); %#ok<SFLD>
o = dd_options(struct());
D = dd_tolevels(C, VH);
Mm = C.model(1);
X = dd_X(Mm, D);
d1 = struct('R', 1, 'B', Mm.b0);
I = sum(dd_rowinfo(Mm, X, d1), 3);
H = numHess(@(b) ellMNL(X, Mm.b0, b), Mm.b0);
T = add(T, 'information', 'MNL information = minus Hessian of the expected log-likelihood', relErr(-H, I) < 1e-5, sprintf('rel. error %.1e', relErr(-H, I)));
bq = Mm.b0 .* [1.3 0.6 1.2 0.8 1.1 0.9 1]' + [0 0 0 0 0 0 0.004]';
I2 = sum(dd_rowinfo(Mm, X, struct('R', 1, 'B', bq)), 3);
H2 = numHess(@(b) ellMNL(X, bq, b), bq);
T = add(T, 'information', 'MNL information at another parameter point (with the interaction)', relErr(-H2, I2) < 1e-5, sprintf('rel. error %.1e', relErr(-H2, I2)));

% mixed logit: normal, lognormal, uniform and triangular coefficients
spM = struct();
spM.alternatives = {'a', 'b', 'sq'};
spM.rows = 6;
spM.attributes.x = [0 1 2];
spM.attributes.q = [0 1];
spM.attributes.cost = [1 2 3];
spM.balance = 'none';
spM.dominance = false;
spM.model.type = 'rp';
spM.model.utility.a = 'b_x*x + b_q*q + b_c*cost';
spM.model.utility.b = 'b_x*x + b_q*q + b_c*cost';
spM.model.utility.sq = 'asc';
spM.model.priors = struct('b_x', 0.5, 'b_q', 0.4, 'b_c', -0.6, 'asc', 0.3);
spM.model.random = struct('b_x', 'normal(0.4)', 'b_q', 'uniform(0.5)', 'b_c', '-lognormal(0.5)', 'asc', 'triangular(0.6)');
Dm = [1 1 1 3 2 3; 2 2 2 1 1 1; 3 1 3 2 2 1; 1 2 2 3 1 3; 2 1 1 1 2 2; 3 2 3 2 1 2];
Cm = dd_compile(spM);
dm = dd_draws(Cm, 1, 300, 5, 'sobol', 0);
Mx = Cm.model(1);
Xm = dd_X(Mx, Dm);
th0 = [Mx.b0; Mx.s0(Mx.ridx)];
Irp = sum(dd_rowinfo(Mx, Xm, dm), 3);
Hrp = numHess(@(th) ellRP(Mx, Xm, dm.eta, th0, th), th0);
T = add(T, 'information', 'Cross-sectional mixed logit information (four distributions)', relErr(-Hrp, Irp) < 1e-4, sprintf('rel. error %.1e', relErr(-Hrp, Irp)));

Cp = Cm;
Cp.model(1).type = 'rppanel';
dp = dd_draws(Cp, 1, 300, 5, 'sobol', 4000);
Ipan = dd_panelinfo(Cp.model(1), Xm, dp, 1, 1e6);
Hpan = numHess(@(th) ellPanel(Mx, Xm, dp.eta, th0, th), th0);
T = add(T, 'information', 'Panel mixed logit information, exact over choice sequences', relErr(-Hpan, Ipan) < 1e-4, sprintf('rel. error %.1e', relErr(-Hpan, Ipan)));
[Isim, how] = dd_panelinfo(Cp.model(1), Xm, dp, 1, 1);
dEx = det(Ipan)^(-1 / Mx.Pt);
dSim = det(Isim)^(-1 / Mx.Pt);
T = add(T, 'information', 'Panel information from 4000 simulated respondents close to the exact one', ...
    strcmp(how, 'simulated') && abs(dSim / dEx - 1) < 0.08, sprintf('D-error %.4f vs %.4f', dSim, dEx));
% a panel of one task is the cross-sectional model
I1 = zeros(Mx.Pt);
for s = 1:size(Xm, 3)
    I1 = I1 + dd_panelinfo(Cp.model(1), Xm(:, :, s), dp, 1, 1e6);
end
T = add(T, 'information', 'Panel information with one task per respondent = cross-sectional information', relErr(I1, Irp) < 1e-10, sprintf('rel. error %.1e', relErr(I1, Irp)));
% the default integration (importance sampling around the mode of each choice
% sequence, dd_panelis) against the plain simulation with very many draws
dpL = dd_draws(Cp, 1, 2^17, 6, 'sobol', 0);
dRef = det(dd_panelinfo(Cp.model(1), Xm, dpL, 1, 1e6))^(-1 / Mx.Pt);
dpI = dd_draws(Cp, 1, 16, 5, 'sobol', 4000, 128);
[Iis, ~, dgi] = dd_panelinfo(Cp.model(1), Xm, dpI, 1, 1e6);
dIs = det(Iis)^(-1 / Mx.Pt);
T = add(T, 'information', 'Panel information by importance sampling (128 draws per sequence, four distributions) = plain simulation with 131072 draws', ...
    abs(dIs / dRef - 1) < 0.02 && dgi.noise < 0.05, sprintf('D-error %.4f vs %.4f', dIs, dRef));
dIsSim = det(dd_panelinfo(Cp.model(1), Xm, dpI, 1, 1))^(-1 / Mx.Pt);
T = add(T, 'information', 'The same from 4000 simulated respondents', abs(dIsSim / dRef - 1) < 0.06, sprintf('D-error %.4f vs %.4f', dIsSim, dRef));
% a long panel: 16 tasks per respondent
spL = spM;
spL.rows = 16;
spL.model.type = 'rppanel';
spL.model.random = struct('b_x', 'normal(0.6)', 'b_q', 'normal(0.8)', 'b_c', 'normal(0.5)', 'asc', 'normal(1)');
CL = dd_compile(spL);
rsL = RandStream('mt19937ar', 'Seed', 3);
DL = [randi(rsL, 3, 16, 1) randi(rsL, 2, 16, 1) randi(rsL, 3, 16, 1) randi(rsL, 3, 16, 1) randi(rsL, 2, 16, 1) randi(rsL, 3, 16, 1)];
ML = CL.model(1);
XL = dd_X(ML, DL);
dLong = @(d) det(dd_panelinfo(ML, XL, d, 1, 1))^(-1 / ML.Pt);
dL32 = dLong(dd_draws(CL, 1, 16, 5, 'sobol', 2000, 32));
dL512 = dLong(dd_draws(CL, 1, 16, 5, 'sobol', 2000, 512));
dLplain = dLong(dd_draws(CL, 1, 250, 5, 'sobol', 2000));
dLplainMany = dLong(dd_draws(CL, 1, 2^16, 5, 'sobol', 2000));
T = add(T, 'information', 'Long panel (16 tasks): importance sampling with 32 and with 512 draws per sequence agree', ...
    abs(dL32 / dL512 - 1) < 0.03, sprintf('D-error %.4f, %.4f', dL32, dL512));
T = add(T, 'information', 'Long panel: plain simulation understates the D-error with 250 draws and approaches importance sampling with 65536', ...
    dLplain < 0.96 * dL512 && abs(dLplainMany / dL512 - 1) < 0.02, sprintf('plain %.4f (250 draws), %.4f (65536 draws); importance sampling %.4f', dLplain, dLplainMany, dL512));
% criterion with excluded parameters
It = I;
Mex = Mm;
Mex.eff = [true(6, 1); false];
Mex.Peff = 6;
Om = inv(It);
T = add(T, 'information', 'D-error of a parameter subset (nuisance parameters)', ...
    abs(dd_crit(It, Mex, 'd') - det(Om(1:6, 1:6))^(1 / 6)) < 1e-12, '');
Mex.exMode = 'known';
T = add(T, 'information', 'D-error of a parameter subset (others known)', ...
    abs(dd_crit(It, Mex, 'd') - det(inv(It(1:6, 1:6)))^(1 / 6)) < 1e-12, '');
T = add(T, 'information', 'A-error = trace of the covariance / number of parameters', abs(dd_crit(It, Mm, 'a') - trace(Om) / 7) < 1e-12, '');
A = randn(RandStream('mt19937ar', 'Seed', 3), 6, 6, 40);
A = pagemtimes(A, 'transpose', A, 'none') + eye(6);
ld = dd_logdet(A);
ldRef = arrayfun(@(i) log(det(A(:, :, i))), (1:40)');
T = add(T, 'information', 'Batched log-determinant', max(abs(ld - ldRef)) < 1e-10 && dd_logdet(zeros(3)) == -Inf, '');

% Bayesian draws
spB = specH;
spB.model.priors.b_price = 'N(-0.02,0.005)';
spB.model.priors.b_pool = 'U(0.1,0.5)';
Cb = dd_compile(spB);
db1 = dd_draws(Cb, 4096, 1, 1, 'sobol', 0);
ip = find(strcmp(Cb.model(1).parNames, 'b_price'));
iq = find(strcmp(Cb.model(1).parNames, 'b_pool[1]'));
T = add(T, 'information', 'Prior draws: normal and uniform moments', ...
    abs(mean(db1.B(ip, :)) + 0.02) < 2e-5 && abs(std(db1.B(ip, :)) - 0.005) < 5e-5 && ...
    abs(mean(db1.B(iq, :)) - 0.3) < 1e-4 && min(db1.B(iq, :)) > 0.1 && max(db1.B(iq, :)) < 0.5, '');
db2 = dd_draws(Cb, 4096, 1, 1, 'sobol', 0);
T = add(T, 'information', 'Prior draws: same seed, same draws', isequal(db1.B, db2.B), '');
rB = dce_evaluate(spB, VH, 'Verbose', 0, 'EvalBayesDraws', 512);
frs = zeros(1, 512);
dbe = dd_draws(Cb, 512, 1, 2, 'sobol', 0);
for q = 1:512
    Iq = sum(dd_rowinfo(Cb.model(1), dd_X(Cb.model(1), D), struct('R', 1, 'B', dbe.B(:, q))), 3);
    frs(q) = det(Iq)^(-1 / 7);
end
T = add(T, 'information', 'Bayesian D-error = mean of the D-errors over the prior draws', abs(rB.models.Derror - mean(frs)) < 1e-12, sprintf('%.6f', rB.models.Derror));

% =========================== 3. Monte Carlo ===================================
if ~quick
    Nresp = 300;
    reps = 300;
    rsm = RandStream('mt19937ar', 'Seed', 11);
    bt = Mm.b0;
    bt(7) = 0.004;
    Itrue = sum(dd_rowinfo(Mm, X, struct('R', 1, 'B', bt)), 3);
    sePred = sqrt(diag(inv(Itrue)) / Nresp);
    est = zeros(7, reps);
    Xs = permute(X, [3 2 1]);
    dX = Xs(:, :, 1) - Xs(:, :, 2);
    p1 = 1 ./ (1 + exp(-dX * bt));
    for rr = 1:reps
        yv = rand(rsm, size(dX, 1), Nresp) < p1;
        b = zeros(7, 1);
        for it = 1:30
            pp = 1 ./ (1 + exp(-dX * b));
            g = dX' * sum(yv - pp, 2);
            Hs = Nresp * (dX' * (dX .* (pp .* (1 - pp))));
            step = Hs \ g;
            b = b + step;
            if max(abs(step)) < 1e-10, break; end
        end
        est(:, rr) = b;
    end
    ratio = std(est, 0, 2) ./ sePred;
    T = add(T, 'monte carlo', sprintf('MNL: s.d. of %d re-estimations vs predicted s.e. (N=%d)', reps, Nresp), ...
        all(ratio > 0.85 & ratio < 1.18), sprintf('ratios %.2f-%.2f', min(ratio), max(ratio)));
    T = add(T, 'monte carlo', 'MNL: estimates centred on the true parameters', ...
        all(abs(mean(est, 2) - bt) < 4 * sePred / sqrt(reps) + 0.03 * abs(bt) + 1e-4), '');
end

% =========================== 4. the search ====================================
% global optimum by enumeration: 4 tasks, 2 alternatives, three 2-level attributes
spE = struct();
spE.alternatives = {'a', 'b'};
spE.rows = 4;
spE.attributes.x1 = [0 1];
spE.attributes.x2 = [0 1];
spE.attributes.x3 = [0 1];
spE.dominance = false;
spE.unique = false;
u = 'b1*x1 + b2*x2 + b3*x3';
spE.model.utility = struct('a', u, 'b', u);
spE.model.priors = struct('b1', 0.6, 'b2', -0.4, 'b3', 0.9);
Ce = dd_compile(spE);
arr = unique(perms([1 1 2 2]), 'rows');
na = size(arr, 1);
best = Inf;
idx = cell(1, 6);
[idx{:}] = ndgrid(1:na);
allIdx = cell2mat(cellfun(@(v) v(:), idx, 'UniformOutput', false));
for i = 1:size(allIdx, 1)
    De = arr(allIdx(i, :), :)';
    Ie = sum(dd_rowinfo(Ce.model(1), dd_X(Ce.model(1), De), struct('R', 1, 'B', Ce.model(1).b0)), 3);
    f = dd_crit(Ie, Ce.model(1), 'd');
    if f < best, best = f; end
end
rE = dce_design(spE, 'Verbose', 0, 'Restarts', 20, 'Kicks', 2);
T = add(T, 'search', sprintf('Global optimum of a small problem (%d designs enumerated)', size(allIdx, 1)), ...
    abs(rE.models.Derror - best) < 1e-10, sprintf('%.6f vs %.6f', rE.models.Derror, best));

% hotel: constraints, frequencies, reproducibility, running objective
specS = specH;
specS.blocks = 3;
specS.reject = {'hotel1.price == 120 & hotel1.stars == 5', 'hotel2.price == 120 & hotel2.stars == 5'};
r1 = dce_design(specS, 'Verbose', 0, 'Restarts', 3, 'Kicks', 1);
r2 = dce_design(specS, 'Verbose', 0, 'Restarts', 3, 'Kicks', 1);
r3 = dce_design(specS, 'Verbose', 0, 'Restarts', 3, 'Kicks', 1, 'Seed', 7);
T = add(T, 'search', 'Same seed and specification give the same design', isequal(r1.design, r2.design) && r1.objective == r2.objective, '');
T = add(T, 'search', 'Another seed gives another design of similar quality', ~isequal(r1.design, r3.design) && abs(r3.objective / r1.objective - 1) < 0.03, ...
    sprintf('%.5f vs %.5f', r1.objective, r3.objective));
Vd = r1.design{:, 3:end};
T = add(T, 'search', 'Reject expressions hold in every task', ...
    ~any(Vd(:, 2) == 120 & Vd(:, 3) == 5) && ~any(Vd(:, 6) == 120 & Vd(:, 7) == 5) && r1.checks.constraintsOK, '');
T = add(T, 'search', 'Level balance is exact', all(all(r1.levelCounts(1:3, [2 3 4 6 7 8]) == 6)) && all(all(r1.levelCounts(1:2, [1 5 9]) == 9)), '');
T = add(T, 'search', 'No duplicate tasks, no dominated alternatives', isempty(r1.checks.duplicateTasks) && isempty(r1.checks.violations), '');
T = add(T, 'search', 'Blocks of equal size with every level in every block', ...
    isequal(accumarray(r1.blocks, 1)', [6 6 6]) && r1.blocking.absentAvoidable == 0, sprintf('Psi %.1f%%', r1.blocking.Psi));
re = dce_evaluate(specS, r1.design, 'Verbose', 0);
T = add(T, 'search', 'Evaluating the generated design returns the same D-error', abs(re.models.Derror - r1.models.Derror) < 1e-12, '');
T = add(T, 'search', 'Better than the Ngene design printed by Bliemer et al. (2026)', r1.models.Derror < 0.0462587 * 1.005, sprintf('%.5f vs 0.04626', r1.models.Derror));
Cs = dd_compile(specS);
os = dd_options(struct(), 'Restarts', 1, 'Kicks', 0);
DRs = dd_draws(Cs, 1, 1, 1, 'sobol', 0);
rs = RandStream('mrg32k3a', 'Seed', 1);
D0 = dd_initial(Cs, rs);
[Dl, fl] = dd_localsearch(Cs, D0, DRs, os, rs);
ffresh = dd_crit(sum(dd_rowinfo(Cs.model(1), dd_X(Cs.model(1), Dl), DRs), 3), Cs.model(1), 'd');
T = add(T, 'search', 'Running objective of the local search = fresh evaluation', abs(fl - ffresh) < 1e-12, sprintf('%.2e', abs(fl - ffresh)));
T = add(T, 'search', 'Random start satisfies constraints and frequencies', all(dd_rowok(Cs, D0)) && ~any(dd_dups(Cs, D0)), '');
% frequency bounds and level changes
spF = specH;
spF.frequency.price = struct('min', 4, 'max', 8);
rF = dce_design(spF, 'Verbose', 0, 'Restarts', 2, 'Kicks', 0);
cF = rF.levelCounts(1:3, [2 6]);
T = add(T, 'search', 'Frequency bounds are respected and used', all(cF(:) >= 4 & cF(:) <= 8) && rF.models.Derror <= r1.models.Derror * 1.02, mat2str(cF'));
% Bayesian search with screening equals a design that passes its own checks
rBs = dce_design(spB, 'Verbose', 0, 'Restarts', 1, 'Kicks', 0, 'BayesDraws', 64, 'EvalBayesDraws', 256);
rBf = dce_design(spB, 'Verbose', 0, 'Restarts', 1, 'Kicks', 0, 'BayesDraws', 64, 'EvalBayesDraws', 256, 'FullCheck', true);
T = add(T, 'search', 'Bayesian search: screening loses nothing against the full check', rBs.models.Derror <= rBf.models.Derror * 1.003, ...
    sprintf('%.5f vs %.5f', rBs.models.Derror, rBf.models.Derror));
% two models, normalised
spT = specH;
m2 = spT.model;
m2.name = 'B';
m2.priors.b_price = -0.04;
m1 = spT.model;
m1.name = 'A';
spT = rmfield(spT, 'model');
spT.models = {m1, m2};
rT = dce_design(spT, 'Verbose', 0, 'Restarts', 3, 'Kicks', 0);
T = add(T, 'search', 'Two models: both identified, normalised objective used', ...
    all([rT.models.identified]) && rT.search.normalized && all(rT.search.reference > 0), sprintf('%.4f %.4f', rT.models.Derror));
% multinomial logit + panel mixed logit, respondents answer one block
spP = struct();
spP.alternatives = {'a', 'b', 'none'};
spP.rows = 12;
spP.blocks = 2;
spP.attributes.x = [0 1 2];
spP.attributes.q = [0 1];
spP.attributes.cost = [1 2 3];
u = 'b_x.dummy*x + b_q*q + b_c*cost';
Up = struct('a', u, 'b', u, 'none', 'asc');
pm1 = struct('name', 'MNL', 'type', 'mnl', 'utility', Up, 'priors', struct('b_x', [0.3 0.6], 'b_q', 0.4, 'b_c', -0.5, 'asc', -0.5));
pm2 = struct('name', 'MXL', 'type', 'rppanel', 'utility', Up, 'priors', struct('b_x', [0.4 0.8], 'b_q', 0.5, 'b_c', -0.8, 'asc', -1), ...
    'random', struct('b_q', 'normal(0.5)', 'b_c', 'normal(0.4)', 'asc', 'normal(1)'));
spP.models = {pm1, pm2};
rP0 = dce_design(spP, 'Verbose', 0, 'Restarts', 2, 'Kicks', 0, 'RandomDraws', 64, 'EvalRandomDraws', 200, 'PanelPolish', 0);
rP = dce_design(spP, 'Verbose', 0, 'Restarts', 2, 'Kicks', 0, 'RandomDraws', 64, 'EvalRandomDraws', 200, 'PanelPolish', 3);
ps = rP.search.panelStage;
cntBlock = @(r) cell2mat(arrayfun(@(b) reshape(sort(r.design{r.blocks == b, 3:end}), 1, []), (1:2)', 'UniformOutput', false));
T = add(T, 'search', 'Panel stage: criterion evaluated per block and not worse than after the surrogate stage', ...
    strcmp(rP.models(2).method, 'exact over choice sequences') && all(ps.after <= ps.before * (1 + 1e-9) | ps.after(2) < ps.before(2)) && ...
    sum(ps.after ./ ps.before) <= 2 + 1e-9 && rP.checks.constraintsOK && isempty(rP.checks.duplicateTasks), ...
    sprintf('panel D-error %.4f -> %.4f, %d moves', ps.before(2), ps.after(2), ps.moves));
T = add(T, 'search', 'Panel stage keeps the level counts of every block', isequal(cntBlock(rP0), cntBlock(rP)) && isequal(rP0.blocks, rP.blocks), '');
spP1 = spP;
spP1.panelUnit = 'design';
rPd = dce_evaluate(spP1, rP.design, 'Verbose', 0, 'EvalPanelRespondents', 500);
T = add(T, 'search', 'Panel criterion for a respondent who answers every row (Ngene''s convention) is a different number', ...
    abs(rPd.models(2).Derror / rP.models(2).Derror - 1) > 0.005 && strcmp(rPd.models(2).method, '500 simulated respondents'), ...
    sprintf('%.4f vs %.4f by block', rPd.models(2).Derror, rP.models(2).Derror));

% status quo with fixed levels, a constant, and an excluded constant
spQ = struct();
spQ.alternatives = {'sq', 'a', 'b'};
spQ.rows = 12;
spQ.attributes.q = [0 1 2];
spQ.attributes.r = [0 1];
spQ.attributes.cost = [0 10 20 30];
spQ.altLevels.a.cost = [10 20 30];
spQ.altLevels.b.cost = [10 20 30];
spQ.fixed.sq = struct('q', 0, 'r', 0, 'cost', 0);
spQ.model.utility.sq = 'asc_sq + b_q.dummy*q + b_r*r + b_c*cost';
spQ.model.utility.a = 'b_q.dummy*q + b_r*r + b_c*cost';
spQ.model.utility.b = 'b_q.dummy*q + b_r*r + b_c*cost';
spQ.model.priors = struct('asc_sq', 0.3, 'b_q', [0.4 0.7], 'b_r', 0.3, 'b_c', -0.05);
spQ.reject = {'a.q == 0 & a.cost > sq.cost + 25'};
rQ = dce_design(spQ, 'Verbose', 0, 'Restarts', 2, 'Kicks', 0);
spQ2 = spQ;
spQ2.exclude = {'asc_sq'};
rQ2 = dce_evaluate(spQ2, rQ.design, 'Verbose', 0);
Omq = rQ.models.avc;
T = add(T, 'search', 'Status quo with fixed levels and a constant', rQ.models.identified && numel(rQ.columns) == 6 && rQ.checks.constraintsOK, '');
T = add(T, 'search', 'Excluding the constant gives the D-error of the other parameters', ...
    abs(rQ2.models.Derror - det(Omq(2:5, 2:5))^(1 / 4)) < 1e-10 && ~rQ2.models.inCriterion(1), '');
spQ3 = rmfield(spQ, 'attributes');
spQ3.attributes.q = [0 1 2];
spQ3.attributes.cost = [0 10 20 30];
spQ3.fixed.sq = struct('q', 0, 'cost', 0);
spQ3.model.utility = struct('sq', 'asc_sq + b_q.dummy*q + b_c*cost', 'a', 'b_q.dummy*q + b_c*cost', 'b', 'b_q.dummy*q + b_c*cost');
spQ3.model.priors = struct('asc_sq', 0.3, 'b_q', [0.4 0.7], 'b_c', -0.05);
T = add(T, 'errors', 'Fewer distinct feasible tasks than rows', throws(@() dce_design(spQ3, 'Verbose', 0, 'Restarts', 1), 'dce:search:infeasible'), '');

% =========================== 5. Ngene, export, errors ==========================
ng = ['design' newline ';alts = hotel1*, hotel2*' newline ';rows = 18' newline ';eff = (mnl,d)' newline ...
    ';cond:' newline 'if(hotel1.purpose = 0, hotel2.purpose = 0), if(hotel1.purpose = 1, hotel2.purpose = 1)' newline ...
    ';model:' newline 'U(hotel1) = b1[-0.02] * price[120,180,240] ? price' newline ...
    ' + b2.dummy[0.4|0.7] * stars[3,5,1] + b3.dummy[0.3|0.8] * meals[1,2,0] + b4.dummy[0.3] * pool[1,0]' newline ...
    ' + b5 * price * purpose[1,0]' newline '/' newline ...
    'U(hotel2) = b1 * price + b2.dummy * stars + b3.dummy * meals + b4.dummy * pool + b5 * price * purpose' newline '$'];
sN = dce_ngene(ng);
Vn = [VH(:, 1) VH(:, 2:5) VH(:, 1) VH(:, 6:9)];
Tn = array2table(Vn, 'VariableNames', {'hotel1.purpose', 'hotel1.price', 'hotel1.stars', 'hotel1.meals', 'hotel1.pool', ...
    'hotel2.purpose', 'hotel2.price', 'hotel2.stars', 'hotel2.meals', 'hotel2.pool'});
rN = dce_evaluate(sN, Tn, 'Verbose', 0);
rH = dce_evaluate(specH, VH, 'Verbose', 0);
T = add(T, 'ngene', 'Ngene script of Bliemer et al. (2026): same D-error as the native specification', ...
    abs(rN.models.Derror - rH.models.Derror) < 1e-12 && rN.checks.constraintsOK, sprintf('%.7f', rN.models.Derror));
T = add(T, 'ngene', 'Base level = last level listed; starred alternatives checked for dominance', ...
    sN.base.stars == 1 && sN.base.meals == 0 && isequal(sN.dominance, {'hotel1', 'hotel2'}) && numel(sN.reject) == 2, '');
ng2 = ['Design ;alts(M1) = a, b, sq ;alts(M2) = a, b, sq ;rows = 24 ;block = 2, minmax ;bdraws = sobol(300) ;rep = 250' ...
    ' ;eff = 2*M1(mnl,d,median) + 0.5*M2(rppanel,d) ;cond: if(a.x = [1,2], a.cost <> 10)' ...
    ' ;model(M1): U(a) = c1.dummy[(n,0.2,0.1)|(u,0.1,0.5)] * x[2,1,0](8,8,8) + c2[(n,-0.03,0.01)] * cost[10,20,30](6-10,6-10,6-10) /' ...
    ' U(b) = c1.dummy * x + c2 * cost / U(sq) = c0[0.5]' ...
    ' ;model(M2): U(a) = c1.dummy[n,0.2,0.3|n,0.4,0.3] * x[2,1,0] + c2[u,-0.05,-0.01] * cost[10,20,30] /' ...
    ' U(b) = c1.dummy * x + c2 * cost / U(sq) = c0[n,0.5,(u,0.5,1.5)] $'];
s3 = dce_ngene(ng2);
okD = isequal(s3.weights, [1 1]) && ~isfield(s3, 'exclude') && ~isfield(s3, 'constants') && ~isfield(s3, 'panelUnit') && ...
    ~(isfield(s3, 'search') && isfield(s3.search, 'Normalize')) && strcmp(s3.models{1}.stat, 'median');
T = add(T, 'ngene', 'Default translation uses the toolbox conventions (constants in, respondent = block, models scaled)', okD, '');
s4 = dce_ngene(strrep(ng2, 'M1(mnl,d,median)', 'M1(mnl,d)'));
s5 = dce_ngene(strrep(ng2, 'M1(mnl,d,median)', 'M1(mnl,d)'), 'Emulate', true);
T = add(T, 'ngene', 'Bayesian priors without mean/median: averaged by default, at the prior means when emulating Ngene', ...
    strcmp(s4.models{1}.stat, 'mean') && strcmp(s5.models{1}.stat, 'fixed'), '');
s6 = dce_ngene('design ;alts = a, b ;rows = 12 ;cond: if(a.x = 1, a.z = [0,1]) ;model: U(a) = b1 * x[0,1] + b2 * z[0,1,2] / U(b) = b1 * x + b2 * z $');
T = add(T, 'ngene', 'An attribute nested by ;cond gets frequency bounds instead of exact balance', ...
    isfield(s6, 'frequency') && isequal([s6.frequency.z.min s6.frequency.z.max], [2 6]) && ~isfield(s6.frequency, 'x'), '');
[s2, n2] = dce_ngene(ng2, 'Emulate', true);
okN = numel(s2.models) == 2 && isequal(s2.weights, [2 0.5]) && strcmp(s2.models{1}.stat, 'median') && ...
    strcmp(s2.models{2}.type, 'rppanel') && isequal(s2.frequency.x.exact, [8 8 8]) && isequal(s2.frequency.cost.min, [6 6 6]) && ...
    s2.blocks == 2 && isequal(s2.constants, false) && strcmp(s2.panelUnit, 'design') && ...
    strcmp(s2.reject{1}, '((a.x == 1 | a.x == 2)) & ~(a.cost ~= 10)') && abs(s2.models{2}.priors.c2 + 0.03) < 1e-12 && ...
    strcmp(s2.models{2}.random.c2, 'uniform(0.02)') && s2.search.EvalPanelRespondents == 2000 && s2.search.BayesDraws == 256;
T = add(T, 'ngene', 'Composite script: weights, statistics, priors, frequencies, constraints, constants', okN, sprintf('%d notes', numel(n2)));
C2 = dd_compile(s2);
T = add(T, 'ngene', 'Composite script compiles', C2.M == 2 && C2.model(2).Pr == 4 && C2.model(1).isBayes && all(~C2.model(1).eff(strcmp(C2.model(1).thBase, 'c0'))), '');
T = add(T, 'ngene', 'Unsupported Ngene features raise an error', ...
    throws(@() dce_ngene('design ;alts = a, b ;rows = 8 ;orth = ood ;model: U(a) = b1 * x[0,1] / U(b) = b1 * x $'), 'dce:ngene:unsupported') && ...
    throws(@() dce_ngene('design ;alts = a, b ;rows = 8 ;eff = (ecpanel,d) ;model: U(a) = b1 * x[0,1] / U(b) = b1 * x $'), 'dce:ngene:unsupported'), '');

% errors
T = add(T, 'errors', 'Unknown field', throws(@() dce_design(setfield(specH, 'row', 3), 'Verbose', 0), 'dce:spec:unknownField'), ''); %#ok<SFLD>
bad = specH;
bad.model.priors.b_stars = 0.4;
T = add(T, 'errors', 'Wrong number of priors for a dummy coded attribute', throws(@() dce_design(bad, 'Verbose', 0), 'dce:spec:priors'), '');
bad = specH;
bad.reject = {'hotel1.price > 0'};
T = add(T, 'errors', 'Constraints that cannot be met', throws(@() dce_design(bad, 'Verbose', 0, 'Restarts', 1), 'dce:search:infeasible'), '');
bad = specH;
bad.reject = {'hotel3.price > 0'};
T = add(T, 'errors', 'Constraint with an unknown name', throws(@() dce_design(bad, 'Verbose', 0), 'dce:spec:constraint'), '');
bad = specH;
bad.rows = 3;
T = add(T, 'errors', 'Too few rows for the parameters', throws(@() dce_design(bad, 'Verbose', 0), 'dce:spec:rows'), '');
bad = struct();
bad.alternatives = {'a', 'b'};
bad.rows = 12;
bad.attributes.purpose = [0 1];
bad.attributes.price = [1 2 3];
bad.scenario = {'purpose'};
bad.dominance = false;
bad.model.utility.a = 'b_price*price + b_main*purpose';
bad.model.utility.b = 'b_price*price + b_main*purpose';
T = add(T, 'errors', 'Parameter that the design cannot identify', throws(@() dce_design(bad, 'Verbose', 0, 'Restarts', 1, 'Kicks', 0), 'dce:search:notIdentified'), '');
T = add(T, 'errors', 'Unknown option', throws(@() dce_design(specH, 'Restart', 2), 'dce:options:unknown'), '');

% =========================== 6. methods, other design types, advice ============
% search method: chosen from the level-frequency rules
T = add(T, 'method', 'Exact level frequencies: exchanges within columns only', strcmp(r1.search.method, 'swap') && r1.search.rowMoves == 0, '');
rFs = dce_design(spF, 'Verbose', 0, 'Restarts', 3, 'Kicks', 0, 'Method', 'swap');
rFb = dce_design(spF, 'Verbose', 0, 'Restarts', 3, 'Kicks', 0);
T = add(T, 'method', 'Frequency bounds: row exchange is added and does not end worse than exchanges alone', ...
    strcmp(rFb.search.method, 'both') && rFb.search.candidatePool > 0 && rFb.models.Derror <= rFs.models.Derror * (1 + 1e-9), ...
    sprintf('%.5f vs %.5f, %d task replacements', rFb.models.Derror, rFs.models.Derror, rFb.search.rowMoves));
spN = specH;
spN.balance = 'none';
rNb = dce_design(spN, 'Verbose', 0, 'Restarts', 3, 'Kicks', 0);
T = add(T, 'method', 'No level balance: lower D-error than with exact balance, constraints kept', ...
    rNb.models.Derror < r1.models.Derror && rNb.checks.constraintsOK && isempty(rNb.checks.duplicateTasks), sprintf('%.5f', rNb.models.Derror));
% exhaustive search on a small problem, and the default search against it
spX = struct();
spX.alternatives = {'a', 'b'};
spX.rows = 5;
spX.attributes.x1 = [0 1];
spX.attributes.x2 = [0 1 2];
spX.dominance = false;
spX.model.utility = struct('a', 'b1*x1 + b2*x2', 'b', 'b1*x1 + b2*x2');
spX.model.priors = struct('b1', 0.5, 'b2', -0.4);
rXe = dce_design(spX, 'Verbose', 0, 'Method', 'exhaustive');
rXs = dce_design(spX, 'Verbose', 0);
T = add(T, 'method', sprintf('Exhaustive search (%d designs) and the default search reach the same optimum', rXe.search.exhaustive.designs), ...
    abs(rXe.objective - rXs.objective) < 1e-10 && rXe.checks.frequenciesOK, sprintf('%.6f vs %.6f', rXe.objective, rXs.objective));
T = add(T, 'method', 'Exhaustive search refuses a problem that is too large', throws(@() dce_design(specH, 'Verbose', 0, 'Method', 'exhaustive'), 'dce:search:exhaustive'), '');
% constants switch
spC = spQ;
spC.constants = false;
rC = dce_evaluate(spC, rQ.design, 'Verbose', 0);
T = add(T, 'method', 'spec.constants = false leaves the constants out of the criterion', ...
    abs(rC.models.Derror - rQ2.models.Derror) < 1e-12 && ~rC.models.inCriterion(1) && all(rQ.models.inCriterion), '');
% full factorial
rFull = dce_full(spE, 'Verbose', 0, 'Constraints', false);
spE2 = spE;
spE2.dominance = 'auto';
spE2.unique = true;
rFeas = dce_full(spE2, 'Verbose', 0);
bt = [0.6 -0.4 0.9];
[q1, q2, q3, q4, q5, q6] = ndgrid(0:1);
pa = [q1(:) q2(:) q3(:)];
pb = [q4(:) q5(:) q6(:)];
du = (pa - pb) .* bt;
keepTask = ~(all(du >= 0, 2) | all(du <= 0, 2));
T = add(T, 'designs', 'Full factorial: all 64 tasks; feasible and distinct ones counted independently', ...
    rFull.rows == 64 && rFeas.rows == nnz(keepTask) / 2 && rFeas.checks.constraintsOK, sprintf('%d feasible distinct', rFeas.rows));
rFr1 = dce_full(spE2, 'Verbose', 0, 'Rows', 6);
rFr2 = dce_design(setfield(setfield(spE2, 'type', 'full'), 'rows', 6), 'Verbose', 0); %#ok<SFLD>
T = add(T, 'designs', 'Random fraction of the full factorial is reproducible, also through spec.type', ...
    rFr1.rows == 6 && isequal(rFr1.design, rFr2.design) && all(ismember(rFr1.levels, rFeas.levels, 'rows')), '');
% orthogonal arrays
rso = RandStream('mt19937ar', 'Seed', 5);
oaCases = {8, 2 * ones(1, 7); 9, 3 * ones(1, 4); 16, 4 * ones(1, 5); 12, 2 * ones(1, 11); 18, [2 3 3 3 3 3 3 3]; 36, [3 3 3 3 2 2]};
okOA = true;
txt = '';
for i = 1:size(oaCases, 1)
    [Ao, isOA] = dd_oa(oaCases{i, 1}, oaCases{i, 2}, 30, rso);
    okOA = okOA && isOA && isequal(size(Ao), [oaCases{i, 1} numel(oaCases{i, 2})]);
    txt = [txt sprintf('L%d:%d ', oaCases{i, 1}, isOA)]; %#ok<AGROW>
end
T = add(T, 'designs', 'Orthogonal arrays found: L8, L9, L16, L12 (2^11), L18 (2 x 3^7), L36 (3^4 x 2^2)', okOA, strtrim(txt));
spO = struct();
spO.alternatives = {'a', 'b'};
spO.rows = 12;
spO.attributes.x1 = [0 1 2];
spO.attributes.x2 = [1 2 3];
spO.attributes.x3 = [0 1];
spO.dominance = false;
spO.model.utility = struct('a', 'b1*x1 + b2*x2 + b3*x3', 'b', 'b1*x1 + b2*x2 + b3*x3');
spO.model.priors = struct('b1', 0.3, 'b2', -0.2, 'b3', 0.5);
rO = dce_orthogonal(spO, 'Verbose', 0);
crO = abs(rO.correlations - eye(6));
T = add(T, 'designs', 'Simultaneous orthogonal design: all columns uncorrelated and balanced, rows raised to 36', ...
    rO.search.orthogonal && rO.rows == 36 && max(crO(:)) < 1e-12 && all(all(rO.levelCounts(1:3, [1 2 4 5]) == 12)) && rO.models.identified, ...
    sprintf('%d rows, D-error %.4f', rO.rows, rO.models.Derror));
rOs = dce_orthogonal(spO, 'Verbose', 0, 'Type', 'sequential');
crS = abs(rOs.correlations);
T = add(T, 'designs', 'Sequential orthogonal design: uncorrelated within each alternative, 18 rows', ...
    rOs.search.orthogonal && rOs.rows == 18 && max(max(crS(1:3, 1:3) - eye(3))) < 1e-12 && max(max(crS(4:6, 4:6) - eye(3))) < 1e-12 && rOs.models.identified, ...
    sprintf('D-error %.4f', rOs.models.Derror));
spOb = spO;
spOb.blocks = 2;
spOb.type = 'orthogonal';
rOb = dce_design(spOb, 'Verbose', 0);
T = add(T, 'designs', 'Orthogonal design with blocks as a column of the array: every block perfectly balanced', ...
    rOb.search.blockInArray && abs(rOb.blocking.Psi - 100) < 1e-9 && isequal(accumarray(rOb.blocks, 1)', [18 18]), '');
rEff = dce_design(setfield(spO, 'rows', 36), 'Verbose', 0, 'Restarts', 3, 'Kicks', 0); %#ok<SFLD>
T = add(T, 'designs', 'Efficient design of the same size has a lower D-error than the orthogonal one', rEff.models.Derror < rO.models.Derror, ...
    sprintf('%.4f vs %.4f', rEff.models.Derror, rO.models.Derror));
sF = dce_ngene('design ;alts = a, b ;rows = all ;fact ;model: U(a) = b1 * x[0,1] + b2 * z[0,1,2] / U(b) = b1 * x + b2 * z $');
sO = dce_ngene('design ;alts = a, b ;rows = 18 ;orth = seq ;model: U(a) = b1 * x[0,1] + b2 * z[0,1,2] / U(b) = b1 * x + b2 * z $');
T = add(T, 'designs', 'Ngene ;fact and ;orth are translated', strcmp(sF.type, 'full') && ~isfield(sF, 'rows') && strcmp(sO.type, 'orthogonal-sequential'), '');
% advice
aH = dce_advise(specS);
a20 = dce_advise(setfield(specH, 'rows', 20)); %#ok<SFLD>
a9 = dce_advise(setfield(setfield(specH, 'rows', 18), 'blocks', 9)); %#ok<SFLD>
T = add(T, 'advice', 'Size of the problem and the time an exhaustive search would take', ...
    any(strcmp(aH.topic, 'size')) && contains(aH.text{find(strcmp(aH.topic, 'size'), 1)}, 'possible designs'), '');
T = add(T, 'advice', 'Rows that cannot balance a 3-level attribute: nearest numbers that can', ...
    any(strcmp(a20.topic, 'rows') & contains(a20.text, 'multiple of 6') & contains(a20.text, '18 or 24')), '');
T = add(T, 'advice', 'Blocks of two tasks: every level cannot appear, fewer blocks suggested', ...
    any(strcmp(a9.topic, 'blocks') & strcmp(a9.level, 'warning')) && any(strcmp(a9.topic, 'blocks') & contains(a9.text, 'Fewer blocks')), '');
aR = struct2table(r1.advice);
ie = find(strcmp(aR.topic, 'efficiency'), 1);
tk = regexp(aR.text{ie}, 'at least ([\d\.]+)%', 'tokens', 'once');
T = add(T, 'advice', 'Lower bound on the D-error of any design, and the design''s efficiency against it', ...
    ~isempty(ie) && str2double(tk{1}) > 80 && str2double(tk{1}) <= 100, sprintf('at least %s%% efficient', tk{1}));
T = add(T, 'advice', 'Sample size, search method and agreement of the starts are reported', ...
    any(strcmp(aR.topic, 'parameters')) && any(strcmp(aR.topic, 'search') & contains(aR.text, 'Search method')), '');
spZ = specH;
spZ.model = rmfield(spZ.model, 'priors');
spZ.dominance = false;
aZ = dce_advise(spZ);
T = add(T, 'advice', 'All priors zero: warning', any(strcmp(aZ.topic, 'priors') & strcmp(aZ.level, 'warning')), '');

% export and command-line entry
tmp = tempname;
mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
f = dce_export(r1, tmp, 'Name', 'hotel_test');
js = jsondecode(fileread(f.json));
Tl = readtable(f.long);
T = add(T, 'export', 'Files written and readable', all(cellfun(@isfile, struct2cell(f))) && js.rows == 18 && height(Tl) == 36, '');
T = add(T, 'export', 'Existing files are not overwritten', throws(@() dce_export(r1, tmp, 'Name', 'hotel_test'), 'dce:export:exists'), '');
sj = fullfile(tmp, 'spec.json');
fid = fopen(sj, 'w');
fprintf(fid, '%s', jsonencode(specS));
fclose(fid);
txt = evalc('sm = dce_cli(sj, fullfile(tmp, ''cli''), ''Restarts'', 3, ''Kicks'', 1);');
T = add(T, 'export', 'JSON specification through dce_cli gives the same design', ...
    sm.ok && contains(txt, 'DCE_RESULT {') && abs(sm.objective - r1.objective) < 1e-12 && isfile(sm.files.xlsx), '');
clear cleanup

% ---- report ------------------------------------------------------------------------
report = struct2table(T);
nFail = nnz(~report.pass);
fprintf('\n%d checks, %d failed, %.0f s\n', height(report), nFail, toc(t0));
if nFail > 0
    error('dce:selftest:failed', '%d check(s) failed.', nFail);
end
if nargout == 0
    clear report
end
end

% =============================================================================
function T = add(T, group, check, pass, detail)
T(end + 1) = struct('group', group, 'check', check, 'pass', logical(pass), 'detail', detail);
fprintf('%s  %-12s %s%s\n', ternary(pass, 'PASS', 'FAIL'), group, check, ternary(isempty(detail), '', ['  [' detail ']']));
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end

function tf = throws(f, id)
tf = false;
try
    f();
catch err
    tf = strcmp(err.identifier, id);
    if ~tf
        fprintf('      (expected %s, got %s: %s)\n', id, err.identifier, err.message);
    end
end
end

function e = relErr(A, B)
e = norm(A - B, 'fro') / norm(B, 'fro');
end

function H = numHess(f, x)
n = numel(x);
H = zeros(n);
h = 1e-4 * max(abs(x), 0.1);
for a = 1:n
    for b = a:n
        ea = zeros(n, 1); ea(a) = h(a);
        eb = zeros(n, 1); eb(b) = h(b);
        H(a, b) = (f(x + ea + eb) - f(x + ea - eb) - f(x - ea + eb) + f(x - ea - eb)) / (4 * h(a) * h(b));
        H(b, a) = H(a, b);
    end
end
end

function v = ellMNL(X, b0, b)
p0 = softmaxP(X, b0);
p = softmaxP(X, b);
v = sum(p0(:) .* log(p(:)));
end

function p = softmaxP(X, B)
V = pagemtimes(X, reshape(B, size(X, 2), 1, 1, []));
V = V - max(V, [], 1);
p = exp(V);
p = p ./ sum(p, 1);
end

function v = ellRP(Mm, X, eta, th0, th)
P0 = mean(softmaxP(X, betas(Mm, th0, eta)), 4);
P = mean(softmaxP(X, betas(Mm, th, eta)), 4);
v = sum(P0(:) .* log(P(:)));
end

function v = ellPanel(Mm, X, eta, th0, th)
L0 = seqProb(X, betas(Mm, th0, eta));
L = seqProb(X, betas(Mm, th, eta));
v = sum(L0 .* log(L));
end

function B = betas(Mm, th, eta)
B = dd_betas(Mm, th(1:Mm.P), th(Mm.P + 1:end), eta);
end

function L = seqProb(X, B)
[J, ~, T] = size(X);
p = softmaxP(X, B);
lp = reshape(log(p), J, T, []);
n = J^T;
idx = (0:n - 1)';
LP = zeros(n, size(lp, 3));
for t = 1:T
    y = mod(floor(idx / J^(t - 1)), J) + 1;
    LP = LP + reshape(lp(y, t, :), n, []);
end
L = mean(exp(LP), 2);
end

function [spec, V] = hotel()
% Hotel example and the Ngene design of Bliemer, Hess, Collins & Amaris (2026), Table 2.
spec = struct();
spec.name = 'hotel';
spec.alternatives = {'hotel1', 'hotel2'};
spec.rows = 18;
spec.attributes.purpose = [0 1];
spec.attributes.price = [240 180 120];
spec.attributes.stars = [1 3 5];
spec.attributes.meals = [0 1 2];
spec.attributes.pool = [0 1];
spec.scenario = {'purpose'};
u = 'b_price*price + b_stars.dummy*stars + b_meals.dummy*meals + b_pool.dummy*pool + b_pp*price*purpose';
spec.model.utility.hotel1 = u;
spec.model.utility.hotel2 = u;
spec.model.priors.b_price = -0.02;
spec.model.priors.b_stars = [0.4 0.7];
spec.model.priors.b_meals = [0.3 0.8];
spec.model.priors.b_pool = 0.3;
T2 = [1 0 1 2 0 1 2 1 1; 0 1 1 2 0 0 2 1 1; 1 1 1 2 1 2 0 0 0; 1 0 2 0 1 1 0 2 0; 1 0 0 1 1 1 1 0 0; 1 0 1 1 1 2 0 0 0
    0 2 0 0 1 1 1 1 0; 0 2 0 0 1 2 2 1 0; 0 1 2 0 0 0 1 2 1; 1 2 1 1 1 1 2 2 0; 1 2 0 0 0 0 2 2 1; 0 2 1 0 0 0 2 2 1
    0 1 2 1 0 0 0 2 1; 0 1 2 2 1 2 0 1 0; 1 0 2 2 0 1 1 1 1; 0 0 2 2 0 2 0 0 1; 0 1 0 1 1 0 1 0 0; 1 2 0 1 0 2 1 0 1];
price = [240 180 120];
stars = [1 3 5];
V = T2;
V(:, [2 6]) = price(T2(:, [2 6]) + 1);
V(:, [3 7]) = stars(T2(:, [3 7]) + 1);
end
