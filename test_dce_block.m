% TEST_DCE_BLOCK  Checks dce_block against Bliemer, Hess, Collins & Amaris
% (2026), ITLS-WP-26-22, and against exhaustive enumeration on small designs.
%
% Run from this folder:  test_dce_block
% or:  matlab -batch "cd('<this folder>'); test_dce_block"

clear; clc;
nFail = 0;
tStart = tic;

%% Hotel design of the paper (Table 2, design coding) --------------------------
% Columns: task, blocking N (naive), O (Ngene near-orthogonal), B1 (balanced,
% within), B2 (balanced, across), Purpose, Hotel 1 Price Stars Meals Pool,
% Hotel 2 Price Stars Meals Pool.
T2 = [
 1 1 3 3 3 1 0 1 2 0 1 2 1 1
 2 1 1 1 2 0 1 1 2 0 0 2 1 1
 3 1 2 2 2 1 1 1 2 1 2 0 0 0
 4 1 2 1 1 1 0 2 0 1 1 0 2 0
 5 1 2 1 3 1 0 0 1 1 1 1 0 0
 6 1 1 2 1 1 0 1 1 1 2 0 0 0
 7 2 1 2 1 0 2 0 0 1 1 1 1 0
 8 2 3 3 2 0 2 0 0 1 2 2 1 0
 9 2 3 3 3 0 1 2 0 0 0 1 2 1
10 2 2 3 1 1 2 1 1 1 1 2 2 0
11 2 1 2 2 1 2 0 0 0 0 2 2 1
12 2 1 1 1 0 2 1 0 0 0 2 2 1
13 3 3 2 1 0 1 2 1 0 0 0 2 1
14 3 2 1 3 0 1 2 2 1 2 0 1 0
15 3 2 2 2 1 0 2 2 0 1 1 1 1
16 3 1 3 3 0 0 2 2 0 2 0 0 1
17 3 3 3 2 0 1 0 1 1 0 1 0 0
18 3 3 1 3 1 2 0 1 0 2 1 0 1];
blk = T2(:, 2:5);
X = T2(:, 6:14);
names = {'Purpose', 'Price1', 'Stars1', 'Meals1', 'Pool1', 'Price2', 'Stars2', 'Meals2', 'Pool2'};
across = [1 2 3 4 5 2 3 4 5];

%% 1. The balance measure reproduces the paper's figures -------------------------
fprintf('\n== 1. Balance of the paper''s blocking columns ==\n');
labelsBlk = {'N (naive)', 'O (Ngene)', 'B1 (balanced within)', 'B2 (balanced across)'};
expW = [55 67 88.89 70];          % Section 6.1 and Table 9
expA = [78 79 NaN 100];
for j = 1:4
    [~, oW] = dce_block(X, 3, 'Evaluate', blk(:, j), 'Names', names, 'Verbose', false);
    [~, oA] = dce_block(X, 3, 'Evaluate', blk(:, j), 'Names', names, 'Groups', across, 'Verbose', false);
    fprintf('  %-22s within %6.2f%% (paper %g)   across %6.2f%% (paper %g)\n', ...
        labelsBlk{j}, oW.Psi, expW(j), oA.Psi, expA(j));
    nFail = nFail + check(abs(oW.Psi - expW(j)) < 0.5, sprintf('%s within-alternative Psi', labelsBlk{j}));
    if ~isnan(expA(j))
        nFail = nFail + check(abs(oA.Psi - expA(j)) < 0.5, sprintf('%s across-alternative Psi', labelsBlk{j}));
    end
end
[~, oB1] = dce_block(X, 3, 'Evaluate', blk(:, 3), 'Verbose', false);
nFail = nFail + check(abs(oB1.h^2 - 5/9) < 1e-12, 'B1: h^2 = 0.5556 (paper''s global optimum)');

% Ngene criterion (eq. 5): naive 0.221, Ngene's O column 0
[~, oN] = dce_block(X, 3, 'Evaluate', blk(:, 1), 'Verbose', false);
[~, oO] = dce_block(X, 3, 'Evaluate', blk(:, 2), 'Verbose', false);
fprintf('  Ngene criterion: naive %.3f (paper 0.221), O %.3g (paper 0), B1 %.3f (paper 0.071)\n', ...
    oN.corrNgene, oO.corrNgene, oB1.corrNgene);
nFail = nFail + check(abs(oN.corrNgene - 0.221) < 5e-4, 'Ngene criterion of naive blocking');
nFail = nFail + check(oO.corrNgene < 1e-12, 'Ngene criterion of Ngene''s blocking is zero');

% Relabelling blocks changes Ngene's criterion (paper: 0.221 to 0.424) but not Psi
pm = perms(1:3);
cr = zeros(6, 1);
ps = zeros(6, 1);
for j = 1:6
    yy = pm(j, blk(:, 1))';
    [~, oo] = dce_block(X, 3, 'Evaluate', yy, 'Verbose', false);
    cr(j) = oo.corrNgene;
    ps(j) = oo.Psi;
end
fprintf('  Six relabellings of naive blocking: corr %.3f-%.3f, Psi %.2f-%.2f\n', min(cr), max(cr), min(ps), max(ps));
nFail = nFail + check(abs(min(cr) - 0.221) < 5e-4 && abs(max(cr) - 0.424) < 5e-4, 'Ngene criterion range 0.221-0.424');
nFail = nFail + check(max(ps) - min(ps) < 1e-10, 'Psi invariant to block labels');

% Absent levels flagged as in the paper: naive blocking excludes a level from
% at least one block in five of nine columns; Ngene's O never shows Price2 = 0 in block 2
nCols = numel(unique(oN.absent.Column));
fprintf('  Naive blocking: %d columns miss a level in some block (paper: 5)\n', nCols);
nFail = nFail + check(nCols == 5, 'naive blocking: five columns with an absent level');
nFail = nFail + check(any(oO.absent.Block == 2 & strcmp(oO.absent.Column, 'X6') & oO.absent.Level == 0), ...
    'Ngene blocking: Price (Hotel 2) level 0 absent from block 2');

%% 2. Search finds the paper's optima --------------------------------------------
fprintf('\n== 2. Search on the hotel design ==\n');
[yW, oW] = dce_block(X, 3, 'Names', names);
nFail = nFail + check(oW.Psi >= 88.888, 'within-alternative search reaches 88.89% (global optimum)');
nFail = nFail + check(isequal(sort(accumarray(yW, 1))', [6 6 6]), 'equal block sizes');
[yA, oA] = dce_block(X, 3, 'Names', names, 'Groups', across);
nFail = nFail + check(oA.Psi > 99.999, 'across-alternative search reaches 100%');
fprintf('  Across-alternative optimum evaluated within alternatives: %.2f%% (paper''s B2: 70%%)\n', oA.PsiWithin);
[yW2, ~] = dce_block(X, 3, 'Names', names, 'Verbose', false);
nFail = nFail + check(isequal(yW, yW2), 'same seed gives the same blocking');

%% 3. Global optimum on small designs (exhaustive enumeration) ---------------------
fprintf('\n== 3. Exhaustive check on random level-balanced designs ==\n');
cfgs = {12, 3, [2 3 3 4 2 3]; 12, 2, [2 3 4 3 2 4 3]; 12, 4, [3 2 3 2 3]};
for c = 1:size(cfgs, 1)
    [S, B, Lk] = cfgs{c, :};
    Pall = allPartitions(S, B);
    Sb = (S / B) * ones(1, B);
    for seed = 1:3
        rng(100 * c + seed, 'twister');
        Xr = zeros(S, numel(Lk));
        for k = 1:numel(Lk)
            Xr(:, k) = repmat((1:Lk(k))', S / Lk(k), 1);
            Xr(:, k) = Xr(randperm(S), k);
        end
        h2all = zeros(size(Pall, 2), 1);
        mall = h2all;
        for j = 1:size(Pall, 2)
            [h2all(j), mall(j)] = refObjective(Xr, Pall(:, j), B, Sb);
        end
        hmax = sqrt(B * sum(Lk - 1));
        psiOpt = 100 * (1 - sqrt(min(h2all)) / hmax);
        mOpt = min(mall);
        psiLex = 100 * (1 - sqrt(min(h2all(mall == mOpt))) / hmax);
        [~, o1] = dce_block(Xr, B, 'RequireAll', false, 'Restarts', 50, 'Iterations', 500, 'Verbose', false);
        [y2, o2] = dce_block(Xr, B, 'Restarts', 50, 'Iterations', 500, 'Verbose', false);
        [h2c, mc] = refObjective(Xr, y2, B, Sb);
        fprintf('  S=%d B=%d design %d: %d partitions; optimum %.3f%%, search %.3f%%; with absent-level priority %.3f%%/%d absent, search %.3f%%/%d\n', ...
            S, B, seed, size(Pall, 2), psiOpt, o1.Psi, psiLex, mOpt, o2.Psi, mc);
        nFail = nFail + check(abs(o1.Psi - psiOpt) < 1e-9, 'search reaches global optimum of h^2');
        nFail = nFail + check(mc == mOpt && abs(o2.Psi - psiLex) < 1e-9, 'search reaches lexicographic optimum');
        nFail = nFail + check(abs(100 * (1 - sqrt(h2c) / hmax) - o2.Psi) < 1e-9, 'objective matches independent implementation');
    end
end

%% 4. MNL information diagnostics (Ngene priors of Appendix A) ----------------------
fprintf('\n== 4. MNL information ==\n');
% Parameters: price, stars ***, stars *****, breakfast, all-inclusive, pool, price x business
priceVal = [240 180 120];
Xi = zeros(18, 2, 7);
for j = 1:2
    off = 1 + 4 * (j - 1);
    pr = reshape(priceVal(X(:, off + 1) + 1), [], 1);
    st = X(:, off + 2);
    me = X(:, off + 3);
    po = X(:, off + 4);
    Xi(:, j, :) = reshape([pr, st == 1, st == 2, me == 1, me == 2, po == 1, pr .* (X(:, 1) == 1)], 18, 1, 7);
end
info = struct('X', Xi, 'beta', [-0.02 0.4 0.7 0.3 0.8 0.3 0]');
[~, oI] = dce_block(X, 3, 'Evaluate', blk(:, 3), 'Info', info, 'Verbose', false);
% numerical check: Fisher information = minus Hessian of the expected log-likelihood
b0 = info.beta;
Hn = zeros(7);
e = 1e-4;
for a = 1:7
    for b = 1:7
        ea = zeros(7, 1); ea(a) = e;
        eb = zeros(7, 1); eb(b) = e;
        Hn(a, b) = (expLL(Xi, b0, b0 + ea + eb) - expLL(Xi, b0, b0 + ea - eb) ...
            - expLL(Xi, b0, b0 - ea + eb) + expLL(Xi, b0, b0 - ea - eb)) / (4 * e^2);
    end
end
relErr = norm(-Hn - oI.info.I) / norm(oI.info.I);
fprintf('  D-error %.4g; shares %s; block ranks %s (P = 7)\n', oI.info.Derror, ...
    sprintf('%.3f ', oI.info.share), sprintf('%d ', oI.info.blockRank));
fprintf('  Information matrix vs numerical Hessian: relative error %.1e\n', relErr);
nFail = nFail + check(relErr < 1e-5, 'MNL information matrix matches numerical Hessian');
nFail = nFail + check(abs(sum(oI.info.share) - 1) < 1e-10, 'information shares sum to one');
nFail = nFail + check(all(isinf(oI.info.DerrorBlock)) && all(isfinite(oI.info.DerrorLOBO)), ...
    'single blocks are rank deficient (6 tasks, 7 parameters); leave-one-out D-errors are finite');

% information-based choice among near-best blockings
[~, oSel0] = dce_block(X, 3, 'Info', info, 'Tolerance', 0, 'Verbose', false);
[~, oSel2] = dce_block(X, 3, 'Info', info, 'Tolerance', 2);
[~, oBal] = dce_block(X, 3, 'Info', info, 'Select', 'balance', 'Verbose', false);
fprintf('  Worst leave-one-block-out D-error ratio: balance only %.4f, info (tol 0) %.4f, info (tol 2) %.4f at Psi %.2f%%\n', ...
    max(oBal.info.LOBOratio), max(oSel0.info.LOBOratio), max(oSel2.info.LOBOratio), oSel2.Psi);
nFail = nFail + check(max(oSel0.info.LOBOratio) <= max(oBal.info.LOBOratio) + 1e-12 && oSel0.Psi >= oBal.Psi - 1e-9, ...
    'info selection at tolerance 0 keeps the best balance and does not worsen information');
nFail = nFail + check(oSel2.Psi >= oW.Psi - 2 - 1e-9 && max(oSel2.info.LOBOratio) <= max(oSel0.info.LOBOratio) + 1e-12, ...
    'info selection at tolerance 2 stays within 2 points and improves or keeps information');

%% 5. Options: extremes, unequal blocks, tables -------------------------------------
fprintf('\n== 5. Options ==\n');
[yE, oE] = dce_block(X, 3, 'Names', names, 'Extremes', [2 6], 'Verbose', false);
okE = true;
for b = 1:3
    for k = [2 6]
        z = X(yE == b, k);
        okE = okE && any(z == 0) && any(z == 2);
    end
end
nFail = nFail + check(okE && all(oE.extremes.HasLowest & oE.extremes.HasHighest), 'both price extremes in every block');
[yU, oU] = dce_block(X, 3, 'BlockSizes', [7 6 5], 'Verbose', false);
nFail = nFail + check(isequal(accumarray(yU, 1)', [7 6 5]), 'unequal block sizes respected');
Xt = array2table(X, 'VariableNames', names);
Xt.Purpose = categorical(Xt.Purpose, [0 1], {'leisure', 'business'});
[yT, oT] = dce_block(Xt, 3, 'Verbose', false);
nFail = nFail + check(abs(oT.Psi - oW.Psi) < 1e-9, 'table input with a categorical column');
fprintf('  Extremes: Psi %.2f%%; unequal sizes [7 6 5]: Psi %.2f%%\n', oE.Psi, oU.Psi);

%% Summary ------------------------------------------------------------------------
fprintf('\n%d check(s) failed; %.1f s\n', nFail, toc(tStart));
if nFail > 0
    error('test_dce_block:failed', '%d check(s) failed.', nFail);
end

% ===================================================================================
function f = check(ok, msg)
if ok
    fprintf('  PASS  %s\n', msg);
    f = 0;
else
    fprintf('  FAIL  %s\n', msg);
    f = 1;
end
end

function P = allPartitions(S, B)
% All partitions of tasks 1..S into B unlabelled blocks of equal size.
P = enumRec(zeros(S, 1), (1:S)', 1, S / B);
end

function P = enumRec(y, rest, b, m)
if isempty(rest)
    P = y;
    return
end
first = rest(1);
others = rest(2:end);
combs = nchoosek(others, m - 1);
P = [];
for i = 1:size(combs, 1)
    members = [first; combs(i, :)'];
    yy = y;
    yy(members) = b;
    P = [P, enumRec(yy, setdiff(rest, members), b + 1, m)]; %#ok<AGROW>
end
end

function [h2, miss] = refObjective(Xc, y, B, Sb)
% Independent implementation of eqs. (9)-(16): squared CV of relative counts.
h2 = 0;
miss = 0;
for k = 1:size(Xc, 2)
    L = max(Xc(:, k));
    for b = 1:B
        n = accumarray(Xc(y == b, k), 1, [L 1])';
        Rr = L * n / Sb(b);
        h2 = h2 + mean((Rr - 1).^2);
        if Sb(b) >= L
            miss = miss + nnz(n == 0);
        end
    end
end
end

function v = expLL(Xi, b0, b)
% Expected log-likelihood at b when choices follow the MNL at b0.
[S, J, P] = size(Xi);
v = 0;
for s = 1:S
    Xs = reshape(Xi(s, :, :), J, P);
    p0 = exp(Xs * b0); p0 = p0 / sum(p0);
    u = Xs * b;
    lp = u - max(u) - log(sum(exp(u - max(u))));
    v = v + p0' * lp;
end
end
