function [It, how, pdg] = dd_designinfo(C, m, D, y, d, o, rsub)
%DD_DESIGNINFO Fisher information of a whole design for one model. Internal.
%   It is Pt-by-Pt-by-1-by-R (R prior draws): the information from one
%   replication of the design, i.e. S choice observations.
%   mnl, rp:  sum of the task informations.
%   rppanel:  sum over blocks of the information of one respondent who answers
%             that block (C.panelUnit 'block', y gives the blocks), or the
%             information of one respondent who answers every row ('design').
%   rsub (optional) restricts the prior draws that are evaluated.
%   pdg (rppanel only, otherwise []): how the expectation over the coefficients
%   was taken: .integration 'importance' or 'prior', .tasks (most tasks per
%   respondent), and for importance sampling .ess, .noise and .pooled (see
%   dd_panelis; mean, largest and any over blocks and prior draws).

Mm = C.model(m);
how = 'analytical';
pdg = [];
if nargin < 7 || isempty(rsub)
    rsub = 1:d.R;
end
X = dd_X(Mm, D);
if ~strcmp(Mm.type, 'rppanel')
    if numel(rsub) < d.R
        d = subDraws(d, rsub);
    end
    % rows in chunks, so that large designs with many prior draws fit in memory
    n = size(X, 3);
    chunk = max(1, floor(2e7 / (d.R * Mm.Pt^2)));
    It = zeros(Mm.Pt, Mm.Pt, 1, d.R);
    for c0 = 1:chunk:n
        It = It + sum(dd_rowinfo(Mm, X(:, :, c0:min(n, c0 + chunk - 1)), d), 3);
    end
    if strcmp(Mm.type, 'rp')
        how = 'simulated probabilities';
    end
    return
end
if strcmp(C.panelUnit, 'block') && ~isempty(y)
    groups = arrayfun(@(b) find(y == b), unique(y(:))', 'UniformOutput', false);
else
    groups = {(1:size(D, 1))'};
end
It = zeros(Mm.Pt, Mm.Pt, 1, numel(rsub));
used = {};
ess = [];
noise = [];
pooled = false;
Ipool = zeros(Mm.Pt, Mm.Pt, numel(rsub));
for i = 1:numel(rsub)
    for g = 1:numel(groups)
        [Ib, meth, dg] = dd_panelinfo(Mm, X(:, :, groups{g}), d, rsub(i), o.PanelExactMax);
        It(:, :, 1, i) = It(:, :, 1, i) + Ib;
        used{end + 1} = meth; %#ok<AGROW>
        if ~isempty(dg)
            ess(end + 1) = dg.ess; %#ok<AGROW>
            noise(end + 1) = dg.noise; %#ok<AGROW>
            Ipool(:, :, i) = Ipool(:, :, i) + dg.pooled;
        end
    end
    % The noise-corrected estimate is a difference; if the sum over blocks is
    % not positive definite although the pooled estimate is, use the latter.
    if ~isempty(ess) && ~isfinite(dd_logdet(It(:, :, 1, i))) && isfinite(dd_logdet(Ipool(:, :, i)))
        It(:, :, 1, i) = Ipool(:, :, i);
        pooled = true;
    end
end
pdg = struct('integration', 'prior', 'tasks', max(cellfun(@numel, groups)), 'ess', NaN, 'noise', NaN, 'pooled', false);
if ~isempty(ess)
    pdg.integration = 'importance';
    pdg.ess = mean(ess);
    pdg.noise = max(noise);
    pdg.pooled = pooled;
end
if all(strcmp(used, 'exact'))
    how = 'exact over choice sequences';
elseif all(strcmp(used, 'simulated'))
    how = sprintf('%d simulated respondents', size(d.etaSim, 2));
else
    how = sprintf('exact or %d simulated respondents, by block', size(d.etaSim, 2));
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
