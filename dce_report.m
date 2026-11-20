function txt = dce_report(result, file)
%DCE_REPORT Print (or return, or save) a plain-text summary of a design.
%
%   DCE_REPORT(result) prints the summary of a result of dce_design or
%   dce_evaluate. txt = DCE_REPORT(result) returns it as text instead.
%   DCE_REPORT(result, file) also writes it to a text file.
%
%   See also dce_design, dce_evaluate, dce_export.

L = {};
add = @(varargin) sprintf(varargin{:});
S = result.rows;
J = numel(result.alternatives);
nb = max(result.blocks);
L{end + 1} = add('Design "%s"  (%s, seed %d, spec %s)', result.name, result.toolbox, result.seed, result.specHash);
if nb > 1
    bs = accumarray(result.blocks, 1)';
    if all(bs == bs(1))
        btxt = add(', %d blocks of %d', nb, bs(1));
    else
        btxt = add(', %d blocks (%s)', nb, strtrim(sprintf('%d ', bs)));
    end
else
    btxt = '';
end
L{end + 1} = add('%d choice tasks x %d alternatives (%s)%s', S, J, strjoin(result.alternatives, ', '), btxt);
cname = 'D-error';
if strcmp(result.criterion, 'a'), cname = 'A-error'; end
if numel(result.models) > 1
    L{end + 1} = add('Objective (weighted %s): %.6g, weights %s', cname, result.objective, mat2str(result.weights, 4));
end
for m = 1:numel(result.models)
    r = result.models(m);
    L{end + 1} = ''; %#ok<*AGROW>
    L{end + 1} = add('Model %s (%s): D-error %.6g, A-error %.6g', r.name, r.type, r.Derror, r.Aerror);
    if ~isempty(r.bayes)
        b = r.bayes;
        L{end + 1} = add('  Bayesian %s over %d prior draws; at the prior means %.6g; median %.6g; 5%%-95%% %.4g-%.4g%s', ...
            r.stat, b.draws, r.DerrorFixed, b.median, b.q05, b.q95, ...
            ternary(b.shareNotIdentified > 0, add('; NOT identified for %.1f%% of draws', 100 * b.shareNotIdentified), ''));
    end
    if ~strcmp(r.type, 'mnl')
        how = r.method;
        if isfield(r, 'panelIntegration') && ~isempty(r.panelIntegration)
            if strcmp(r.panelIntegration.integration, 'importance')
                how = sprintf('%s; coefficients integrated by importance sampling for each sequence (%.0f%% effective draws)', how, 100 * r.panelIntegration.ess);
            else
                how = sprintf('%s; coefficients integrated with the same draws for all sequences', how);
            end
        end
        L{end + 1} = add('  Information: %s', how);
    end
    if ~r.identified
        L{end + 1} = '  NOT IDENTIFIED at the prior means: the information matrix is singular.';
        continue
    end
    L{end + 1} = add('  %-28s %10s %10s %8s %10s', 'Parameter', 'Prior', 's.e.(N=1)', 't(N=1)', 'N: t=1.96');
    for p = 1:numel(r.parameters)
        nm = r.parameters{p};
        if ~r.inCriterion(p), nm = [nm ' (not in criterion)']; end
        if isnan(r.Sestimate(p)), stxt = '-'; else, stxt = add('%.1f', r.Sestimate(p)); end
        L{end + 1} = add('  %-28s %10.5g %10.4g %8.3f %10s', nm, r.prior(p), r.se(p), r.tRatio(p), stxt);
    end
    if nb > 1
        L{end + 1} = add('  N counts replications of the whole design; multiply by %d for respondents (one block each).', nb);
    end
end
L{end + 1} = '';
sh = arrayfun(@(j) add('%s %.0f%%', result.alternatives{j}, 100 * result.meanShare(j)), 1:J, 'UniformOutput', false);
L{end + 1} = add('Predicted shares at the prior means (model 1): %s', strjoin(sh, ', '));
L{end + 1} = add('Utility balance %.0f%%; highest choice probability in a task: median %.2f, max %.2f', ...
    result.utilityBalance, median(result.maxProbability), max(result.maxProbability));
L{end + 1} = 'Level counts (column: count of each level, in the order of its levels):';
for k = 1:numel(result.columns)
    c = result.levelCounts(:, k);
    L{end + 1} = add('  %-28s %s', result.columns{k}, strtrim(sprintf('%d ', c(~isnan(c)))));
end
cr = result.correlations;
if ~any(isnan(cr(:))) && size(cr, 1) > 1
    off = abs(cr - diag(diag(cr)));
    [mx, ix] = max(off(:));
    [i1, i2] = ind2sub(size(off), ix);
    L{end + 1} = add('Largest correlation between design columns: %.2f (%s, %s)', mx, result.columns{i1}, result.columns{i2});
end
if ~isempty(result.blocking)
    b = result.blocking;
    L{end + 1} = add('Blocks: level balance Psi %.1f%% (upper bound %.1f%%); levels absent from a block: %d (%d avoidable)', ...
        b.Psi, b.PsiBound, b.absentLevels, b.absentAvoidable);
    if ~isempty(b.extremeColumns)
        L{end + 1} = add('  Lowest and highest level in every block (%s): %s', strjoin(b.extremeColumns, ', '), ternary(b.extremesOK, 'yes', 'NO'));
    end
    if ~isempty(b.leaveOneBlockOut)
        L{end + 1} = add('  D-error without one block, relative to all blocks: %s', strtrim(sprintf('%.3f ', b.leaveOneBlockOut)));
    end
end
ck = result.checks;
L{end + 1} = add('Checks: constraints %s; duplicate tasks %s; level frequencies %s', ...
    ternary(ck.constraintsOK, 'OK', add('VIOLATED in %d task(s)', size(ck.violations, 1))), ...
    ternary(isempty(ck.duplicateTasks), 'none', add('%d', numel(ck.duplicateTasks))), ...
    ternary(ck.frequenciesOK, 'OK', 'outside the bounds'));
if ~ck.constraintsOK
    for i = 1:min(10, size(ck.violations, 1))
        L{end + 1} = add('  task %d: %s', ck.violations(i, 1), ck.violationRules{ck.violations(i, 2)});
    end
end
if ~isempty(result.search) && isfield(result.search, 'startObjectives')
    s = result.search;
    best = min(s.startObjectives);
    L{end + 1} = add('Search: %d starts, %d evaluations, %.1f s; %d distinct optima; %d of %d starts within 0.1%% of the best', ...
        s.starts, s.evaluations, s.seconds, s.distinctOptima, nnz(s.startObjectives <= best * 1.001), s.starts);
elseif ~isempty(result.search) && isfield(result.search, 'construction')
    L{end + 1} = add('Construction: %s', result.search.construction);
end
for i = 1:numel(result.notes)
    L{end + 1} = ['Note: ' result.notes{i}];
end
if isfield(result, 'advice') && ~isempty(result.advice)
    L{end + 1} = '';
    L{end + 1} = 'Advice';
    for i = 1:numel(result.advice)
        a = result.advice(i);
        mark = '  - ';
        if strcmp(a.level, 'warning'), mark = '  ! '; end
        L{end + 1} = [mark a.text];
    end
end
out = strjoin(L, newline);
if nargin > 1 && ~isempty(file)
    fid = fopen(file, 'w', 'n', 'UTF-8');
    if fid < 0
        error('dce:export:file', 'Cannot write %s.', file);
    end
    fprintf(fid, '%s\n', out);
    fclose(fid);
end
if nargout > 0
    txt = out;
else
    fprintf('\n%s\n\n', out);
end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
