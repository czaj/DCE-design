% EXAMPLE_5_TYPES_AND_ADVICE  Full factorial, orthogonal and efficient designs
% for the same specification, and the advice the toolbox gives.
%
% Two alternatives with two 3-level attributes and one 2-level attribute.
% Takes a few seconds.

addpath(fileparts(fileparts(mfilename('fullpath'))));

spec = struct();
spec.name = 'three_designs';
spec.alternatives = {'a', 'b'};
spec.rows = 18;
spec.attributes.time    = [10 20 30];
spec.attributes.comfort = [0 1 2];
spec.attributes.wifi    = [0 1];
u = 'b_time*time + b_comfort.dummy*comfort + b_wifi*wifi';
spec.model.utility = struct('a', u, 'b', u);
spec.model.priors = struct('b_time', -0.04, 'b_comfort', [0.3 0.5], 'b_wifi', 0.25);

% what the toolbox can say before any search: size of the problem, rows and
% blocks against level balance, priors, constraints
dce_advise(spec);

% all feasible choice tasks (dominated and repeated ones are left out)
full = dce_full(spec, 'Verbose', 0);

% orthogonal within each alternative: 18 rows
orth = dce_orthogonal(spec, 'Type', 'sequential', 'Verbose', 0);

% efficient: uses the priors, respects the constraints
eff = dce_design(spec, 'Verbose', 0);

fprintf('\n%-34s %6s %10s %s\n', 'design', 'rows', 'D-error', 'tasks with a dominated alternative');
fprintf('%-34s %6d %10.4f %d\n', 'all feasible tasks', full.rows, full.models.Derror, size(full.checks.violations, 1));
fprintf('%-34s %6d %10.4f %d\n', 'orthogonal (sequential)', orth.rows, orth.models.Derror, size(orth.checks.violations, 1));
fprintf('%-34s %6d %10.4f %d\n', 'efficient', eff.rows, eff.models.Derror, size(eff.checks.violations, 1));
fprintf('The D-error falls with the number of rows, so compare designs of the same size.\n\n');

% remarks on the finished design: sample size, how close it is to the best
% possible, which search method was used and why
dce_advise(eff);
