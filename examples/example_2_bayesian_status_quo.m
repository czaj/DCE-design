% EXAMPLE_2_BAYESIAN_STATUS_QUO  Bayesian design with a status quo, constraints and blocks.
%
% A status quo and two programmes. Priors come from a pilot (mean and about
% 1.5 x standard error). A programme must improve on the status quo, and a
% programme that is at least as good as the other on every attribute must
% cost more. 24 tasks in 3 blocks. Takes about a minute.

addpath(fileparts(fileparts(mfilename('fullpath'))));

spec = struct();
spec.name = 'water_main';
spec.alternatives = {'sq', 'a', 'b'};
spec.rows = 24;
spec.blocks = 3;

spec.attributes.quality = [0 1 2];
spec.attributes.access  = [0 1];
spec.attributes.cost    = [0 20 40 80 120];
spec.altLevels.a.cost = [20 40 80 120];
spec.altLevels.b.cost = [20 40 80 120];
spec.fixed.sq = struct('quality', 0, 'access', 0, 'cost', 0);

u = 'b_q.dummy*quality + b_acc*access + b_cost*cost';
spec.model.utility.sq = ['asc_sq + ' u];
spec.model.utility.a = u;
spec.model.utility.b = u;
spec.model.priors.asc_sq = 'N(0.2, 0.25)';
spec.model.priors.b_q    = {'N(0.35, 0.15)', 'N(0.60, 0.15)'};
spec.model.priors.b_acc  = 'N(0.30, 0.12)';
spec.model.priors.b_cost = 'N(-0.012, 0.004)';

spec.reject = {'a.quality == 0 & a.access == 0', 'b.quality == 0 & b.access == 0'};
% dominance between a and b is checked automatically (same utility structure)

result = dce_design(spec);

% the design that is being replaced can be evaluated under the same priors:
%   old = readtable('old_design.csv');
%   rOld = dce_evaluate(spec, old);
%   fprintf('D-error %.4f -> %.4f\n', rOld.models.Derror, result.models.Derror);
