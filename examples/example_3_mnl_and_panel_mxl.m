% EXAMPLE_3_MNL_AND_PANEL_MXL  One design for two models: MNL and panel mixed logit.
%
% The design should work both for a multinomial logit (priors from an MNL
% estimated on the previous wave) and for a panel mixed logit with normal
% coefficients (priors from an MXL). Each model's D-error is divided by the
% best value found for that model alone, so equal weights mean equal
% importance. 24 tasks in 3 blocks of 8: a respondent answers one block, which
% is what the panel model is evaluated on. Takes a couple of minutes.

addpath(fileparts(fileparts(mfilename('fullpath'))));

spec = struct();
spec.name = 'food_labels';
spec.alternatives = {'a', 'b', 'none'};
spec.rows = 24;
spec.blocks = 3;
spec.attributes.organic = [0 1];
spec.attributes.label   = [0 1 2];
spec.attributes.price   = [0.8 1.0 1.2 1.4 1.6 1.8];

u = 'b_org*organic + b_lab.dummy*label + b_price*price';
U = struct('a', u, 'b', u, 'none', 'asc_none');

mnl = struct('name', 'MNL', 'type', 'mnl', 'utility', U);
mnl.priors = struct('b_org', 'N(0.37, 0.15)', 'b_lab', {{'N(0.21, 0.10)', 'N(0.50, 0.15)'}}, ...
    'b_price', 'N(-1.9, 0.5)', 'asc_none', 'N(-2.3, 0.6)');

mxl = struct('name', 'MXL', 'type', 'rppanel', 'utility', U);
mxl.priors = struct('b_org', 0.5, 'b_lab', [0.3 0.6], 'b_price', -3.5, 'asc_none', -5);
mxl.random = struct('b_org', 'normal(0.6)', 'b_lab', {{'normal(0.8)', 'normal(0.8)'}}, ...
    'b_price', 'normal(2)', 'asc_none', 'normal(3)');

spec.models = {mnl, mxl};
spec.weights = [1 1];

result = dce_design(spec);
fprintf('MNL D-error %.4f; panel MXL D-error %.4f (%s)\n', result.models(1).Derror, result.models(2).Derror, result.models(2).method);
