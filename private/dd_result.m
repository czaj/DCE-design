function R = dd_result(C, D, y, o, meas)
%DD_RESULT Assemble the result structure of a design. Internal.
%   Rows are sorted by block. meas is the output of dd_measures for (D, y).

S = size(D, 1);
if isempty(y)
    y = ones(S, 1);
end
V = zeros(S, C.K);
for k = 1:C.K
    V(:, k) = C.LV{k}(D(:, k));
end
vn = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(strrep(C.colName, '.', '_')), {'Task', 'Block'});
T = array2table([(1:S)' y V], 'VariableNames', [{'Task', 'Block'}, vn]);

R = struct();
R.name = C.name;
R.toolbox = ['DCE-design ' C.version];
R.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
R.specHash = C.hash;
R.seed = o.Seed;
R.design = T;
R.columns = C.colName;
R.levels = D;
R.blocks = y;
R.rows = S;
R.alternatives = C.altNames;
R.criterion = C.criterion;
R.weights = C.w';
R.objective = meas.objective;
R.models = meas.models;
R.probabilities = meas.probabilities;
R.meanShare = meas.meanShare;
R.utilityBalance = meas.utilityBalance;
R.maxProbability = meas.maxProbability;
R.levelCounts = meas.levelCounts;
R.correlations = meas.correlations;
R.checks = meas.checks;
R.blocking = [];
R.search = [];
R.advice = [];
R.notes = C.notes(:);
R.spec = C.spec;
end
