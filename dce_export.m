function files = dce_export(result, folder, varargin)
%DCE_EXPORT Write a design and its evaluation to files.
%
%   files = DCE_EXPORT(result, folder) writes, into folder (created if needed):
%     <name>_design.csv   one row per choice task: Task, Block, and one column
%                         per attribute of each alternative, in level values
%     <name>_long.csv     one row per task and alternative: Task, Block,
%                         Alternative, one column per attribute (fixed levels
%                         of status-quo alternatives filled in; empty where an
%                         alternative does not have the attribute)
%     <name>.xlsx         sheets Design, Long, Labels (when spec.labels is
%                         given), Models, Probabilities, LevelCounts, Notes
%     <name>.json         everything in machine-readable form, including the
%                         specification, so that the design can be reproduced
%     <name>.mat          the result struct
%     <name>_report.txt   the text of dce_report
%   <name> is result.name with characters other than letters, digits, - and _
%   replaced. files is a struct with the paths written.
%
%   Options:
%     'Formats'    cell array out of {'csv','xlsx','json','mat','txt'} (default all)
%     'Name'       file name stem (default from result.name)
%     'Overwrite'  false (default): an existing file raises an error
%
%   See also dce_design, dce_evaluate, dce_report.

p = inputParser;
p.FunctionName = 'dce_export';
addParameter(p, 'Formats', {'csv', 'xlsx', 'json', 'mat', 'txt'});
addParameter(p, 'Name', '');
addParameter(p, 'Overwrite', false);
parse(p, varargin{:});
fmt = lower(cellstr(p.Results.Formats));
stem = char(p.Results.Name);
if isempty(stem)
    stem = regexprep(result.name, '[^A-Za-z0-9_\-]', '_');
end
if ~isfolder(folder)
    mkdir(folder);
end
C = dd_compile(evalSpec(result));
S = result.rows;
J = numel(result.alternatives);

% ---- long format -----------------------------------------------------------------
att = C.attrNames;
L = nan(S * J, numel(att));
for j = 1:J
    rows = (0:S - 1)' * J + j;
    for a = 1:numel(att)
        if C.scenCol(a) > 0
            L(rows, a) = C.LV{C.scenCol(a)}(result.levels(:, C.scenCol(a)));
        elseif C.colIndex(j, a) > 0
            k = C.colIndex(j, a);
            L(rows, a) = C.LV{k}(result.levels(:, k));
        elseif isfield(C.fixed, result.alternatives{j}) && isfield(C.fixed.(result.alternatives{j}), att{a})
            L(rows, a) = C.fixed.(result.alternatives{j}).(att{a});
        end
    end
end
long = table(repelem((1:S)', J), repelem(result.blocks, J), repmat(result.alternatives(:), S, 1), ...
    'VariableNames', {'Task', 'Block', 'Alternative'});
long = [long array2table(L, 'VariableNames', att)];

files = struct();
if any(strcmp(fmt, 'csv'))
    files.design = target(folder, [stem '_design.csv'], p.Results.Overwrite);
    writetable(result.design, files.design);
    files.long = target(folder, [stem '_long.csv'], p.Results.Overwrite);
    writetable(long, files.long);
end
if any(strcmp(fmt, 'txt'))
    files.report = target(folder, [stem '_report.txt'], p.Results.Overwrite);
    reportText = dce_report(result, files.report); %#ok<NASGU>  (assigned so that nothing is printed)
end
if any(strcmp(fmt, 'mat'))
    files.mat = target(folder, [stem '.mat'], p.Results.Overwrite);
    save(files.mat, 'result');
end
if any(strcmp(fmt, 'json'))
    files.json = target(folder, [stem '.json'], p.Results.Overwrite);
    fid = fopen(files.json, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s', jsonencode(jsonView(result), 'PrettyPrint', true));
    fclose(fid);
end
if any(strcmp(fmt, 'xlsx'))
    files.xlsx = target(folder, [stem '.xlsx'], p.Results.Overwrite);
    writetable(result.design, files.xlsx, 'Sheet', 'Design');
    writetable(long, files.xlsx, 'Sheet', 'Long');
    if ~isempty(fieldnames(C.labels))
        lab = long;
        for a = 1:numel(att)
            if isfield(C.labels, att{a})
                lv = C.globalLevels.(att{a});
                txt = repmat({''}, height(long), 1);
                [tf, loc] = ismember(L(:, a), lv);
                txt(tf) = C.labels.(att{a})(loc(tf));
                lab.(att{a}) = txt;
            end
        end
        writetable(lab, files.xlsx, 'Sheet', 'Labels');
    end
    M = table();
    for m = 1:numel(result.models)
        r = result.models(m);
        n = numel(r.parameters);
        M = [M; table(repmat({r.name}, n, 1), repmat({r.type}, n, 1), r.parameters(:), r.prior(:), r.se(:), r.tRatio(:), ...
            r.Sestimate(:), r.inCriterion(:), repmat(r.Derror, n, 1), repmat(r.Aerror, n, 1), ...
            'VariableNames', {'Model', 'Type', 'Parameter', 'Prior', 'SE_N1', 't_N1', 'N_for_t196', 'InCriterion', 'Derror', 'Aerror'})]; %#ok<AGROW>
    end
    writetable(M, files.xlsx, 'Sheet', 'Models');
    Pt = array2table([(1:S)' result.blocks result.probabilities], 'VariableNames', ...
        [{'Task', 'Block'}, strcat('P_', result.alternatives(:)')]);
    writetable(Pt, files.xlsx, 'Sheet', 'Probabilities');
    lc = result.levelCounts';
    lcT = array2table(lc, 'VariableNames', arrayfun(@(i) sprintf('Level%d', i), 1:size(lc, 2), 'UniformOutput', false));
    lcT = [table(result.columns(:), 'VariableNames', {'Column'}) lcT];
    writetable(lcT, files.xlsx, 'Sheet', 'LevelCounts');
    info = [{sprintf('Design %s', result.name)}; {result.toolbox}; {['Created ' result.created]}; ...
        {sprintf('Seed %d, specification hash %s', result.seed, result.specHash)}; ...
        {sprintf('Objective %.8g', result.objective)}; result.notes(:)];
    writetable(table(info, 'VariableNames', {'Notes'}), files.xlsx, 'Sheet', 'Notes');
end
end

% =============================================================================
function f = target(folder, name, overwrite)
f = fullfile(folder, name);
if isfile(f)
    if overwrite
        delete(f);
    else
        error('dce:export:exists', '%s already exists. Choose another Name or folder, or set Overwrite to true.', f);
    end
end
end

function spec = evalSpec(result)
% Specification compiled as in dce_evaluate (level frequencies not enforced).
spec = result.spec;
spec.rows = result.rows;
if isfield(spec, 'frequency'), spec = rmfield(spec, 'frequency'); end
spec.balance = 'none';
if max(result.blocks) > 1
    spec.blocks = max(result.blocks);
elseif isfield(spec, 'blocks')
    spec = rmfield(spec, 'blocks');
end
end

function J = jsonView(result)
J = struct();
keep = {'name', 'toolbox', 'created', 'specHash', 'seed', 'rows', 'alternatives', 'columns', 'criterion', 'weights', ...
    'objective', 'meanShare', 'utilityBalance', 'notes'};
for i = 1:numel(keep)
    J.(keep{i}) = result.(keep{i});
end
J.design = table2struct(result.design);
J.blocks = result.blocks(:)';
J.probabilities = result.probabilities;
J.levelCounts = result.levelCounts';
J.checks = result.checks;
J.blocking = result.blocking;
mm = result.models;
for m = 1:numel(mm)
    mm(m).information = squeeze(mm(m).information);
end
J.models = mm;
if ~isempty(result.search)
    s = result.search;
    s.candidates = table2struct(s.candidates);
    J.search = s;
end
sp = result.spec;
if isfield(sp, 'start') && ~(ischar(sp.start) || isstring(sp.start))
    sp = rmfield(sp, 'start');
end
J.spec = sp;
end
