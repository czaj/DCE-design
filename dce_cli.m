function summary = dce_cli(specFile, outFolder, varargin)
%DCE_CLI Generate (or evaluate) a design from a file and write all outputs.
%
%   DCE_CLI(specFile, outFolder) reads a specification (.json, or an Ngene
%   .ngs script), runs dce_design, writes every output of dce_export into
%   outFolder and prints one line of JSON that starts with DCE_RESULT, e.g.
%
%     DCE_RESULT {"ok":true,"name":"wave2","objective":0.0862,...}
%
%   On failure it prints DCE_RESULT {"ok":false,"error":"dce:spec:...",
%   "message":"..."} and rethrows, so that  matlab -batch  exits with a
%   non-zero status. This is the entry point meant for scripts and AI agents:
%
%     matlab -batch "addpath('<toolbox>'); dce_cli('spec.json','out')"
%
%   DCE_CLI(specFile, outFolder, 'Evaluate', designFile) evaluates an existing
%   design (.csv, .xlsx or .ngd) instead of searching.
%   Other name-value options are passed to dce_design / dce_evaluate, and
%   'Name' and 'Overwrite' to dce_export.
%
%   See also dce_design, dce_evaluate, dce_export.

evalFile = '';
exportArgs = {};
keep = true(1, numel(varargin));
for i = 1:2:numel(varargin)
    switch lower(char(varargin{i}))
        case 'evaluate'
            evalFile = char(varargin{i + 1});
            keep([i i + 1]) = false;
        case {'name', 'overwrite'}
            exportArgs = [exportArgs varargin(i:i + 1)]; %#ok<AGROW>
            keep([i i + 1]) = false;
    end
end
args = varargin(keep);
try
    if isempty(evalFile)
        r = dce_design(specFile, args{:}, 'Verbose', 0);
    else
        r = dce_evaluate(specFile, evalFile, args{:}, 'Verbose', 0);
    end
    files = dce_export(r, outFolder, exportArgs{:});
    summary = struct('ok', true, 'name', r.name, 'rows', r.rows, 'blocks', max(r.blocks), ...
        'criterion', r.criterion, 'objective', r.objective, ...
        'models', {arrayfun(@(m) struct('name', m.name, 'type', m.type, 'Derror', m.Derror, 'DerrorAtPriorMeans', m.DerrorFixed, ...
        'Aerror', m.Aerror, 'identified', m.identified), r.models)}, ...
        'constraintsOK', r.checks.constraintsOK, 'duplicateTasks', numel(r.checks.duplicateTasks), ...
        'frequenciesOK', r.checks.frequenciesOK, 'seed', r.seed, 'specHash', r.specHash, ...
        'notes', {r.notes}, 'files', files);
    summary.advice = arrayfun(@(a) sprintf('%s (%s): %s', a.topic, a.level, a.text), r.advice, 'UniformOutput', false);
    if ~isempty(r.search) && isfield(r.search, 'starts')
        summary.seconds = r.search.seconds;
        summary.starts = r.search.starts;
        summary.method = r.search.method;
    elseif ~isempty(r.search) && isfield(r.search, 'construction')
        summary.construction = r.search.construction;
    end
    if ~isempty(r.blocking)
        summary.blockBalancePsi = r.blocking.Psi;
        summary.levelsAbsentFromABlock = r.blocking.absentLevels;
    end
    fprintf('DCE_RESULT %s\n', jsonencode(summary));
catch err
    summary = struct('ok', false, 'error', err.identifier, 'message', err.message);
    fprintf('DCE_RESULT %s\n', jsonencode(summary));
    rethrow(err);
end
end
