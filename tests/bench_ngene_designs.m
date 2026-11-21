function T = bench_ngene_designs(files, outFile, varargin)
%BENCH_NGENE_DESIGNS Compare designs saved by Ngene with designs from dce_design.
%
%   T = BENCH_NGENE_DESIGNS(files) takes Ngene .ngd files (each holds a design
%   and the script that produced it), translates the script with dce_ngene
%   (Ngene's conventions kept), evaluates Ngene's design, generates a design
%   for the same specification and evaluates it in the same way. files is a
%   cell array of paths or a folder (searched recursively).
%   BENCH_NGENE_DESIGNS(files, outFile) also writes the table as CSV.
%   Further arguments are options of dce_design.
%
%   Columns: file, rows, blocks, models, status, D-error (or weighted
%   objective) of Ngene's design and of the generated one, gain in %, seconds,
%   and whether Ngene's design keeps this specification's constraints.

if ischar(files) || isstring(files)
    if isscalar(string(files)) && isfolder(files)
        d = dir(fullfile(char(files), '**', '*.ngd'));
        files = fullfile({d.folder}, {d.name});
    else
        files = cellstr(files);
    end
end
n = numel(files);
T = table('Size', [n 12], 'VariableTypes', {'string', 'double', 'double', 'string', 'string', 'double', 'double', 'double', 'double', 'logical', 'string', 'string'}, ...
    'VariableNames', {'File', 'Rows', 'Blocks', 'Models', 'Status', 'NgeneObjective', 'NewObjective', 'GainPercent', 'Seconds', 'NgeneKeepsConstraints', 'NgeneModelValues', 'NewModelValues'});
for i = 1:n
    [~, nm] = fileparts(files{i});
    T.File(i) = nm;
    T.Status(i) = "ok";
    try
        [des, script] = dce_readngd(files{i});
        if isempty(script)
            error('dce:bench:noscript', 'no script in the file');
        end
        spec = dce_ngene(script, 'Emulate', true);   % Ngene's own criterion, so that both designs are judged as Ngene judged its own
        spec.name = nm;
        if isfield(spec, 'start'), spec = rmfield(spec, 'start'); end
        rN = dce_evaluate(spec, des, varargin{:}, 'Verbose', 0);
        t = tic;
        r = dce_design(spec, 'Verbose', 0, varargin{:});
        T.Seconds(i) = toc(t);
        T.Rows(i) = r.rows;
        T.Blocks(i) = max(r.blocks);
        T.Models(i) = strjoin(arrayfun(@(m) sprintf('%s:%s', m.type, m.stat), r.models, 'UniformOutput', false), ' + ');
        T.NgeneObjective(i) = rN.objective;
        T.NewObjective(i) = r.objective;
        T.GainPercent(i) = 100 * (1 - r.objective / rN.objective);
        T.NgeneKeepsConstraints(i) = rN.checks.constraintsOK && isempty(rN.checks.duplicateTasks) && rN.checks.frequenciesOK;
        T.NgeneModelValues(i) = strjoin(arrayfun(@(m) sprintf('%.5g', m.criterion), rN.models, 'UniformOutput', false), ' ');
        T.NewModelValues(i) = strjoin(arrayfun(@(m) sprintf('%.5g', m.criterion), r.models, 'UniformOutput', false), ' ');
    catch err
        T.Status(i) = string(sprintf('%s: %s', err.identifier, err.message));
    end
    fprintf('%-45s %s  Ngene %.5g  new %.5g  gain %.1f%%  %.0f s\n', nm, extractBefore(T.Status(i) + "          ", 11), ...
        T.NgeneObjective(i), T.NewObjective(i), T.GainPercent(i), T.Seconds(i));
end
if nargin > 1 && ~isempty(outFile)
    writetable(T, outFile);
end
end
