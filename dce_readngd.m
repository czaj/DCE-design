function [design, script, all] = dce_readngd(file, which)
%DCE_READNGD Read a design saved by Ngene (.ngd file).
%
%   [design, script] = DCE_READNGD(file) returns the design as a table (one
%   row per choice situation; variables named as in Ngene, 'alt.attribute',
%   plus 'Block' when present) and the Ngene script stored in the same file.
%   An .ngd file can hold several designs; DCE_READNGD(file, n) returns the
%   n-th (default 1). The third output is the table of all designs with the
%   variable 'Design'.
%
%   Typical use:
%     [d, s] = dce_readngd('my design.ngd');
%     r = dce_evaluate(dce_ngene(s), d);      % evaluate Ngene's design
%
%   See also dce_ngene, dce_evaluate.

if nargin < 2
    which = 1;
end
if ~isfile(file)
    error('dce:ngd:file', 'File not found: %s', file);
end
txt = fileread(file);
txt = strrep(txt, sprintf('\r'), '');
lines = strsplit(txt, newline);
sep = find(~cellfun(@isempty, regexp(lines, '^\|+\s*$', 'once')), 1);
if isempty(sep)
    tabLines = lines;
    script = '';
else
    tabLines = lines(1:sep - 1);
    script = strtrim(strjoin(lines(sep + 1:end), newline));
end
tabLines = tabLines(~cellfun(@(s) isempty(strtrim(s)), tabLines));
if numel(tabLines) < 2
    if isempty(script)
        error('dce:ngd:format', '%s contains neither a design table nor a script.', file);
    end
    design = table();
    all = table();
    return
end
head = strtrim(strsplit(strtrim(tabLines{1}), sprintf('\t')));
head = head(~cellfun(@isempty, head));
n = numel(tabLines) - 1;
A = nan(n, numel(head));
for i = 1:n
    v = sscanf(strrep(tabLines{i + 1}, ',', '.'), '%f')';
    if numel(v) ~= numel(head)
        error('dce:ngd:format', '%s: row %d of the design has %d values, the header has %d columns.', file, i, numel(v), numel(head));
    end
    A(i, :) = v;
end
all = array2table(A, 'VariableNames', head);
iDes = find(strcmpi(head, 'Design'), 1);
iSit = find(strcmpi(head, 'Choice situation'), 1);
rows = true(n, 1);
if ~isempty(iDes)
    ids = unique(A(:, iDes), 'stable');
    if which > numel(ids)
        error('dce:ngd:design', '%s holds %d design(s); design %d was requested.', file, numel(ids), which);
    end
    rows = A(:, iDes) == ids(which);
end
keep = true(1, numel(head));
keep([iDes iSit]) = false;
design = all(rows, keep);
end
