function result = dce_evaluate(spec, design, varargin)
%DCE_EVALUATE Evaluate an existing design under a specification.
%
%   result = DCE_EVALUATE(spec, design) computes for a given design everything
%   that dce_design reports: D-error and A-error of every model (at the prior
%   means and over the prior draws), information and covariance matrices,
%   standard errors and S-estimates for one replication, predicted choice
%   shares, level counts, correlations, and whether the design satisfies the
%   constraints of spec. Nothing is searched.
%
%   spec    a struct, a .json file or an Ngene .ngs script (see dce_design).
%   design  level values, one row per choice task:
%           - a numeric matrix with the columns in the order of
%             result.columns (scenario variables first, then the attributes of
%             each alternative in the order of spec.attributes), or
%           - a table with variables named 'alt_attr' (or 'alt.attr'); a
%             variable 'Block' is used as the blocking; others are ignored, or
%           - a result struct of dce_design, or a .csv / .xlsx / .ngd file.
%
%   Options (name-value, as in dce_design): 'Blocks' (vector), 'Seed',
%   'EvalBayesDraws', 'EvalRandomDraws', 'PanelRespondents', 'Verbose', ...
%   spec.rows is replaced by the number of rows of the design.
%
%   Use it to compare designs under the same priors, to check a design made
%   elsewhere (for example in Ngene), or to re-evaluate a fielded design under
%   updated priors.
%
%   See also dce_design, dce_block, dce_report.

spec = dd_readspec(spec);
blocks = [];
keep = true(1, numel(varargin));
for i = 1:2:numel(varargin)
    if strcmpi(char(varargin{i}), 'Blocks')
        blocks = varargin{i + 1};
        keep([i i + 1]) = false;
    end
end
varargin = varargin(keep);

if ischar(design) || (isstring(design) && isscalar(design))
    file = char(design);
    [~, ~, ext] = fileparts(file);
    if strcmpi(ext, '.ngd')
        design = dce_readngd(file);
    else
        design = readtable(file, 'VariableNamingRule', 'preserve');
    end
elseif isstruct(design) && isfield(design, 'design')
    design = design.design;
end
if istable(design) && isempty(blocks)
    hit = find(strcmpi(design.Properties.VariableNames, 'Block'), 1);
    if ~isempty(hit)
        blocks = design{:, hit};
    end
end
S = size(design, 1);
spec.rows = S;
requestedBlocks = [];
if isfield(spec,'blocks'), requestedBlocks = spec.blocks; end
if ~isempty(blocks)
    if ~isnumeric(blocks) || ~isreal(blocks) || ~isvector(blocks) || any(~isfinite(blocks))
        dd_err('design:blocks', 'Block labels must be a vector of finite numeric values.');
    end
    blocks = double(blocks(:));
    if numel(blocks) ~= S
        dd_err('design:blocks', 'Blocks must have one entry per row of the design (%d).', S);
    end
    [~, ~, blocks] = unique(blocks);
    spec.blocks = max(blocks);
elseif isfield(spec, 'blocks')
    spec = rmfield(spec, 'blocks');
end
% Report violations even if the supplied design has too few rows to meet the bounds.
C = dd_compile(spec, true);
panelByBlock = strcmp(C.panelUnit,'block') && any(strcmp({C.model.type},'rppanel'));
if ~isempty(requestedBlocks) && requestedBlocks > 1 && panelByBlock && isempty(blocks)
    dd_err('design:blocks', 'This panel specification assigns one block per respondent. Supply a Block column or the Blocks vector.');
end
o = dd_options(C.spec, varargin{:});
D = dd_tolevels(C, design, 'The design');
meas = dd_measures(C, D, blocks, o);
result = dd_result(C, D, blocks, o, meas);
if ~result.checks.frequenciesOK
    result.notes{end + 1} = 'Observed level frequencies are outside the specification bounds; compare levelCounts with the specified frequencies.';
end
if panelByBlock && isempty(blocks)
    result.notes{end + 1} = 'Panel mixed logit: no blocks were given, so a respondent is assumed to answer every row.';
end
if ~isempty(blocks) && ~isempty(requestedBlocks) && requestedBlocks ~= max(blocks)
    result.notes{end+1} = sprintf('The file has %d block(s); the specification requested %d. Evaluation uses the block labels in the file.',max(blocks),requestedBlocks);
end
result.advice = dd_advice(C, result);
if o.Verbose >= 1
    dce_report(result);
end
end
