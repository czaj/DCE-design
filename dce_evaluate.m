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
if ~isempty(blocks)
    blocks = double(blocks(:));
    if numel(blocks) ~= S
        dd_err('design:blocks', 'Blocks must have one entry per row of the design (%d).', S);
    end
    [~, ~, blocks] = unique(blocks);
    spec.blocks = max(blocks);
elseif isfield(spec, 'blocks')
    spec = rmfield(spec, 'blocks');
end
% level frequencies of the specification do not bind an existing design
freqSpec = isfield(spec, 'frequency') || ~isfield(spec, 'balance') || strcmpi(char(spec.balance), 'strict');
specEval = spec;
if isfield(specEval, 'frequency'), specEval = rmfield(specEval, 'frequency'); end
specEval.balance = 'none';
C = dd_compile(specEval);
o = dd_options(C.spec, varargin{:});
D = dd_tolevels(C, design, 'The design');
meas = dd_measures(C, D, blocks, o);
result = dd_result(C, D, blocks, o, meas);
if freqSpec
    result.notes{end + 1} = 'Evaluation only: the level-frequency rules of the specification are not applied; see levelCounts.';
end
if any(strcmp({C.model.type}, 'rppanel')) && isempty(blocks) && strcmp(C.panelUnit, 'block')
    result.notes{end + 1} = 'Panel mixed logit: no blocks were given, so a respondent is assumed to answer every row.';
end
result.advice = dd_advice(C, result);
if o.Verbose >= 1
    dce_report(result);
end
end
