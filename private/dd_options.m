function o = dd_options(spec, varargin)
%DD_OPTIONS Search and evaluation options. Internal.
%   Defaults, overridden by spec.search (a struct), overridden by name-value
%   arguments. Names are not case sensitive. Unknown names raise an error.

o = struct( ...
    'Seed', 1, ...                 % one seed fixes starts and draws
    'Restarts', [], ...            % independent random starts ([] = 12; 4 for costly problems: Bayesian priors, rp models, 25+ parameters, rows x columns >= 600)
    'Kicks', [], ...               % perturb-and-reoptimise rounds per start ([] = 3; 1 for costly problems)
    'Method', 'auto', ...          % 'auto', 'swap', 'rowexchange', 'both' or 'exhaustive' (see README: how the search works)
    'CandidateMax', 3000, ...      % candidate tasks for row exchange (all feasible tasks if there are no more)
    'ExhaustiveMax', 2e6, ...      % largest number of designs that Method 'exhaustive' will evaluate
    'KickSize', 4, ...             % random swaps per perturbation
    'MaxPasses', 200, ...          % passes of the local search
    'Tol', 1e-9, ...               % relative improvement that counts
    'MaxTime', Inf, ...            % seconds; when reached, the best design so far is returned
    'Parallel', 'auto', ...        % 'auto' uses an open parallel pool; true opens one; false never
    'BayesDraws', 128, ...         % prior draws during the search
    'EvalBayesDraws', 2048, ...    % prior draws for the reported values
    'ScreenDraws', 8, ...          % prior draws on which moves are ranked before full evaluation
    'ScreenTop', 4, ...            % moves per step that are evaluated on all prior draws
    'FullCheck', false, ...        % finish with a pass that evaluates every move on all prior draws
    'RandomDraws', 128, ...        % random-coefficient draws during the search
    'EvalRandomDraws', 1000, ...   % random-coefficient draws for the reported values
    'DrawType', 'sobol', ...       % 'sobol' or 'halton'
    'PanelSurrogate', 'mnl', ...   % how rppanel models enter the first search stage: 'mnl' or 'rp'
    'PanelPolish', 3, ...          % passes of the second stage on the panel criterion itself (0 = none)
    'PanelIntegration', 'importance', ... % panel information: 'importance' (draws around the mode of each choice sequence) or 'prior' (the same draws for all sequences)
    'PanelDraws', 32, ...          % draws per choice sequence during the search (importance sampling)
    'EvalPanelDraws', 128, ...     % draws per choice sequence for the reported values
    'PanelExactMax', 5000, ...     % enumerate choice sequences up to this many per respondent
    'PanelRespondents', 500, ...   % simulated respondents during the search, when sequences are not enumerated
    'EvalPanelRespondents', 2000, ... % simulated respondents for the reported values
    'PanelBayesDraws', 16, ...     % prior draws used for panel evaluations
    'Normalize', [], ...           % scale each model by its own best value ([] = when several models)
    'Candidates', 3, ...           % best distinct designs that are blocked and fully evaluated
    'Block', true, ...             % assign blocks when spec.blocks > 1
    'BlockRestarts', 30, ...       % random starts of the blocking search (dce_block)
    'Start', [], ...               % starting design (level values, S-by-K) or a result struct
    'Checkpoint', '', ...          % file that receives the best design after every start
    'Verbose', 1);

names = fieldnames(o);
if isfield(spec, 'search') && ~isempty(spec.search)
    if ~isstruct(spec.search)
        dd_err('options:type', 'spec.search must be a struct of options.');
    end
    f = fieldnames(spec.search);
    for i = 1:numel(f)
        o = setOne(o, names, f{i}, spec.search.(f{i}));
    end
end
if mod(numel(varargin), 2) ~= 0
    dd_err('options:type', 'Options must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    o = setOne(o, names, char(varargin{i}), varargin{i + 1});
end

if ~isempty(o.Restarts) && (~isnumeric(o.Restarts) || ~isscalar(o.Restarts) || o.Restarts < 1 || o.Restarts ~= round(o.Restarts))
    dd_err('options:value', 'Option Restarts must be a positive integer.');
end
intPos = {'KickSize', 'MaxPasses', 'BayesDraws', 'EvalBayesDraws', 'RandomDraws', ...
    'EvalRandomDraws', 'PanelExactMax', 'PanelBayesDraws', 'Candidates', 'BlockRestarts', ...
    'ScreenDraws', 'ScreenTop', 'PanelDraws', 'EvalPanelDraws'};
for i = 1:numel(intPos)
    v = o.(intPos{i});
    if ~isnumeric(v) || ~isscalar(v) || v < 1 || v ~= round(v)
        dd_err('options:value', 'Option %s must be a positive integer.', intPos{i});
    end
end
if ~isempty(o.Kicks) && (~isnumeric(o.Kicks) || ~isscalar(o.Kicks) || o.Kicks < 0 || o.Kicks ~= round(o.Kicks))
    dd_err('options:value', 'Option Kicks must be a non-negative integer.');
end
if ~isnumeric(o.PanelPolish) || ~isscalar(o.PanelPolish) || o.PanelPolish < 0 || o.PanelPolish ~= round(o.PanelPolish)
    dd_err('options:value', 'Option PanelPolish must be a non-negative integer.');
end
if ~isnumeric(o.PanelRespondents) || ~isscalar(o.PanelRespondents) || o.PanelRespondents < 0 || ...
        ~isnumeric(o.EvalPanelRespondents) || ~isscalar(o.EvalPanelRespondents) || o.EvalPanelRespondents < 0
    dd_err('options:value', 'Options PanelRespondents and EvalPanelRespondents must be non-negative integers.');
end
if ~isnumeric(o.Seed) || ~isscalar(o.Seed) || o.Seed < 0 || o.Seed ~= round(o.Seed)
    dd_err('options:value', 'Option Seed must be a non-negative integer.');
end
o.DrawType = lower(char(o.DrawType));
if ~any(strcmp(o.DrawType, {'sobol', 'halton'}))
    dd_err('options:value', 'Option DrawType must be ''sobol'' or ''halton''.');
end
o.Method = lower(char(o.Method));
if ~any(strcmp(o.Method, {'auto', 'swap', 'rowexchange', 'both', 'exhaustive'}))
    dd_err('options:value', 'Option Method must be ''auto'', ''swap'', ''rowexchange'', ''both'' or ''exhaustive''.');
end
if ~isnumeric(o.CandidateMax) || ~isscalar(o.CandidateMax) || o.CandidateMax < 10 || o.CandidateMax ~= round(o.CandidateMax)
    dd_err('options:value', 'Option CandidateMax must be an integer of at least 10.');
end
o.method = o.Method;
o.PanelIntegration = lower(char(o.PanelIntegration));
if ~any(strcmp(o.PanelIntegration, {'importance', 'prior'}))
    dd_err('options:value', 'Option PanelIntegration must be ''importance'' or ''prior''.');
end
if o.PanelDraws < 16 || o.EvalPanelDraws < 16
    dd_err('options:value', 'Options PanelDraws and EvalPanelDraws must be at least 16.');
end
o.PanelSurrogate = lower(char(o.PanelSurrogate));
if ~any(strcmp(o.PanelSurrogate, {'mnl', 'rp'}))
    dd_err('options:value', 'Option PanelSurrogate must be ''mnl'' or ''rp''.');
end
if ~(isnumeric(o.MaxTime) && isscalar(o.MaxTime) && o.MaxTime > 0)
    dd_err('options:value', 'Option MaxTime must be a positive number of seconds.');
end
o.deadline = posixtime(datetime('now')) + o.MaxTime;
end

function o = setOne(o, names, name, value)
hit = find(strcmpi(name, names));
if isempty(hit)
    dd_err('options:unknown', 'Unknown option ''%s''. Options: %s.', name, strjoin(names', ', '));
end
o.(names{hit}) = value;
end
