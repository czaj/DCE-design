function q = dd_parseprior(v, n, where, isRandom)
%DD_PARSEPRIOR Parse prior values of one parameter (n coefficients). Internal.
%   Fixed or Bayesian priors: a number, 'N(mean,sd)' / 'normal(mean,sd)',
%   'U(lo,hi)' / 'uniform(lo,hi)', or Ngene style '(n,mean,sd)', '(u,lo,hi)'.
%   Several coefficients (dummy or effects coding): a numeric vector, a cell
%   array, or text with entries separated by |.
%   With isRandom = true the entries are distributions of random coefficients:
%   'normal(s)', 'lognormal(s)', '-lognormal(s)', 'uniform(s)', 'triangular(s)'
%   or 'fixed', where the spread s is a number or a Bayesian prior as above.
%
%   Output fields (n-by-1): dist (0 fixed, 1 normal, 2 uniform), a, b, mean;
%   for random coefficients also rdist (0 fixed, 1 normal, 2 lognormal,
%   3 uniform, 4 triangular) and rsign.

if nargin < 4
    isRandom = false;
end
if isnumeric(v) || islogical(v)
    ent = num2cell(double(v(:)'));
elseif ischar(v) || (isstring(v) && isscalar(v))
    ent = strtrim(strsplit(char(v), '|'));
elseif iscell(v)
    ent = v(:)';
elseif isstring(v)
    ent = cellstr(v(:))';
else
    dd_err('spec:priors', '%s: cannot read the value.', where);
end
if numel(ent) ~= n
    dd_err('spec:priors', '%s: %d value(s) given, %d needed (one per coefficient; for dummy or effects coding one per non-base level).', ...
        where, numel(ent), n);
end
q = struct('dist', zeros(n, 1), 'a', zeros(n, 1), 'b', zeros(n, 1), 'mean', zeros(n, 1), ...
    'rdist', zeros(n, 1), 'rsign', ones(n, 1));
for i = 1:n
    e = ent{i};
    if isRandom
        if isnumeric(e) || ~(ischar(e) || isstring(e))
            dd_err('spec:random', '%s: a random coefficient is written as distribution(spread), e.g. ''normal(0.3)''.', where);
        end
        s = strtrim(char(e));
        if strcmpi(s, 'fixed')
            continue
        end
        tk = regexp(s, '^(-?)\s*(normal|n|lognormal|ln|uniform|u|triangular|t)\s*\((.*)\)$', 'tokens', 'once', 'ignorecase');
        if isempty(tk)
            dd_err('spec:random', ['%s: cannot read ''%s''. Use normal(s), lognormal(s), -lognormal(s), uniform(s), ' ...
                'triangular(s) or fixed.'], where, s);
        end
        switch lower(tk{2})
            case {'normal', 'n'}, q.rdist(i) = 1;
            case {'lognormal', 'ln'}, q.rdist(i) = 2;
            case {'uniform', 'u'}, q.rdist(i) = 3;
            case {'triangular', 't'}, q.rdist(i) = 4;
        end
        if ~isempty(tk{1})
            if q.rdist(i) ~= 2
                dd_err('spec:random', '%s: a leading minus is only meaningful for lognormal coefficients.', where);
            end
            q.rsign(i) = -1;
        end
        [q.dist(i), q.a(i), q.b(i), q.mean(i)] = parseValue(tk{3}, where);
        lowest = q.a(i);
        if q.dist(i) == 1
            lowest = q.mean(i);
        end
        if ~(lowest >= 0) || ~(q.mean(i) > 0)
            dd_err('spec:random', '%s: the spread of a random coefficient must be positive.', where);
        end
    else
        [q.dist(i), q.a(i), q.b(i), q.mean(i)] = parseValue(e, where);
    end
end
end

function [dist, a, b, mu] = parseValue(e, where)
num = '([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)';
if isnumeric(e)
    if ~isscalar(e) || ~isfinite(e)
        dd_err('spec:priors', '%s: a prior must be a finite number.', where);
    end
    dist = 0; a = double(e); b = 0; mu = a;
    return
end
s = strtrim(char(e));
x = str2double(s);
if ~isnan(x) && isfinite(x)
    dist = 0; a = x; b = 0; mu = x;
    return
end
tk = regexp(s, ['^(normal|n|uniform|u)\s*\(\s*' num '\s*,\s*' num '\s*\)$'], 'tokens', 'once', 'ignorecase');
if isempty(tk)
    tk = regexp(s, ['^\(\s*(normal|n|uniform|u)\s*,\s*' num '\s*,\s*' num '\s*\)$'], 'tokens', 'once', 'ignorecase');
end
if isempty(tk)
    dd_err('spec:priors', '%s: cannot read ''%s''. Use a number, N(mean,sd) or U(lo,hi).', where, s);
end
a = str2double(tk{2});
b = str2double(tk{3});
if any(strcmpi(tk{1}, {'normal', 'n'}))
    dist = 1;
    mu = a;
    if ~(b > 0)
        dd_err('spec:priors', '%s: the standard deviation in ''%s'' must be positive.', where, s);
    end
else
    dist = 2;
    mu = (a + b) / 2;
    if ~(b > a)
        dd_err('spec:priors', '%s: the bounds in ''%s'' must be increasing.', where, s);
    end
end
end
