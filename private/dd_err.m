function dd_err(id, fmt, varargin)
%DD_ERR Throw an error with identifier dce:<id>. Internal.
if isempty(varargin)
    error(['dce:' id], '%s', fmt);
else
    error(['dce:' id], fmt, varargin{:});
end
end
