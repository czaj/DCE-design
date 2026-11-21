function spec = dd_readspec(spec)
%DD_READSPEC Accept a specification as a struct, a JSON file or an Ngene script. Internal.
if isstruct(spec)
    return
end
if ~(ischar(spec) || (isstring(spec) && isscalar(spec)))
    dd_err('spec:type', 'The specification must be a struct, a .json file or an .ngs file.');
end
file = char(spec);
if ~isfile(file)
    dd_err('spec:file', 'Specification file not found: %s', file);
end
[~, ~, ext] = fileparts(file);
switch lower(ext)
    case '.json'
        try
            spec = jsondecode(fileread(file));
        catch err
            dd_err('spec:json', 'Cannot read %s as JSON (%s).', file, err.message);
        end
        if ~isstruct(spec) || ~isscalar(spec)
            dd_err('spec:json', '%s must contain one JSON object.', file);
        end
    case {'.ngs', '.ngd'}
        spec = dce_ngene(file);
    otherwise
        dd_err('spec:file', 'Unknown specification file type ''%s''. Use .json, .ngs or .ngd with an embedded script.', ext);
end
spec.source = file;
end
