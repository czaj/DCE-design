function h = dd_spechash(spec)
%DD_SPECHASH Identify the design specification independently of search options.
try
    spec = rmfield(spec, intersect(fieldnames(spec), {'search','notes','source'}));
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(unicode2native(jsonencode(spec), 'UTF-8'));
    h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
    h = h(1:16);
catch
    h = '';
end
end
