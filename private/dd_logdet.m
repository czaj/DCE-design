function ld = dd_logdet(A)
%DD_LOGDET Log-determinants of a stack of symmetric matrices. Internal.
%   A is P-by-P-by-(anything). ld has size(A) with the first two dimensions
%   removed (as a column for one page dimension). Pages that are not positive
%   definite return -Inf. Uses an LDL' elimination vectorised over pages.

sz = size(A);
P = sz(1);
N = numel(A) / (P * P);
A = reshape(A, P, P, N);
ld = zeros(1, 1, N);
bad = false(1, 1, N);
scale = zeros(1, 1, N);
for j = 1:P
    scale = max(scale, abs(A(j, j, :)));
end
tol = 1e-12 * scale;
for j = 1:P
    d = A(j, j, :);
    bad = bad | ~(d > tol);
    d(bad) = 1;
    ld = ld + log(d);
    if j < P
        v = A(j + 1:P, j, :);
        A(j + 1:P, j + 1:P, :) = A(j + 1:P, j + 1:P, :) - pagemtimes(v, 'none', v, 'transpose') ./ d;
    end
end
ld(bad) = -Inf;
if numel(sz) <= 3
    ld = reshape(ld, N, 1);
else
    ld = reshape(ld, [sz(3:end) 1]);
end
end
