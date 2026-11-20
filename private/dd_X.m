function X = dd_X(Mm, Dr, tl)
%DD_X Coded utility regressors of design rows. Internal.
%   X = DD_X(Mm, Dr) returns J-by-P-by-n: for each of the n rows of Dr (level
%   indices, n-by-K), the regressor of every coefficient in every alternative.
%   X = DD_X(Mm, Dr, tl) sums only the terms listed in tl.

n = size(Dr, 1);
X = zeros(Mm.J, Mm.P, n);
if nargin < 3
    tl = 1:Mm.T;
end
for t = tl
    v = Mm.tConst(t);
    for f = 1:Mm.tNF(t)
        tab = Mm.tTab{t, f};
        v = v .* tab(Dr(:, Mm.tCol(t, f)));
    end
    j = Mm.tAlt(t);
    p = Mm.tPar(t);
    X(j, p, :) = X(j, p, :) + reshape(v, 1, 1, []);
end
end
