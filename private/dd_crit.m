function [f, fr] = dd_crit(It, Mm, crit)
%DD_CRIT Efficiency criterion from information matrices. Internal.
%   It is Pt-by-Pt-by-n-by-R (n candidate designs, R prior draws).
%   f is n-by-1: the mean over draws of the D-error (crit 'd') or A-error
%   (crit 'a'); fr is n-by-R with the value for each draw. A draw for which the
%   parameters are not identified gives Inf.
%
%   D-error = det(Omega)^(1/P), Omega = inv(I), for one respondent who answers
%   every choice task. Parameters with Mm.eff false are left out: as nuisance
%   parameters (Omega is the block of the full inverse) or as known constants
%   (their rows and columns are dropped from I).

sz = size(It);
n = size(It, 3);
R = size(It, 4);
eff = Mm.eff;
allIn = all(eff);
if ~allIn && strcmp(Mm.exMode, 'known')
    It = It(eff, eff, :, :);
    allIn = true;
end
Pe = Mm.Peff;
switch crit
    case 'd'
        ld = reshape(dd_logdet(It), n, R);
        if ~allIn
            ld = ld - reshape(dd_logdet(It(~eff, ~eff, :, :)), n, R);
        end
        fr = exp(-ld / Pe);
    case 'a'
        ld = reshape(dd_logdet(It), n, R);
        fr = inf(n, R);
        ok = isfinite(ld);
        if any(ok(:))
            A = reshape(It, size(It, 1), size(It, 2), n * R);
            Om = pageinv(A(:, :, ok(:)));
            if allIn
                e = true(size(It, 1), 1);
            else
                e = eff;
            end
            tr = zeros(1, 1, size(Om, 3));
            for i = find(e)'
                tr = tr + Om(i, i, :);
            end
            fr(ok) = tr(:) / Pe;
        end
end
if strcmp(Mm.stat, 'median')
    f = median(fr, 2);
else
    f = mean(fr, 2);
end
if numel(sz) < 3
    f = f(1);
end
end
