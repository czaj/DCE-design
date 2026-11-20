function Cs = dd_searchmodel(C, o)
%DD_SEARCHMODEL Version of the compiled specification used inside the search. Internal.
%   The search needs models whose information adds up over choice tasks.
%   mnl and rp models are used as they are. A panel mixed logit (rppanel) is
%   replaced by a surrogate:
%     'mnl'  the multinomial logit at the prior means of the coefficients
%            (criterion over the means only), or
%     'rp'   the cross-sectional mixed logit with the same parameters.
%   The panel model itself is used when the candidate designs are evaluated.

Cs = C;
for m = 1:C.M
    Mm = C.model(m);
    if ~strcmp(Mm.type, 'rppanel')
        continue
    end
    if strcmp(o.PanelSurrogate, 'rp')
        Mm.type = 'rp';
    else
        Mm.type = 'mnl';
        Mm.eff = Mm.eff(1:Mm.P);
        if ~any(Mm.eff)
            Mm.eff(:) = true;
        end
        Mm.Peff = nnz(Mm.eff);
        Mm.Pt = Mm.P;
        Mm.thNames = Mm.parNames;
        Mm.thBase = Mm.parBase;
        Mm.ridx = zeros(1, 0);
        Mm.Pr = 0;
        Mm.rdist(:) = 0;
        Mm.sdist(:) = 0;
        Mm.isBayes = any(Mm.bdist > 0);
        if ~Mm.isBayes
            Mm.stat = 'fixed';
        end
    end
    Cs.model(m) = Mm;
end
end
