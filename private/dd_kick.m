function D = dd_kick(C, D, rs, n)
%DD_KICK Perturb a design by n random feasible swaps within columns. Internal.
S = C.S;
K = C.K;
done = 0;
tries = 0;
while done < n && tries < 200 * n
    tries = tries + 1;
    k = randi(rs, K);
    ij = randperm(rs, S, 2);
    if D(ij(1), k) == D(ij(2), k)
        continue
    end
    Dn = D;
    Dn(ij, k) = D(fliplr(ij), k);
    if all(dd_rowok(C, Dn(ij, :))) && ~any(dd_dups(C, Dn))
        D = Dn;
        done = done + 1;
    end
end
end
