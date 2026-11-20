function [ok, why] = dd_rowok(C, Dr)
%DD_ROWOK Row-level feasibility of choice tasks. Internal.
%   Dr is n-by-K (level indices). ok is n-by-1. why (optional) is n-by-1 with
%   the first violated rule: 0 none, 1 reject expression, 2 require expression,
%   3 dominated alternative, 4 identical alternatives, 5 maxProb.

n = size(Dr, 1);
ok = true(n, 1);
wantWhy = nargout > 1;
if wantWhy
    why = zeros(n, 1);
end
for q = 1:numel(C.rejectFun)
    r = logical(C.rejectFun{q}(Dr, C.LV));
    if wantWhy, why(ok & r) = 1; end
    ok = ok & ~r;
end
for q = 1:numel(C.requireFun)
    r = logical(C.requireFun{q}(Dr, C.LV));
    if wantWhy, why(ok & ~r) = 2; end
    ok = ok & r;
end
for i = 1:numel(C.dom)
    d = C.dom(i);
    nq = numel(d.ka);
    dl = zeros(n, nq);
    eq = false(n, nq);
    for q = 1:nq
        la = Dr(:, d.ka(q));
        lb = Dr(:, d.kb(q));
        dl(:, q) = d.ua{q}(la) - d.ub{q}(lb);
        eq(:, q) = C.LV{d.ka(q)}(la) == C.LV{d.kb(q)}(lb);
    end
    ident = all(eq, 2);
    unknown = any(~eq(:, ~d.comp), 2);
    dc = dl(:, d.comp);
    tol = 1e-12;
    aDom = all(dc >= -tol, 2) & any(dc > tol, 2) & ~unknown;
    bDom = all(dc <= tol, 2) & any(dc < -tol, 2) & ~unknown;
    if wantWhy
        why(ok & (aDom | bDom)) = 3;
        why(ok & ident & ~(aDom | bDom)) = 4;
    end
    ok = ok & ~(aDom | bDom | ident);
end
if ~isempty(C.maxProb) && any(ok)
    Mm = C.model(1);
    X = dd_X(Mm, Dr(ok, :));
    V = pagemtimes(X, Mm.b0);
    V = V - max(V, [], 1);
    p = exp(V);
    p = p ./ sum(p, 1);
    big = reshape(max(p, [], 1), [], 1) > C.maxProb;
    idx = find(ok);
    if wantWhy, why(idx(big)) = 5; end
    ok(idx(big)) = false;
end
end
