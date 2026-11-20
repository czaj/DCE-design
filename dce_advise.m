function adv = dce_advise(x)
%DCE_ADVISE Remarks on a specification or on a finished design.
%
%   DCE_ADVISE(spec) prints what the toolbox can say about a specification
%   before any search: how many choice tasks and designs are possible and how
%   long evaluating all of them would take (and whether that is an option),
%   whether the number of rows and blocks allows level balance in the design
%   and inside every block and which numbers would, how many independent
%   observations there are per parameter, priors that are all zero or that
%   leave the sign open, and how much the constraints exclude.
%
%   DCE_ADVISE(result), with a result of dce_design or dce_evaluate, adds what
%   depends on the design: the weakest parameter and the sample size it asks
%   for, tasks with a near-certain choice, an alternative that takes most
%   choices, a lower bound on the D-error of any design for this problem and
%   how close this design is to it, which search method was used and why, and
%   whether the random starts agree.
%
%   adv = DCE_ADVISE(...) returns a table with the variables topic, level
%   ('info', 'suggestion' or 'warning') and text, and prints nothing.
%
%   Every remark rests on a number computed for this problem. The same remarks
%   are part of the output of dce_check, dce_design and dce_evaluate.
%
%   See also dce_check, dce_design, dce_evaluate.

if isstruct(x) && isfield(x, 'design') && isfield(x, 'advice')
    a = x.advice;
else
    a = dd_advice(dd_compile(dd_readspec(x)));
end
if isempty(a)
    T = table(cell(0, 1), cell(0, 1), cell(0, 1), 'VariableNames', {'topic', 'level', 'text'});
else
    T = struct2table(a, 'AsArray', true);
end
if nargout > 0
    adv = T;
    return
end
if isempty(a)
    fprintf('No remarks.\n');
end
for i = 1:numel(a)
    fprintf('%-12s %-11s %s\n', a(i).topic, a(i).level, a(i).text);
end
end
