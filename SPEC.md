# Design specification

A specification is a MATLAB struct or a JSON file with the same fields. `dce_check(spec)` validates it and prints how it is read; `dce_design(spec)` searches for the design; `dce_evaluate(spec, design)` evaluates an existing one. An Ngene script can be used instead (see [Ngene scripts](#ngene-scripts)).

Unknown fields are errors, so a typo cannot be ignored silently.

## Minimal example

```matlab
spec.name         = 'water_main';
spec.alternatives = {'sq', 'a', 'b'};
spec.rows         = 24;
spec.blocks       = 3;

spec.attributes.quality = [0 1 2];           % level values, in any order
spec.attributes.access  = [0 1];
spec.attributes.cost    = [0 20 40 80 120];

spec.altLevels.a.cost = [20 40 80 120];      % the programmes never cost 0
spec.altLevels.b.cost = [20 40 80 120];
spec.fixed.sq = struct('quality', 0, 'access', 0, 'cost', 0);   % status quo

u = 'b_q.dummy*quality + b_acc*access + b_cost*cost';
spec.model.type       = 'mnl';
spec.model.utility.sq = ['asc_sq + ' u];
spec.model.utility.a  = u;
spec.model.utility.b  = u;
spec.model.priors.asc_sq = 'N(0.2, 0.25)';
spec.model.priors.b_q    = {'N(0.35, 0.15)', 'N(0.60, 0.15)'};   % levels 1 and 2; base 0
spec.model.priors.b_acc  = 0.3;
spec.model.priors.b_cost = 'N(-0.012, 0.004)';

spec.reject = {'a.quality == 0 & a.access == 0', 'b.quality == 0 & b.access == 0'};
```

The same as JSON:

```json
{
  "name": "water_main",
  "alternatives": ["sq", "a", "b"],
  "rows": 24,
  "blocks": 3,
  "attributes": {"quality": [0, 1, 2], "access": [0, 1], "cost": [0, 20, 40, 80, 120]},
  "altLevels": {"a": {"cost": [20, 40, 80, 120]}, "b": {"cost": [20, 40, 80, 120]}},
  "fixed": {"sq": {"quality": 0, "access": 0, "cost": 0}},
  "model": {
    "type": "mnl",
    "utility": {
      "sq": "asc_sq + b_q.dummy*quality + b_acc*access + b_cost*cost",
      "a": "b_q.dummy*quality + b_acc*access + b_cost*cost",
      "b": "b_q.dummy*quality + b_acc*access + b_cost*cost"
    },
    "priors": {"asc_sq": "N(0.2, 0.25)", "b_q": ["N(0.35, 0.15)", "N(0.60, 0.15)"], "b_acc": 0.3, "b_cost": "N(-0.012, 0.004)"}
  },
  "reject": ["a.quality == 0 & a.access == 0", "b.quality == 0 & b.access == 0"]
}
```

## Fields

| Field | Required | Meaning |
|---|---|---|
| `name` | no | Used in reports and file names. |
| `type` | no | `'efficient'` (default), `'full'` (all feasible choice tasks; with `rows`, a random fraction of that size), `'orthogonal'` or `'orthogonal-sequential'`. See [Other design types](#other-design-types). |
| `alternatives` | yes | Names of the alternatives (valid identifiers). |
| `rows` | yes | Number of choice tasks in the design. Not needed for `type = 'full'`; a lower bound for orthogonal designs. |
| `blocks` | no | Number of blocks; each respondent answers one. Default 1. |
| `attributes` | yes | `name -> vector of levels` (numbers). An attribute with one level is a constant. |
| `altLevels` | no | `alternative -> attribute -> levels`, when an alternative uses other levels than `attributes`. |
| `fixed` | no | `alternative -> attribute -> value`: the attribute is constant in that alternative (status quo, opt-out). |
| `scenario` | no | Attributes that describe the choice context: one design column shared by all alternatives. |
| `base` | no | `attribute -> level`: base level for dummy and effects coding. Default: the **first** level listed. |
| `labels` | no | `attribute -> list of labels`, one per level; used in the exported labelled sheet. |
| `model` / `models` | yes | One model, or a list of models to average over (see below). |
| `weights` | no | One weight per model. Default 1 each. |
| `criterion` | no | `'d'` (D-error, default) or `'a'` (A-error). |
| `constants` | no | `true` (default): constants (parameters that multiply no attribute, such as a status-quo constant) are part of the criterion. `false`: they stay in the model but are left out of the criterion. |
| `exclude` | no | Other parameters left out of the criterion (they are still in the model). |
| `excludeMode` | no | `'nuisance'` (default): excluded parameters are estimated, the criterion uses the covariance block of the others. `'known'`: they are treated as known. |
| `reject` | no | Expressions; a task is infeasible if any is true. |
| `require` | no | Expressions; a task is infeasible unless all are true. |
| `dominance` | no | `'auto'` (default), `true`, `false`, or a list of alternatives. See [Dominance](#dominance-identical-alternatives-repeated-tasks). |
| `maxProb` | no | Reject tasks in which one alternative's choice probability at the prior means exceeds this. |
| `balance` | no | `'strict'` (default): every level of a column appears equally often (as equally as the number of rows allows). `'none'`: no restriction. |
| `frequency` | no | `attribute -> struct` with `exact`, or `min` and/or `max` (one number, or one per level): overrides `balance` for that attribute. |
| `unique` | no | `true` (default): no choice task may be repeated. |
| `mirror` | no | `true` (default): a task that equals another after exchanging two interchangeable alternatives counts as a repeat. |
| `panelUnit` | no | `'block'` (default): in panel models a respondent answers one block. `'design'`: every row (Ngene's convention). |
| `start` | no | Starting design: a table, a matrix of level values, a result struct, or a `.csv` / `.xlsx` / `.ngd` file. |
| `search` | no | Options of `dce_design` (see README), e.g. `{"Seed": 3, "Restarts": 8}`. |
| `notes` | no | Free text carried into the result. |

## Models

```matlab
spec.model.name    = 'MNL';          % optional
spec.model.type    = 'mnl';          % 'mnl', 'rp' (cross-sectional mixed logit), 'rppanel' (panel mixed logit)
spec.model.utility = struct(...);    % alternative -> utility string
spec.model.priors  = struct(...);    % parameter -> prior
spec.model.random  = struct(...);    % parameter -> distribution (rp, rppanel)
spec.model.stat    = 'mean';         % how Bayesian prior draws are summarised: 'mean' (default), 'median', 'fixed'
```

### Utility strings

A utility is a sum of terms `parameter[.coding] * attribute * attribute ...`.

- A parameter with no attribute is a constant: `asc_sq`.
- The same parameter name in several alternatives makes it generic; different names make it alternative-specific.
- `.dummy` or `.effects` after the parameter codes the **first** attribute of the term: one coefficient per non-base level, named `par[level]`. Write the suffix wherever the parameter appears.
- More than one attribute is an interaction: `b_x*price*purpose`. With `.dummy`, the first attribute is coded and the others multiply it.
- A quadratic term is `b_q2*price*price`.
- Negative effects are written as negative priors; there is no minus sign in a utility string.
- An alternative with no utility (omitted or `''`) has utility zero.
- A prior may be written inline, `b_cost[-0.02]*cost` or `b_q.dummy[0.3|0.6]*quality`, instead of in `priors`.

### Priors

| Form | Meaning |
|---|---|
| `0.3` | fixed prior |
| `'N(0.3, 0.1)'` or `'normal(0.3, 0.1)'` | Bayesian prior, normal with mean 0.3 and standard deviation 0.1 |
| `'U(0.1, 0.5)'` or `'uniform(0.1, 0.5)'` | Bayesian prior, uniform |
| `[0.3 0.6]`, `{'N(0.3,0.1)', 0.6}`, `'0.3|0.6'` | one entry per non-base level of a dummy or effects coded attribute, **in the order of the levels as listed** |

A parameter without a prior gets 0, and the result says so. A design optimised with zero priors treats the attribute as irrelevant to choice; a dominance check cannot rank its levels.

`dce_check` prints the table coefficient / level / base / prior. Read it once for every new specification: a prior attached to the wrong level is the most common error in design work.

### Random coefficients (`rp`, `rppanel`)

`random` lists the random coefficients; the others stay fixed.

| Form | Coefficient |
|---|---|
| `'normal(0.4)'` | mean (from `priors`) + 0.4 × standard normal |
| `'lognormal(0.5)'`, `'-lognormal(0.5)'` | ± exp(m + 0.5 × standard normal), where the prior m is the mean of the underlying normal |
| `'uniform(0.3)'` | mean ± 0.3, uniform |
| `'triangular(0.3)'` | mean ± 0.3, symmetric triangular |
| `'normal(U(0.2, 0.6))'` | Bayesian prior on the spread |
| `'fixed'` | not random (useful for single levels of a coded attribute) |

The model's parameters are the means and the spreads; the D-error covers both.

**Panel mixed logit and blocks.** In an `rppanel` model the choices of one respondent share that respondent's coefficients, so the information depends on which tasks one person answers. With `blocks > 1` a respondent answers one block, and the information of the design is the sum over blocks of the information of one respondent per block (`panelUnit = 'block'`, the default). The information of a block is a sum over all possible choice sequences when their number is at most `PanelExactMax`, otherwise an average over simulated respondents (`PanelRespondents` in the search, `EvalPanelRespondents` for the reported values). The integral over the random coefficients is computed for each sequence by importance sampling around the mode of that sequence (`PanelDraws`, `EvalPanelDraws` draws per sequence). `PanelIntegration = 'prior'` uses the same `EvalRandomDraws` draws for all sequences instead; the documented long-panel test cases gave lower D-errors with small samples of shared draws. `dce_design` works in two stages: a search with a surrogate that adds up over tasks (the multinomial logit at the model's prior means), then a search using the panel model at its prior means, by exchanges within blocks (with simulated sequences, an exchange found on `PanelRespondents` respondents is confirmed on `EvalPanelRespondents` others before it is applied). The second stage's design is kept only if it is also better on a validation sample (other simulated respondents and draws than the search used). The choice among the best designs of stage 1 is made on the same validation sample, and the reported values come from a third sample, so that selecting a design does not flatter its reported D-error.

### Several models

```matlab
m1 = struct('name', 'MNL', 'type', 'mnl',     'utility', U, 'priors', P1);
m2 = struct('name', 'MXL', 'type', 'rppanel', 'utility', U, 'priors', P2, 'random', R2);
spec.models  = {m1, m2};
spec.weights = [1 1];
```

For multiple active models, the default search divides each criterion by the corresponding single-model reference. Stage 1 uses surrogate references. Candidate selection, panel polishing and validation use the same per-model scales stored in `search.evaluationReference`. To compare an existing design, apply those same scales to its evaluated model criteria: `sum(w(active).*criterion(active)./scale(active))`. With `Normalize=false`, each scale is one.

## Constraints

Expressions use MATLAB operators (`& | ~ == ~= < <= > >= + - * /`) and refer to attributes as `alternative.attribute` (scenario variables by name alone). Fixed levels can be referred to in the same way.

```matlab
spec.reject  = {'a.cost > 0 & a.quality == 0 & a.access == 0', ...     % a costly programme must improve something
                'a.quality >= b.quality & a.access >= b.access & a.cost <= b.cost'};
spec.require = {'abs(a.cost - b.cost) <= 80'};
```

Constraints and level frequencies are both hard. If they cannot be met together, `dce_design` stops with an error that says which rule excludes how many tasks, and how many distinct feasible tasks exist.

`dce_evaluate` checks observed counts against the supplied frequency rules, even if an existing design violates them. For a blocked panel, pass the actual block column or `Blocks` vector.

### Dominance, identical alternatives, repeated tasks

With `dominance = 'auto'`, alternatives with the same utility structure and the same attributes are compared: a task is excluded if one of them is at least as good on every attribute and better on one, or if two are identical. The ranking of levels comes from the prior means of the first model (main effects only). An attribute with zero priors cannot be ranked; tasks that differ on it are never treated as dominated. To state only a direction, use a small prior with the right sign.

## Search methods

`dce_design` chooses the method from the level-frequency rules (`Method = 'auto'`), and `result.search.method` with `methodReason` says what it chose.

| Method | Moves | Used by `'auto'` when |
|---|---|---|
| `'swap'` | Two tasks exchange their levels in one column; where frequency bounds leave room, one task changes a level. | Every level frequency is fixed exactly (strict balance with a number of rows that divides evenly, or `frequency.exact`). Nothing else can change such a design without breaking the frequencies. |
| `'both'` | The moves of `'swap'`, in turn with row exchange: a whole task is replaced by the best task from the feasible candidates (a modified Fedorov exchange). | Some level frequency has room (`balance = 'none'`, `frequency` with `min`/`max`, or a number of rows that does not divide evenly). Row exchange changes several attributes of a task at once, which matters when constraints tie attributes together. |
| `'rowexchange'` | Row exchange only. | Never chosen automatically. |
| `'exhaustive'` | Every set of `rows` distinct feasible tasks that keeps the level frequencies is evaluated. | Never chosen automatically; `dce_advise` says when a problem is small enough (at most `ExhaustiveMax` designs, 2 million by default). Exhaustiveness applies to the enumerated stage-1 criterion and its search draws. A panel surrogate or a fresh evaluation sample can give a different ranking. |

Constraints (`reject`, `require`, dominance, no repeated tasks) are respected by every method: a move that would break one is not considered. The candidates for row exchange are all feasible tasks when there are at most `CandidateMax` (3000), otherwise a random sample of that size.

With Bayesian priors, moves are first ranked on a few prior draws and the best are confirmed on all draws. Panel mixed logit models are searched in two stages (see [Models](#models)).

## Other design types

```matlab
spec.type = 'full';                    r = dce_design(spec);      % or dce_full(spec)
spec.type = 'orthogonal';              r = dce_design(spec);      % or dce_orthogonal(spec)
spec.type = 'orthogonal-sequential';   r = dce_design(spec);      % or dce_orthogonal(spec, 'Type', 'sequential')
```

**Full factorial.** All combinations of the levels of all design columns, without the tasks that break the constraints and without repeats of a task with interchangeable alternatives exchanged. `rows` (or `dce_full(spec, 'Rows', n)`) draws a random fraction. `dce_full(spec, 'Constraints', false)` returns the raw full factorial.

**Orthogonal.** The design columns form an orthogonal array of strength 2: every level equally often, every pair of levels of two columns equally often, so all columns are uncorrelated. Priors are not used. *Simultaneous* puts all columns of all alternatives in one array; *sequential* builds an array for one alternative and gives the others the same array in another row order (orthogonal within alternatives only, with far fewer rows). An array needs a number of rows that is a multiple of Li·Lj for every pair of columns; `rows` is raised to the first such number for which an array is found, and the result says so. With `blocks`, the block is tried as one more column of the array, which balances every block perfectly. Orthogonal designs cannot respect constraints on tasks in general; `result.checks` lists the tasks that break them.

All three types return the same result structure as an efficient design, evaluated under the priors of the specification, so their D-errors can be compared directly.

## Results

`dce_design` and `dce_evaluate` return one struct; the main fields:

| Field | Content |
|---|---|
| `design` | table: `Task`, `Block`, one variable per design column (`alt_attribute`), level values, sorted by block |
| `columns`, `levels`, `blocks` | column names, level indices, block of each row |
| `objective` | weighted sum of the raw model criteria, for both generation and evaluation |
| `search.evaluatedScore`, `search.evaluationReference` | search score on reporting draws and its per-model scales (efficient designs) |
| `models(m)` | `Derror`, `Aerror` (over the prior draws), `DerrorFixed` (at the prior means), `bayes` (mean, median, quantiles, share of draws not identified), `parameters`, `prior`, `information`, `avc`, `se`, `tRatio`, `Sestimate` (all for one replication of the design), `method` |
| `probabilities`, `meanShare`, `utilityBalance` | predicted choice shares at the prior means of model 1 |
| `levelCounts`, `correlations` | per design column |
| `checks` | `constraintsOK`, `violations` (task, rule), `duplicateTasks`, `frequenciesOK`, `identified` |
| `blocking` | block balance Psi, absent levels, extremes, leave-one-block-out D-error |
| `search` | options, the method used and why, starts, time, objective of every start, candidates considered; for full and orthogonal designs how the design was constructed |
| `advice` | remarks on the specification and the design, each with a topic, a level (`info`, `suggestion`, `warning`) and a text (see `dce_advise`) |
| `notes` | assumptions and warnings |
| `spec`, `specHash`, `seed`, `toolbox` | design specification and provenance; retain any starting file, command-line overrides and the draw method |

Conventions: the D-error is `det(Omega)^(1/K)` where `Omega` is the covariance matrix for **one respondent who answers every row** (or, with blocks, one respondent per block) and `K` the number of parameters in the criterion. With Bayesian priors it is the mean (or median) of the D-errors over the prior draws. S-estimates count replications of the whole design; multiply by the number of blocks to get respondents.

## Ngene scripts

`dce_ngene(script)` translates an Ngene script; `dce_design('file.ngs')` and `dce_evaluate('file.ngs', design)` do so implicitly. Names are matched without regard to case, as in Ngene.

Taken over as written:

| Ngene | Here |
|---|---|
| base level = last level listed | `base` set to that level |
| `*` after an alternative | `dominance` lists those alternatives |
| `;cond: if(a, b)` | `reject: (a) & ~(b)` |
| `;reject`, `;require` | `reject`, `require` |
| `x[1,2,3](4,4,4)`, `(2-6,2-6)` | `frequency` |
| `x[x]` in later utilities | `scenario` |
| an attribute with one level, e.g. `cost_sq[0]` | a constant |
| `;start = file` | `start` |

Where Ngene's conventions differ from the toolbox, the translation uses the toolbox's and writes a note for each case. `dce_ngene(script, 'Emulate', true)` keeps Ngene's instead, which is how a value reported by Ngene is reproduced.

| | Ngene | Default translation | `'Emulate', true` |
|---|---|---|---|
| Constants | left out of the efficiency measure unless `;con` | part of the criterion | `constants = false` |
| Bayesian priors | used at their means unless `;eff` says `mean` or `median`; a script with `(n,…)` priors and `;eff = (mnl,d)` is therefore not optimised over the prior distribution | D-error averaged over the prior distribution | `stat = 'fixed'` |
| Panel mixed logit | a respondent answers every row | a respondent answers one block | `panelUnit = 'design'` |
| Several models | the weights of `;eff` multiply raw D-errors | single-model reference scales, equal weights | raw weights, `Normalize = false` |

In both modes:

- An attribute that a `;cond` restricts to some of its levels (`if(a.x = 1, a.z = [0,1])`) cannot stay exactly balanced. As in Ngene it loses exact balance; here its level counts may vary between half and one and a half times the balanced count.
- The numbers of draws in `;bdraws` and `;rdraws` set the search draws, between 128 and 256; the reported values always use the toolbox's evaluation draws. `;rep` sets the simulated respondents for the reported values of panel models, at least 2000.
- The search settings of `;alg` are not used, except a `stop = total(...)` time limit.

- `;fact` gives `type = 'full'` (with `;rows = n`, a random fraction), `;orth = sim` gives `type = 'orthogonal'`, `;orth = seq` gives `type = 'orthogonal-sequential'`.

Unsupported features (error components, `;fisher`, `;orth = ood`, pivot and covariate modifiers, continuous levels, sign-only priors) raise an error that names them. So does an attribute that is used but never given levels. Five archived files contained a misspelt attribute name in this position.
