# DCE-design: instructions for AI agents

Use this toolbox whenever you are asked to generate, check, compare or update an experimental design for a discrete choice experiment. Do not write a one-off search script: this code is tested against published designs and against designs made in Ngene, and it is reproducible.

Field reference: [SPEC.md](SPEC.md). Functions and options: [README.md](README.md). Below, `<toolbox>` is the folder that holds this file.

From the cloned repository, try `matlab -batch "dce_cli('examples/pilot.json','outputs/pilot_v1')"`. This needs a working MATLAB installation in the session that runs the command.

## Workflow

1. **Write the specification** as JSON (or a MATLAB struct) and save it in the project, next to where the design will live. The JSON is the record of how the design was made.
2. **Check it** and read the output: the table coefficient / level / base / prior, and the advice (rows and blocks against level balance, priors, constraints, size of the problem).
   ```
   matlab -batch "addpath('<toolbox>'); dce_check('spec.json')"
   ```
3. **Generate**:
   ```
   matlab -batch "addpath('<toolbox>'); dce_cli('spec.json', 'outputs/design_v1')"
   ```
   The last line printed starts with `DCE_RESULT ` followed by JSON: `ok`, `objective`, per-model `Derror`, the checks, `seconds`, `notes`, `advice` and the paths of the files written (`_design.csv`, `_long.csv`, `.xlsx`, `.json`, `.mat`, `_report.txt`). On failure `ok` is false and `error` holds an identifier (see below).
4. **Report to the user**: D-error of each model (Bayesian mean and at the prior means), predicted choice shares (is the status-quo share plausible?), level counts, block balance, the notes and the advice, and how the design compares with the one it replaces (step 5).
5. **Compare** with an existing design under the same specification:
   ```
   matlab -batch "addpath('<toolbox>'); dce_cli('spec.json', 'outputs/old_design_check', 'Evaluate', 'old_design.csv')"
   ```
   Designs are comparable only under the same priors, coding and parameters.

Within MATLAB the same steps are `dce_check(spec)`, `r = dce_design(spec)`, `dce_evaluate(spec, design)`, `dce_export(r, folder)`. `spec.type` selects a full factorial (`'full'`) or an orthogonal design (`'orthogonal'`, `'orthogonal-sequential'`) instead of an efficient one.

## Rules

- **Never edit a design by hand** (for example to remove a dominated pair). Add the rule to `reject` / `require` / `dominance` and regenerate. Hand edits are how fielded cards end up different from the analysed design.
- **Field the file the toolbox wrote.** Build survey cards from `<name>_long.csv` or the xlsx, and keep the `.json` beside them.
- **Reproduction.** Keep the input specification, any starting cards, command, seed, package and MATLAB versions, and the draw method. Use a fixed number of starts and passes for repeatable searches. A wall-clock limit can stop at different moves.
- **Priors**: take them from a previous wave or a pilot, and say where each came from; mark invented ones. Bayesian priors centred on the estimates with a standard deviation of one to two standard errors are common practice.
- **Base levels**: the base of a dummy or effects coded attribute is the first level listed unless `base` says otherwise (in Ngene scripts: the last level). Check with `dce_check`.
- **Zero priors** switch off the dominance check for that attribute and make the design "utility neutral" for it. If only the direction is known, use a small prior with that sign.
- **Status quo**: give its levels in `fixed`; a status-quo constant goes into its utility. If programme alternatives must not take status-quo levels, restrict them with `altLevels`.
- **Constants** enter the criterion by default. For a Bayesian design use the mean or median over prior draws and check `models(m).stat`.
- **Level balance** is strict by default. Under strong constraints use `frequency` with `min`/`max` instead of switching balance off; otherwise linear attributes collapse to their extreme levels.
- **Search method**: leave `Method` at `'auto'`. With exact level frequencies the search exchanges levels within columns; where frequencies have room it also replaces whole tasks from the feasible candidates. `result.search.method` and the advice say what was used and why.
- **Panel mixed logit** (`rppanel`): supply the actual block assignments when evaluating cards. The default respondent answers one block. Final D-errors integrate the priors requested in the specification. Stage 2 uses the panel prior means to screen exchanges and validates the retained design. Report `models(m).Derror` and its method. Compare designs with the same integration settings. The documented convergence tests compare two implementations inside the toolbox.
- **Several models.** For default normalized searches, compare designs using the same `search.evaluationReference` values and model weights. The evaluation-draw score is `search.evaluatedScore`; `objective` is the weighted sum of raw criteria.
- **Long runs** (Bayesian priors with many starts, large designs): `'Parallel', true` runs the starts on a parallel pool with the same result; `'MaxTime', seconds` returns the best design so far.
- **Ngene files**: `dce_design('script.ngs')` works directly; `dce_readngd('file.ngd')` returns Ngene's design and its script. The translation uses this toolbox's conventions where Ngene's differ (constants in the criterion, Bayesian priors integrated, respondents answer one block, single-model reference scales) and lists each in `notes`; read them. To reproduce a value that Ngene reported, translate with `dce_ngene(script, 'Emulate', true)`.
- **Do not overwrite deliverables**: `dce_export` refuses existing files; use a new `Name` or folder.

## Typical run times (one core of a desktop computer)

| Problem | Time |
|---|---|
| MNL, fixed priors, 18–36 rows | 1–30 s |
| MNL, fixed priors, 60 rows in 10 blocks, or 33 parameters | 1–2 min |
| Bayesian MNL (200 draws), 48 rows, 2–3 alternatives | 1–2 min |
| MNL + panel mixed logit, 20–36 rows in 3–5 blocks | 0.3–3 min |

## Error identifiers

| Identifier | Meaning | What to do |
|---|---|---|
| `dce:spec:unknownField` | misspelt or unsupported field | fix the name (the message lists the known ones) |
| `dce:spec:priors`, `dce:spec:random` | wrong number or form of priors | one entry per non-base level; see SPEC.md |
| `dce:spec:utility` | utility string not understood | `par[.dummy] * attr * attr + ...` |
| `dce:spec:constraint` | expression refers to an unknown name or has bad syntax | use `alternative.attribute` |
| `dce:spec:rows` | fewer observations than parameters | more rows or fewer parameters |
| `dce:spec:frequency` | level counts cannot be met | adjust `frequency` or `rows` |
| `dce:search:infeasible` | constraints (with level balance) leave too few tasks | the message says which rule excludes what; relax it or use `frequency` bounds |
| `dce:search:notIdentified` | a parameter has no variation between alternatives | drop it or change the utility |
| `dce:search:exhaustive` | the problem is too large for `Method = 'exhaustive'` | use the default method |
| `dce:design:blocks` | missing or invalid block assignments | supply finite numeric block labels for each task |
| `dce:design:levels`, `dce:design:columns` | a design to evaluate does not match the specification | check column names and level values |
| `dce:ngene:unsupported` | Ngene feature not supported | write the specification directly |
| `dce:export:exists` | output file exists | new name or folder |

## Checking the installation

```
matlab -batch "addpath('<toolbox>'); dce_selftest"
```
All checks must pass (about a minute). Run it after any change to the toolbox.
