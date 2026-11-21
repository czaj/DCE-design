# DCE-design

Experimental designs for discrete choice experiments in MATLAB: generate a design, block it, evaluate it, compare it with an existing one.

The toolbox finds efficient designs: it minimises the D-error (or A-error) of a multinomial logit, a cross-sectional or a panel mixed logit model, or of a weighted combination of models, with fixed or Bayesian priors, under level-balance rules and constraints on the choice tasks. It also builds orthogonal and full factorial designs, reads Ngene scripts, and says what it notices about a specification and a design. Each design comes with its evaluation and checks. To reproduce the cards, keep the specification, starting design, seed, exact call and software versions. Use the same draw method and a fixed search budget in starts and passes. A wall-clock cutoff can stop a run at a different move.

## Getting started

Clone the public repository, change into it and run the bundled JSON example from a shell with MATLAB on the path.

```sh
git clone https://github.com/czaj/DCE-design.git
cd DCE-design
matlab -batch "dce_check('examples/pilot.json'); dce_cli('examples/pilot.json','outputs/pilot_v1')"
```

The final `DCE_RESULT` line contains JSON. The wide CSV records one row per task and the long CSV one per alternative. For a second run, choose a new output folder. To evaluate the exported cards under the same specification, run:

```sh
matlab -batch "dce_cli('examples/pilot.json','outputs/pilot_check','Evaluate','outputs/pilot_v1/pilot_design.csv')"
```


```matlab
addpath('...\DCE-design')

spec.alternatives = {'a', 'b', 'sq'};
spec.rows   = 24;
spec.blocks = 3;
spec.attributes.quality = [0 1 2];
spec.attributes.access  = [0 1];
spec.attributes.cost    = [0 20 40 80 120];
spec.altLevels.a.cost = [20 40 80 120];
spec.altLevels.b.cost = [20 40 80 120];
spec.fixed.sq = struct('quality', 0, 'access', 0, 'cost', 0);

u = 'b_q.dummy*quality + b_acc*access + b_cost*cost';
spec.model.utility.a  = u;
spec.model.utility.b  = u;
spec.model.utility.sq = ['asc_sq + ' u];
spec.model.priors.b_q    = {'N(0.35, 0.15)', 'N(0.60, 0.15)'};   % Bayesian priors, levels 1 and 2
spec.model.priors.b_acc  = 0.3;                                   % fixed prior
spec.model.priors.b_cost = -0.012;
spec.model.priors.asc_sq = 0.2;

dce_check(spec);                 % how the specification is read, and advice on it
r = dce_design(spec);            % search, block, evaluate, print a summary
dce_export(r, 'outputs');        % csv (wide and long), xlsx, json, mat, report
```

`examples/` has five runnable scripts: a basic design, a Bayesian design with a status quo and constraints, one design for a multinomial logit and a panel mixed logit, an Ngene script, and full factorial, orthogonal and efficient designs side by side. [SPEC.md](SPEC.md) documents every field of the specification. [AGENTS.md](AGENTS.md) is a short guide for AI agents that run the toolbox from the command line.

## Functions

| Function | Purpose |
|---|---|
| `dce_check(spec)` | Validate a specification and print how it is read, with advice. |
| `dce_design(spec, ...)` | Search for an efficient design, assign blocks, evaluate. With `spec.type`, a full factorial or an orthogonal design. |
| `dce_evaluate(spec, design, ...)` | Evaluate an existing design (matrix, table, `.csv`, `.xlsx`, `.ngd`). |
| `dce_advise(spec or result)` | Remarks on a specification or a finished design. |
| `dce_full(spec, ...)` | Full factorial (all feasible choice tasks) or a random fraction. |
| `dce_orthogonal(spec, ...)` | Orthogonal design, simultaneous or sequential. |
| `dce_block(X, B, ...)` | Balanced blocking of any design (used by the functions above). |
| `dce_report(result)` | Text summary. |
| `dce_export(result, folder)` | Write the design and its evaluation to files; never overwrites. |
| `dce_ngene(script)` | Translate an Ngene script into a specification. |
| `dce_readngd(file)` | Read a design saved by Ngene, with its script. |
| `dce_cli(specFile, folder)` | One call for scripts and AI agents: JSON or `.ngs` in, files and one line of JSON out. |
| `dce_selftest` | Test suite. |

## What the toolbox does

**Models.** Multinomial logit (`mnl`), cross-sectional mixed logit (`rp`) and panel mixed logit (`rppanel`) with normal, lognormal, uniform or triangular coefficients. Generic and alternative-specific parameters, constants, linear, dummy and effects coding, interactions, scenario variables, alternatives with fixed levels (status quo) or zero utility (opt-out).

**Priors.** Fixed, or Bayesian (normal, uniform). The Bayesian criterion is the mean (or median) of the D-error over scrambled Sobol draws of the priors; the search uses 128 draws and the reported value 2048.

**Several models.** The search scales each active model by a reference design found in a single-model search. With a panel model, the first stage uses its surrogate and the final selection scores its reference design on the panel model. The reference values and evaluated search score are returned alongside the raw weighted sum.

**Constraints.** Any logical expression on the attributes of a task (`reject`, `require`); tasks with a dominated or identical alternative; repeated tasks, including mirror images; exact, minimum and maximum level frequencies. All hold at the same time. If they cannot, the error says which rule excludes what.

**Blocks.** Every level in every block where the block size allows, the lowest and highest level of linear attributes (cost) in every block, then the best level balance (Bliemer, Hess, Collins & Amaris 2026); the design's information without each block is reported.

**Panel mixed logit.** A respondent answers one block, so the panel information is computed per block: over all choice sequences of a block, or over simulated respondents when the sequences are too many to list. The probability of a sequence is an integral over the random coefficients, and the toolbox computes it separately for every sequence: by importance sampling around the coefficients that the sequence makes most likely, with control variates that leave only the departure from normality to be simulated, and with corrections for the bias of the simulation (`private/dd_panelis.m`). In the long-panel examples tested here, sharing a small set of prior draws across sequences underestimated D-error. The magnitude depends on the design and the draw count, as shown under [How it was checked](#how-it-was-checked). The search has two stages. Stage 1 uses a surrogate that adds up over tasks (the multinomial logit at the model's prior means), blocks the best designs and compares them on the panel criterion. Stage 2 optimises the objective with the panel model itself by exchanging levels within blocks, which leaves the level balance of every block untouched. Where the choice sequences are simulated, an exchange that helps on the respondents of the search is applied only if it also helps on a larger sample. The result of the stage is kept only if it is also better on a validation sample of other respondents and draws. All choices between designs of a panel model are made on that sample, and the reported values come from a third one, so that they are not flattered by the selection. Where the search may replace whole tasks (see Search), stage 1 is run both with and without that move, and the designs of both runs are compared on the panel criterion.

**Search.** From each random feasible start the design is improved until no move helps, then perturbed and improved again; several starts run independently (in parallel if a pool is open) and give the same result either way. The information matrix is updated only from the rows that change, for all prior draws at once. The moves are chosen from the specification (`Method = 'auto'`):

- When every level frequency is fixed exactly, a design can only change by two tasks exchanging their levels in one column; all such exchanges are evaluated.
- Where level frequencies have room (bounds, no balance, or a number of rows that does not divide evenly), single levels can also change, and whole tasks are replaced by the best of the feasible candidate tasks (a modified Fedorov exchange). Replacing a whole task changes several attributes at once, which is what gets a search past constraints that tie attributes together.
- A problem small enough to list every design can be solved exactly with `Method = 'exhaustive'`.

Constraints are respected by every method. With Bayesian priors, moves are first ranked on a few draws and the best are confirmed on all. The result records which method was used and why.

**Other design types.** `dce_full` lists all feasible choice tasks or draws a random fraction. `dce_orthogonal` builds orthogonal arrays of strength 2 with mixed levels (the J2 algorithm of Xu 2002), for all alternatives at once or sequentially, with the blocks as one more orthogonal column where possible. Both are evaluated like any design, so an orthogonal and an efficient design can be compared under the same priors.

**Advice.** `dce_check`, `dce_design` and `dce_evaluate` report what the numbers of the problem imply (`dce_advise`): how many choice tasks and designs are possible and how long evaluating all of them would take, whether the rows and blocks allow level balance and which numbers would, the observations per parameter, priors that are all zero or leave a sign open, how much the constraints exclude, the weakest parameter and the sample size it asks for, tasks with a near-certain choice, whether the random starts agree, and a lower bound on the D-error of any design for the problem (from the optimal approximate design) with the design's efficiency against it.

**Evaluation.** D-error and A-error at the prior means and over the prior draws (with the spread and the share of draws for which the model is not identified), information and covariance matrices, standard errors and the number of design replications needed for t = 1.96 per parameter, predicted choice shares of every task, level counts, correlations, and all checks.

## Options of `dce_design`

Given as name-value pairs or in `spec.search`.

| Option | Default | Meaning |
|---|---|---|
| `Seed` | 1 | fixes the starts and all draws |
| `Restarts` | 12; 4 for costly problems | independent random starts (costly: Bayesian priors, `rp` models, 25 or more parameters, rows × columns of 600 or more) |
| `Kicks` | 3; 1 for costly problems | perturb-and-reoptimise rounds per start |
| `Method` | `'auto'` | `'swap'`, `'rowexchange'`, `'both'` or `'exhaustive'`; see Search above and SPEC.md |
| `CandidateMax`, `ExhaustiveMax` | 3000, 2e6 | candidate tasks for row exchange; largest number of designs an exhaustive search will evaluate |
| `MaxTime` | Inf | seconds; the best design so far is then returned (the result says so) |
| `Parallel` | `'auto'` | `'auto'` uses an open parallel pool, `true` opens one, `false` never |
| `Start` | — | starting design: level values, a table, a result, or a file |
| `BayesDraws`, `EvalBayesDraws` | 128, 2048 | prior draws in the search and for the reported values |
| `RandomDraws`, `EvalRandomDraws` | 128, 1000 | random-coefficient draws of cross-sectional mixed logit models (`rp`) |
| `ScreenDraws`, `ScreenTop`, `FullCheck` | 8, 4, false | Bayesian search: draws used to rank moves, moves confirmed on all draws, and an optional final pass that evaluates every move on all draws |
| `PanelSurrogate` | `'mnl'` | how a panel model enters stage 1 of the search: multinomial logit at the prior means, or `'rp'` |
| `PanelPolish` | 3 | passes of stage 2, on the panel criterion itself (0 switches it off) |
| `PanelDraws`, `EvalPanelDraws` | 32, 128 | panel models: draws per choice sequence in stage 2 of the search and for the reported values |
| `PanelRespondents`, `EvalPanelRespondents` | 500, 2000 | panel models, when a respondent's choice sequences are more than `PanelExactMax` (5000) and are simulated: respondents on which the moves of stage 2 are screened, and respondents on which they are confirmed and from which the values are reported (separate samples) |
| `PanelBayesDraws` | 16 | panel models: prior draws for the reported values |
| `PanelIntegration` | `'importance'` | `'prior'` is the usual simulation, with `EvalRandomDraws` draws shared by all sequences |
| `Normalize` | on with several models | scale each model by its single-model reference design |
| `Candidates` | 3 | best distinct designs that are fully evaluated before one is chosen |
| `Block`, `BlockRestarts` | true, 30 | blocking |
| `Checkpoint` | — | file stem; the result of every start is saved |
| `Verbose` | 1 | 0 silent, 1 summary, 2 progress |

## How it was checked

`dce_selftest` (under a minute) must pass after any change:

- **Published designs.** The four designs of Bliemer, Rose & Chorus (2017) give their D-errors (0.304, 0.076, 0.057, 0.064) and their numbers of tasks with a dominated alternative (8, 5, 2, 0 of 8). Three designs from Hole's `dcreate` give their D-efficiencies to ten digits, with linear and dummy coding. The design of Zwerina, Huber & Kuhfeld gives 27^(−1/2) under effects coding.
- **Information matrices.** For the multinomial logit, the cross-sectional and the panel mixed logit (normal, lognormal, uniform and triangular coefficients), the information matrix equals minus the numerical second derivative of the expected log-likelihood (relative error below 2e-6). The panel information from simulated respondents agrees with the one from all choice sequences.
- **Panel information.** The integration over the random coefficients was checked against plain simulation with so many draws that it has converged. For a design with 18 tasks per respondent and six normal coefficients, the scores of 200 choice sequences from 12.6 million draws give a D-error of 6.01. Importance sampling gives 6.06 and 6.08 with 128 draws per sequence and 6.01 and 6.00 with 1,024 (two seeds each). Plain simulation with one set of draws for all sequences gives 1.9–2.0 with 250 draws, 3.2–3.4 with 1,000, 4.6–4.7 with 4,000, 5.6–5.8 with 16,000 and 5.8–6.0 with 64,000: it approaches the same value from below, and with the few hundred or thousand draws that are usual it reports a D-error that is far too low. In three of our designs with shorter sequences, 1,000 such draws understated the D-error by 1% (4 tasks per respondent), 9% (8 tasks) and 14% (12 tasks). The self-test repeats the comparison on a design of 6 tasks with normal, lognormal, uniform and triangular coefficients (4.867 against 4.856 from 131,072 draws) and on one of 16 tasks (1.2880 with 32 draws per sequence, 1.2882 with 512, 1.2873 from plain simulation with 65,536 draws, 1.199 with 250).
- **MNL standard errors.** Choices are simulated and the multinomial logit is re-estimated 300 times. The spread of the estimates matches the predicted standard errors (ratios 0.97–1.06).
- **Panel mixed logit standard errors.** An independent check fits 600 synthetic samples of 1,500 respondents, with five choices per person and a normally distributed coefficient. The ratios of the empirical standard deviations to the toolbox predictions were 1.002 for the constant, 1.027 for the mean and 0.965 for the spread. To reproduce this optional check, run `python tests/validate_panel_se.py` with MATLAB, NumPy and SciPy available.
- **Search.** On two problems small enough to enumerate (46,656 and 25,344 designs) the default search finds the global optimum. Constraints, level frequencies, duplicates, reproducibility and the running objective are checked on generated designs.
- **Orthogonal arrays.** The construction finds L8, L9, L12 (2^11), L16 (4^5), L18 (2 × 3^7) and L36 (3^4 × 2^2).
- **Ngene.** The script and design printed in Bliemer et al. (2026) give the same D-error through the translator as through a native specification. Two Bayesian designs that Ngene saved together with their scripts evaluate here to D-errors 0.0937 and 0.1433; Ngene had reported 0.0936 and 0.1433.

## Compared with designs made in Ngene

An `.ngd` file holds a design together with the script that produced it. For 15 such files from 14 of our studies (2015–2026) the script was translated with Ngene's own conventions (`'Emulate', true`: constants outside the criterion, Bayesian priors as the script says, a panel respondent answers every row, the script's weights), Ngene's design was evaluated, and a new design was generated for the same specification with the default settings, without a parallel pool, on a desktop computer. Both designs are judged on the criterion that Ngene was asked to minimise, as this toolbox computes it; for several models that is the weighted sum the script specified. `tests/bench_ngene_designs.m` repeats this for any list of files.

| Design | Rows / blocks | Models | Ngene's design | New design | Change | Seconds |
|---|---|---|---|---|---|---|
| A | 48 / 4 | Bayesian MNL | 0.1433 | 0.1384 | −3.4% | 44 |
| B | 48 / 4 | Bayesian MNL | 0.09366 | 0.08619 | −8.0% | 81 |
| C | 32 / 4 | Bayesian MNL + panel MXL | 2.527 | 2.231 | −11.7% | 244 |
| D | 20 / 5 | MNL + panel MXL | 2.666 | 2.304 | −13.5% | 58 |
| E | 18 / 1 | MNL | 0.01019 | 0.00875 | −14.2% | 1 |
| F | 36 / 3 | 2 MNL + 2 panel MXL | 6.066 | 5.907 | −2.6% | 418 |
| G | 18 / 1 | panel MXL | 5.891 | 5.549 | −5.8% | 59 |
| H | 36 / 3 | MNL + panel MXL | 2.041 | 1.849 | −9.4% | 391 |
| I | 24 / 2 | MNL | 0.6139 | 0.5561 | −9.4% | 2 |
| J | 60 / 10 | MNL | 0.06169 | 0.06119 | −0.8% | 23 |
| K | 32 / 4 | MNL | 0.2641 | 0.2607 | −1.3% | 13 |
| L | 32 / 2 | panel MXL | 5.536 | 5.777 | +4.4% | 27 |
| M | 36 / 1 | MNL (33 parameters) | 0.7544 | 0.7375 | −2.2% | 63 |
| N | 36 / 4 | MNL | 0.395655 | 0.395733 | +0.02% | 6 |
| O | 36 / 3 | MNL | 0.1321 | 0.1286 | −2.7% | 16 |

The new design has a lower criterion in 13 of 15 cases. The median reduction is 3.35% and the largest is 14.17%. The criterion rises by 0.02% for N and 4.36% for L. Generation took 1–418 seconds on the original desktop, with a median of 44 seconds. These are comparisons with saved designs. Native Ngene run times and time to a common quality threshold were not recorded. Six archived designs (A, B, E, G, I, M) violate one or more translated level-frequency bounds. Of the nine that meet the row and frequency rules, seven improved, with a median reduction of 2.62% across all nine.

An independent rerun from the [public 0.2.0 source](https://github.com/czaj/DCE-design/tree/8ea83116ca33f500f402218c0b5c22165ae4ac0f) reproduced all 15 pairs of criterion values. A fresh evaluation at seed 73 favoured the same 13 new designs. The original study files are held by the project teams.

For the two designs with only a panel model, the comparison depends on simulated respondents. Across four samples of 2,000 respondents, the new design G improved on the saved design by 0.6% to 5.8% and design L was worse by 1.4% to 4.4%.

Archived Ngene output sheets also permitted a direct check of the reported D-errors. Nine MNL values matched the toolbox to six decimal places on identical cards and priors. Four of these designs also had identical block assignments and a panel mixed logit output. On those four, the toolbox's panel D-errors (importance sampling, 512 draws per sequence, 2,000 respondents, seed 1) ranged from 5.71% below to 5.04% above the archived Ngene values. Plain prior sampling in MATLAB with the original scripts' small draw counts gave substantially different numbers from the native output.

With the toolbox's own conventions (constants in the criterion, a respondent answers one block), the two search stages for panel models can be seen separately. D-errors of the panel model, per block:

| Design | Ngene's design | After stage 1 | After stage 2 |
|---|---|---|---|
| C (4 blocks of 8) | 6.429 | 5.437 | 5.014 |
| H (3 blocks of 12) | 3.769 | 3.377 | 3.088 |
| D (5 blocks of 4) | 1.256 | 1.401 | 1.237 |

In the third case the surrogate search alone gave a design that was worse for the panel model than Ngene's; the stage on the panel criterion itself corrected that. Ngene's designs were made for the criterion of the first table, so this second table shows what the two stages do; it is not a comparison of the two programs.

## Two designs from a project-specific Python search

A separate comparison used two archived Bayesian MNL designs produced with a Python swap search, each with 48 tasks and four blocks. The original Python code reproduced both sets of cards exactly. Both tools started from the earlier survey and used the same prior distributions and task constraints. The final designs were scored in MATLAB and independently in Python using 8,192 fresh Sobol draws in each implementation.

| Alternatives | Python D-error | DCE-design D-error | Python search | DCE-design generation |
|---|---:|---:|---:|---:|
| 2 | 0.144006 | 0.143970 | 139 s | 64 s |
| 3 | 0.079201 | 0.078627 | 135 s | 130 s |

For two alternatives both searches retained exactly the same marginal level counts. D-error fell by 0.025% and the Python run took 2.18 times as long. For three alternatives Python retained the original counts in each programme column, while MATLAB could vary each count between the counts in the two original columns. The 0.73% gain includes that difference in design freedom. Times are individual desktop runs. MATLAB includes blocking and final evaluation, while the Python time covers its optimization loop.

## Conventions

- D-error = det(Ω)^(1/K): Ω is the covariance matrix for one respondent who answers every row (with blocks: one respondent per block, which is the same number of choices), K the number of parameters in the criterion. More rows give a lower D-error. D-errors are comparable only for the same model, coding and priors.
- Bayesian D-error = mean over prior draws of the D-error (as in Ngene and idefix; JMP uses the expected log-determinant instead).
- The base level of a dummy or effects coded attribute is the first level listed, unless `spec.base` names another. In translated Ngene scripts it is the last level, as in Ngene.
- Constants are part of the criterion unless `spec.constants = false` (or they are listed in `spec.exclude`). Ngene leaves them out unless `;con` is given.
- The default final Bayesian criterion integrates over the priors. A fixed-statistic evaluation uses their means. Panel refinement screens exchanges at the panel prior means and validates the selected design on the requested criterion.
- In panel models a respondent answers one block. Ngene assumes that a respondent answers every row.
- Translated Ngene scripts follow the toolbox on the last three points and say so in the notes; `dce_ngene(script, 'Emulate', true)` follows Ngene, to reproduce a value Ngene reported.

[FEATURES.md](FEATURES.md) compares the toolbox feature by feature with Ngene and other software and lists what is not implemented.

## Blocking on its own

```matlab
[y, out] = dce_block(D, 3, 'Names', names, 'Extremes', [2 6]);     % D: S-by-K level values
[~, o]   = dce_block(D, 3, 'Evaluate', existingBlockColumn);       % score an existing blocking
```

`dce_block` works on any design matrix. Its criterion is the within-block coefficient of variation of level frequencies, reported as Psi in % (100 = every block is a miniature replicate of the design). It first minimises the number of levels absent from a block, can require the lowest and highest level of chosen columns in every block, and, given priors, reports each block's share of the information and the D-error without that block. `test_dce_block.m` reproduces the published example and checks the search against exhaustive enumeration.

## Files

| Path | Content |
|---|---|
| `dce_*.m` | public functions |
| `private/` | the engine (`dd_*.m`) |
| `examples/` | runnable examples |
| `tests/` | `run_tests.m`, `bench_ngene_designs.m` (compare with designs saved by Ngene) |
| `SPEC.md`, `AGENTS.md`, `FEATURES.md` | specification reference, guide for AI agents, feature comparison |
| `CITATION.cff`, `LICENSE` | how to cite; licence (CC BY 4.0) |

Requires MATLAB R2022a or later (paged matrix functions) and the Statistics and Machine Learning Toolbox for Sobol draws (Halton draws are used without it). The Parallel Computing Toolbox is optional.

## Citing

If you use the toolbox, please cite it:

> Czajkowski, M. (2026). *DCE-design: experimental designs for discrete choice experiments in MATLAB* (version 0.2.1). https://github.com/czaj/DCE-design

and record the version, seed, specification, command, starting cards (if supplied) and draw method. The methods it builds on are listed in [FEATURES.md](FEATURES.md); the blocking criterion is that of Bliemer, Hess, Collins & Amaris (2026) and the orthogonal arrays are found with the algorithm of Xu (2002).

#
The codes are published under a Creative Commons Attribution 4.0 License. This means that you are free to use, share, or modify these codes for any purpose, even commercially. What we ask in return is that you acknowledge the source of the codes or reference one of our papers (see czaj.org/research for details).

The codes come with no warranty. We test them against published results and against other software, as described above, but errors may remain; if you find one, or have a suggestion, please let us know. Related packages: [DCE](https://github.com/czaj/DCE) (estimation of discrete choice models), [CDM](https://github.com/czaj/CDM) (count data models) and [Tools](https://github.com/czaj/Tools).
