# Features: this toolbox, Ngene and other design software

Version 0.2.1, 2026-11-21. Ngene is software of ChoiceMetrics; this toolbox is independent of it and of the other packages named here. Sources: the Ngene 1.4 manual and Ngene Online syntax reference; the documentation of idefix 1.1.0, spdesign 0.0.8, cbcTools 0.7.1, the SAS `%ChoicEff` / `%MktBlock` macros, Stata `dcreate`, JMP Choice Design, Sawtooth CBC, choicedesign (Python) and DCEDesignSA; and the papers listed at the end.

## Compared with Ngene 1.4

| Area | Ngene 1.4 | DCE-design 0.2 |
|---|---|---|
| **Models** | `mnl`, `rp`, `rppanel`, error components (`ec`, `ecpanel`, `rpec`, `rpecpanel`), random regret | `mnl`, `rp`, `rppanel`. No error components, no random regret. |
| **Random coefficients** | normal, uniform | normal, lognormal (either sign), uniform, triangular |
| **Panel mixed logit** | optimised on a simulated sample of `;rep` respondents who are assumed to answer every row; slow (the manual mentions hours or days for complex designs) | a respondent answers one block (or every row, `panelUnit='design'`); information over all choice sequences of a block, or over simulated respondents; the integral over the coefficients is computed for each sequence by importance sampling around its mode, with control variates and bias corrections, because one set of draws for all sequences overstates the information of long sequences. Two stages: a surrogate search (multinomial logit at the prior means), then a search on the panel criterion itself by exchanges within blocks, accepted only if it holds on a validation sample; reported values come from a third sample. |
| **Criterion** | D, A, B, S, C (WTP), imbalance, combinations | D or A as the objective. S-estimates, utility balance and choice shares are reported. No C-error. |
| **Bayesian priors** | normal, uniform; `fixed` (default), `mean`, `median` | normal, uniform; `mean` (default), `median`, `fixed` |
| **Prior draws** | Halton (default 200), Sobol, MLHS, pseudo-random, Gauss quadrature | scrambled Sobol (default) or randomised Halton. 128 draws in the search, with moves ranked on a few draws first; 2048 for the reported values. |
| **Several models** | weighted sum of raw measures | weighted sum; by default each active criterion is scaled against a single-model reference design, so weights mean relative importance |
| **Constants** | left out of the measure unless `;con` | included unless `constants = false` |
| **Bayesian priors in the objective** | only with `mean` or `median` in `;eff`; otherwise the prior means are used | integrated in the default reported criterion; selection of panel refinement moves starts at the prior means |
| **Coding** | linear, dummy, effects; base = last level | the same; base = first level unless `base` says otherwise; `dce_check` prints the mapping |
| **Interactions** | any order; level-specific terms `x.dummy[2]` | up to three varying attributes per term; no level-specific terms |
| **Search** | `swap` (default; relabel/swap/cycle), `mfederov`, `neldermead`; the user picks one; the default runs until stopped | chosen from the specification: exchanges within columns where level frequencies are fixed exactly; in addition level changes and row exchange over the feasible candidate tasks (modified Fedorov) where they have room; exhaustive for problems small enough. Vectorised batches with incremental information updates; perturb and re-optimise; several starts; stops by itself |
| **Constraints** | `;cond` only with `swap`; `;reject` / `;require` only with `mfederov` (no level balance) | `reject` and `require` with every method and together with level balance or frequency bounds; an error with a diagnosis when they cannot be met |
| **Dominance** | `*` on alternatives: identical, dominated, repeated in another order | automatic for alternatives with the same utility structure; also `maxProb` |
| **Level balance** | by construction in `swap`; count bounds `(2-4,...)`; imbalance penalty | strict by default; `frequency` with exact, minimum and maximum counts |
| **Blocking** | minimises correlation between the block number and the columns | balanced blocking (Bliemer, Hess, Collins & Amaris 2026): no level absent from a block, lowest and highest cost in every block, then level balance; reports each block's information |
| **Status quo, opt-out, scenario variables, labelled alternatives** | yes | yes (`fixed`, empty utility, `scenario`, alternative-specific parameters and `altLevels`) |
| **Evaluate an existing design** | `;eval` | `dce_evaluate`, also for `.ngd` files |
| **Start from a design** | `;start` | `Start` / `spec.start` |
| **Full factorial** | `;fact`, all rows or a random fraction, with `;reject` / `;require`; cannot be blocked | `dce_full`: all tasks or all feasible distinct tasks, or a random fraction; blocked and evaluated |
| **Orthogonal designs** | `;orth = sim`, `seq`, `seq2`, `ood`; orthogonal blocking column | `dce_orthogonal`: simultaneous and sequential, arrays of strength 2 with mixed levels found by search (Xu 2002), rows raised to the first size with an array, blocks as a column of the array where possible. No optimal orthogonal-in-the-differences designs (`ood`). |
| **Advice** | warnings on syntax and on level balance | `dce_advise`: size of the problem and time of an exhaustive search, rows and blocks against level balance, observations per parameter, priors, constraints, weakest parameter and sample size, agreement of starts, lower bound on the D-error and the design's efficiency against it |
| **Not in this toolbox** | pivot designs, covariates and segments (`;fisher`), availability and partial-profile designs through candidate files, optimal orthogonal-in-the-differences designs, continuous levels (Nelder-Mead), fold-over | — |
| **Reproducibility** | draw seeds and algorithm-specific stopping rules | a seed fixes the starts and draws for the selected draw method; store the exact call and versions, and use fixed effort for replay |
| **Automation** | GUI only | functions and a command-line entry (`dce_cli`: JSON in, JSON out) |
| **Self-test** | — | `dce_selftest`: published designs, numerical derivatives, Monte Carlo, enumeration, the panel information against simulation with very many draws |

Ngene scripts are read by `dce_ngene`. It keeps the script's content and the base-level rule, uses the toolbox's conventions for constants, Bayesian priors, the panel unit and model weights, and lists each such case in the result's notes; `'Emulate', true` keeps Ngene's conventions to reproduce Ngene's values. Of the 620 Ngene files from our studies of 2012–2026, 588 (95%) are translated, and 233 designs that were saved together with their scripts are evaluated. Of the other 32 files, 15 use features in the "not in this toolbox" row (optimal orthogonal-in-the-differences designs 7, models with different sets of alternatives 3, continuous levels 2, other syntax 3), 10 are incomplete or cannot be read, and 7 contain inconsistencies detected by the translator or evaluator: an attribute name misspelt in one utility function (5 files), a panel model without random parameters, and a saved design with levels that its script does not list.

## What other tools do, and what was taken from them

| Tool | Worth knowing | Taken |
|---|---|---|
| idefix (R) | modified Fedorov and coordinate exchange; Bayesian DB-error as the mean D-error over supplied draws; reports the share of draws with an infinite D-error; evaluate-only functions | share of draws not identified; evaluate-only mode; search on few draws and re-evaluate on many |
| spdesign (R) | Ngene-like utility strings carrying priors and levels; MNL; RSC, modified Fedorov and random search; exclusions with the Fedorov and random searches | one readable utility syntax; refusing anything that is not understood instead of ignoring it |
| cbcTools (R) | dominance by prior choice probability ("total") and by attributes ("partial"); power by simulation | `maxProb`; the Monte Carlo check of standard errors in the self-test |
| SAS `%ChoicEff`, `%MktBlock` | efficiency for a parameter subset; blocking factor chosen for D-efficiency | `exclude` with nuisance or known parameters |
| JMP, DCEDesignSA | Bayesian D-criterion as the expected log-determinant; radial-spherical quadrature; coordinate exchange; simulated annealing | not taken: the criterion here is the mean D-error, as in Ngene, so that values are comparable with past designs |
| Cuervo, Kessels, Goos & Sörensen (2016) | store each choice set's contribution and update the information matrix from the sets that change, for every prior draw | the incremental update |
| Kessels, Jones, Goos & Vandebroek (2009); Jonker & Bliemer (2019) | a small designed sample of prior points ranks designs as well as a large one | ranking moves on 8 draws, confirming on all |
| Bliemer & Rose (2010); Ngene forum | designs optimised for MNL do well under panel mixed logit, cross-sectional ones do not; optimise for MNL and evaluate under the panel model | the default panel surrogate |
| Bliemer, Rose & Chorus (2017) | dominance from prior signs; row-based search handles task constraints | the dominance rule (tested on their four designs) |

## Possible extensions

In rough order of usefulness for the designs made here:

1. Panel stage: exchanges between blocks as well (now only within blocks, to keep the blocks' level balance), Bayesian priors in that stage (now the prior means), and a second-order control variate for uniform and triangular coefficients (now first order; normal and lognormal coefficients have the exact one).
2. Speed for models with many parameters: exact evaluation of a move through the matrix determinant lemma (a move changes the information by a matrix of low rank), and pruning by the first-order bound that the concavity of the log-determinant gives. GPU evaluation of the batches.
3. C-error (WTP) and S-estimate as objectives.
4. Dominance checks that involve a status-quo alternative; partial profiles and overlap constraints.
5. Optimal orthogonal-in-the-differences designs; orthogonal arrays from a catalogue for sizes where the search does not find one.
6. Error components; pivot designs; covariates and segments.

## References

- Bliemer, Hess, Collins & Amaris (2026). On the blocking of non-orthogonal experimental designs for choice experiments. ITLS-WP-26-22.
- Bliemer & Rose (2010). Construction of experimental designs for mixed logit models allowing for correlation across choice observations. Transportation Research Part B 44, 720–734.
- Bliemer, Rose & Chorus (2017). Detecting dominance in stated choice data and accounting for dominance-based scale differences in logit models. Transportation Research Part B 102, 83–104.
- Bliemer, Rose & Hess (2008). Approximation of Bayesian efficiency in experimental choice designs. Journal of Choice Modelling 1, 98–126.
- ChoiceMetrics (2024). Ngene 1.4 User Manual & Reference Guide.
- Cook & Nachtsheim (1980). A comparison of algorithms for constructing exact D-optimal designs. Technometrics 22, 315–324.
- Cuervo, Kessels, Goos & Sörensen (2016). An integrated algorithm for the optimal design of stated choice experiments with partial profiles. Transportation Research Part B 93, 648–669.
- Jonker & Bliemer (2019). On the optimization of Bayesian D-efficient discrete choice experiment designs for the estimation of QALY tariffs that are corrected for nonlinear time preferences. Value in Health 22, 1162–1169.
- Kessels, Jones, Goos & Vandebroek (2009). An efficient algorithm for constructing Bayesian optimal choice designs. Journal of Business & Economic Statistics 27, 279–291.
- Kuhfeld (2010). Marketing Research Methods in SAS (MR-2010).
- Traets, Sanchez & Vandebroek (2020). Generating optimal designs for discrete choice experiments in R: the idefix package. Journal of Statistical Software 96(3).
- Xu (2002). An algorithm for constructing orthogonal and nearly-orthogonal arrays with mixed levels and small runs. Technometrics 44, 356–368.
- Zwerina, Huber & Kuhfeld (1996, in MR-2010E). A general method for constructing efficient choice designs.
