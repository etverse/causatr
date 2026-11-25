# causatr roadmap

causatr estimates causal effects of point treatments by g-computation, inverse probability weighting (IPW), augmented IPW (AIPW), matching and structural nested models (SNM), and owns what every etverse package shares: the `contrast()` and `diagnose()` generics, the intervention constructors, the variance engines and pooling across imputations. Its time-varying code moves to [dynamatr](https://github.com/etverse/dynamatr). The [etverse design](https://etverse.github.io/design.html) fixes its scope and interface. Each item's status is in the [issues](https://github.com/etverse/causatr/issues), not in this file.

## 0.1.0

The first release comes once every etverse estimator has been checked against independent software, the bugs those checks find are fixed, and the package follows the design. Its work is the [0.1.0 milestone](https://github.com/etverse/causatr/milestone/1):

- `contrast()` and `diagnose()` become generics that the other packages extend ([#18](https://github.com/etverse/causatr/issues/18)).
- Bugs in point IPW with censoring weights ([#19](https://github.com/etverse/causatr/issues/19)), longitudinal AIPW with transport ([#20](https://github.com/etverse/causatr/issues/20)), the ICE sandwich with betareg ([#21](https://github.com/etverse/causatr/issues/21)), SNM with transport or censoring weights ([#22](https://github.com/etverse/causatr/issues/22)) and weighted balance in `diagnose()` ([#23](https://github.com/etverse/causatr/issues/23)).
- The intervention definitions and the rejections they imply ([#24](https://github.com/etverse/causatr/issues/24)).
- The numbered vignette set the design asks for, with the generated table of supported combinations ([#25](https://github.com/etverse/causatr/issues/25)).

## Outcomes and interventions

- Categorical and ordinal outcomes beyond point g-computation ([#27](https://github.com/etverse/causatr/issues/27)).
- Stochastic interventions under IPW and AIPW with a user-supplied density ([#28](https://github.com/etverse/causatr/issues/28)).
- `bernoulli_iv()` and `mixture()` constructors ([#32](https://github.com/etverse/causatr/issues/32)).

## Weights and variance

- `stabilize = "baseline"` for multivariate IPW ([#29](https://github.com/etverse/causatr/issues/29)).
- An analytic sandwich that carries the numerator model under `stabilize = "marginal"` ([#30](https://github.com/etverse/causatr/issues/30)).
- An opt-in guard that effect modifiers are baseline covariates ([#31](https://github.com/etverse/causatr/issues/31)).
- A developer interface for pooling across imputations that other packages can call ([#34](https://github.com/etverse/causatr/issues/34)).

## Diagnostics and speed

- A `diagnose()` branch for point SNM fits ([#33](https://github.com/etverse/causatr/issues/33)).
- A measured speed review of the hot paths ([#26](https://github.com/etverse/causatr/issues/26)).

## Not planned here

- Time-varying treatments: [dynamatr](https://github.com/etverse/dynamatr/issues).
- The target-trial protocol object `target_trial()`: [emulatr](https://github.com/etverse/emulatr/issues/2).
- Time-to-event estimands: survatr. Sampled designs: matchatr. Bias analysis: biasetr and negatr.
- Machine-learning nuisance models and cross-fitting.
