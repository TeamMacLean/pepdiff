# Vignette Implementation Prompt

Use this prompt to start a new session for writing pepdiff vignettes.

---

## Prompt

Implement vignettes for pepdiff according to the plan in `vignette_plan.md`.

### Context

pepdiff v2 is complete with:
- `pepdiff_data` and `pepdiff_results` S3 classes
- `read_pepdiff()` for import
- `compare()` with three methods: `"glm"` (default), `"art"`, `"pairwise"`
- Four pairwise tests: `"wilcoxon"`, `"bootstrap_t"`, `"bayes_t"`, `"rankprod"`
- Plot methods for both classes

### Vignettes to write (in order)

1. `basic_workflow.Rmd` - End-to-end introduction
2. `pairwise_tests.Rmd` - Two-group comparisons, the four tests
3. `diagnostic_plots.Rmd` - Visual QC, what good/bad looks like
4. `glm_analysis.Rmd` - Factorial designs, Gamma GLM, interaction effects
5. `art_analysis.Rmd` - Non-parametric factorial analysis

### Key guidelines

- **Synthetic data only** - generate within each vignette, no external files
- **Audience**: Intelligent but statistically uninformed. Knows basic tests and linear models. Wants practical guidance.
- **Philosophy**: Fit the model, check diagnostics, if it fits then it's appropriate. No prescriptive "use X when Y" rules.
- **No imputation**: Accept analysis failure when data is insufficient. Refer to peppwR for missingness/power.
- **Tone**: Direct, practical, "here's what to do". Acknowledge uncertainty.

### First steps

1. Run `usethis::use_vignette("basic_workflow")` to set up infrastructure
2. Implement `basic_workflow.Rmd`
3. Build and verify with `devtools::build_vignettes()`
4. Continue with remaining vignettes

### Verification

After each vignette:
```r
devtools::build_vignettes()
devtools::check(vignettes = TRUE)
```

See `vignette_plan.md` for detailed outlines of each vignette.
