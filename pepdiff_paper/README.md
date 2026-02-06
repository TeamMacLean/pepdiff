# pepdiff Paper

Bioinformatics Application Note describing the pepdiff package.

## Target

- **Journal**: Bioinformatics (Oxford) - Application Note
- **Preprint**: bioRxiv (same content, different template)
- **Format**: ~2,000 words + 1 figure

## Companion Framing

| Package | Question | Phase |
|---------|----------|-------|
| peppwR | "How many samples do I need?" | Experimental design |
| pepdiff | "What's differentially abundant?" | Analysis |

## Files

| File | Purpose |
|------|---------|
| `outline.md` | Full paper outline with rendering pipeline |
| `references.bib` | BibTeX references (17 citations) |
| `example_application_spec.md` | Detailed spec for supplementary code |
| `paper_prompt.md` | Prompt to kick off paper writing session |
| `supplementary_example.Rmd` | Reproducible code (to be created) |
| `paper.Rmd` | Main manuscript (to be created) |

## Key Messages

1. **No imputation needed**: Gamma GLM handles missing data naturally; conventional imputation inflates FDR with MNAR data
2. **Factorial structure matters**: Stratified comparisons capture experimental design; pooled tests lose ~40% of true positives
3. **Diagnostics guide decisions**: Built-in `plot_fit_diagnostics()` identifies when to use ART instead of GLM

## Figure 1 Panels

- **A**: Workflow diagram
- **B**: Volcano plot (treatment effect at 24h)
- **C**: Method comparison (sensitivity/FDR by approach)
- **D**: Diagnostic QQ plots (good vs poor fit)

## Status

- [x] Outline complete
- [x] References collected (17 citations)
- [x] Example application spec
- [ ] supplementary_example.Rmd
- [ ] paper.Rmd
- [ ] Figure 1

## Rendering

```r
# Default (works on any system with LaTeX)
rmarkdown::render("paper.Rmd")

# bioRxiv preprint
rmarkdown::render("paper.Rmd", output_format = rticles::arxiv_article())

# Bioinformatics journal
rmarkdown::render("paper.Rmd", output_format = rticles::bioinformatics_article())
```
