# pepdiff Paper Rewrite Plan

## Critique of Current Draft (v1)

### Major Issues

**1. No clear "hook" or message upfront**
- The introduction buries the lead - it takes 4 paragraphs before we learn what pepdiff actually does
- A reader skimming should immediately understand: "stratified GLM analysis for factorial proteomics designs without imputation"
- The key insight (pooled analysis completely fails for timepoint-specific effects) is hidden in the results

**2. Code blocks don't belong in a Bioinformatics Application Note**
- Lines 44-46, 64-80, 107-117 have code examples
- These belong in the supplementary or vignettes, not the main text
- Application Notes describe capability, not syntax

**3. Implementation section reads like documentation, not a paper**
- Lists features rather than telling a story
- "Three Analysis Methods" / "Simple and Formula Interfaces" / "Model Diagnostics" reads like a manual
- No motivation for *why* these choices matter

**4. Flow problems**
- Jumps between statistical details and software features
- The comparison with conventional workflow is weak - FDR is actually the same (0.02)!
- The "three key features" in Discussion don't match what the paper emphasised

**5. Abstract is a feature list, not a narrative**
- Too much detail on methods (Wilcoxon, bootstrap-t, Bayes factor, rank products)
- Buries the key result

### What's Missing

- **The biological motivation** - why do phosphoproteomics experiments use factorial designs?
- **A clear antagonist** - what goes wrong when you ignore factorial structure?
- **A memorable number** - "pooled analysis detects 0% of timepoint-specific effects"

---

## Competitive Landscape Analysis

### Is "Pooled Pairwise" a Straw Man?

**The concern:** We compare pepdiff against pooled pairwise Wilcoxon, which completely fails (0% sensitivity). Is this a fair comparison or a straw man that reviewers will criticise?

### What Perseus Actually Supports

Perseus [@tyanova2016perseus] is our main benchmark. It supports:

1. **Two-sample tests** (t-test, Wilcoxon) - the "conventional workflow"
2. **Multiple-sample tests** (ANOVA) - can detect overall effects
3. **Multi-group comparisons** with post-hoc tests
4. **Categorical annotation** - can stratify by factors

**Critical question:** Can Perseus do stratified analysis (treatment effect within each timepoint)?

**Answer:** Yes, but it requires manual setup:
- User must filter data to each timepoint separately
- Run separate t-tests per timepoint
- Manually combine results and apply FDR correction

This is error-prone and most users don't do it correctly.

### What MSstats Supports

MSstats [@choi2014msstats] is more sophisticated:
- Linear mixed models
- Can handle factorial designs with interaction terms
- Proper contrast specification

**But:** Requires statistical expertise to specify contrasts correctly. The barrier is knowledge, not capability.

### What limma Supports

limma [@ritchie2015limma] with design matrices:
- Full factorial model specification
- Interaction terms
- Contrast matrices

**But:** Designed for microarrays/RNA-seq; assumes log-normal after transformation. Not native to proteomics data characteristics.

### Fair Comparisons for the Paper

| Comparison | Fair? | Why |
|------------|-------|-----|
| Pooled pairwise (ignore timepoint) | **YES** | Common user error; default in many workflows |
| Per-timepoint t-test (manual stratification) | **YES** | What careful users do; pepdiff automates this |
| Perseus ANOVA | BORDERLINE | Tests overall effect, not contrasts we want |
| MSstats with correct contrasts | NO (straw man in reverse) | Requires expertise pepdiff aims to eliminate |
| limma with design matrix | NO | Different data assumptions |

### Recommended Comparison Strategy

**Keep:**
1. **Conventional workflow** (log2 + impute + t-test at 24h only) - this is what most users do
2. **Complete cases** - shows cost of excluding missing data
3. **Pooled pairwise** - BUT reframe as "what happens if you ignore factorial structure"

**Add:**
4. **Per-timepoint t-test without imputation** - the "correct manual approach" that pepdiff automates

**Reframe the narrative:**
- Don't claim Perseus CAN'T do stratified analysis
- Claim pepdiff makes it EASY and AUTOMATIC
- The comparison is against common practice, not tool capability

### Key Defensive Points for Reviewers

1. "We compare against common proteomics workflows, not optimal use of existing tools"
2. "pepdiff's contribution is making appropriate analysis accessible, not inventing new statistics"
3. "Users CAN achieve similar results with careful manual analysis; pepdiff reduces expertise required"

---

## Revised Paper Structure

### 1. Introduction (~300 words)

**Hook:** Phosphoproteomics experiments routinely use factorial designs (treatment × timepoint), but standard analysis workflows ignore this structure.

**Problem:** When treatment effects emerge at specific timepoints, pooled analysis dilutes the signal. In our simulation, pooled analysis detected 0% of timepoint-specific effects.

**Solution:** pepdiff provides stratified differential abundance analysis that respects experimental design, using Gamma GLM to handle proteomics data characteristics without imputation.

### 2. Implementation (~400 words)

**Workflow** (prose, no code):
- Import → validate → analyse → visualise
- S3 classes for consistent interface

**Why Gamma GLM:**
- Native to positive, right-skewed abundance data
- Handles missing observations without imputation
- Per-peptide fitting accommodates heterogeneous variance

**Stratified Analysis:**
- `within` parameter produces contrasts at each level of stratifying factor
- Captures timepoint-specific effects that pooled analysis misses
- Automatic—no manual subsetting required

**Method Selection:**
- GLM for most data
- ART when diagnostics indicate poor fit
- Built-in `plot_fit_diagnostics()` guides the choice

### 3. Example Application (~500 words)

**Lead with the key result:**
> "Stratified analysis detected 92% of true positives; pooled analysis detected 0%."

**Simulation design:** (brief)
- 500 peptides, 50 with true effect at 24h only
- Factorial: treatment × timepoint, 6 replicates

**Comparison table:**

| Approach | Sensitivity | Note |
|----------|-------------|------|
| pepdiff (stratified GLM) | 0.92 | Automatic stratification |
| Conventional (t-test at 24h) | 0.84 | Requires knowing where to look |
| Complete cases only | 0.74 | Loses peptides with missing data |
| Pooled (ignore timepoint) | 0.00 | Signal completely diluted |

**Figure 1:** Volcano + method comparison + diagnostics (no workflow diagram—save space)

### 4. Discussion (~200 words)

**One message:** Respect your experimental design. Factorial experiments need factorial analysis.

**Limitations:**
- Per-peptide (no information borrowing)
- Cross-sectional only (no repeated measures)

**Companion:** peppwR for power analysis completes the workflow.

---

## Action Items for Rewrite

1. [ ] Remove all code blocks from main text
2. [ ] Restructure introduction with clear hook → problem → solution
3. [ ] Rewrite Implementation as narrative, not feature list
4. [ ] Lead Example Application with the key result (0% vs 92%)
5. [ ] Simplify Figure 1 (consider dropping workflow panel)
6. [ ] Reframe comparisons as "common practice" not "tool capability"
7. [ ] Tighten abstract to one clear message
8. [ ] Add "per-timepoint t-test" comparison to show pepdiff automates best practice

---

## Session Notes

- Current supplementary_example.html and figure_1.pdf are usable
- Paper numbers match supplementary (verified)
- Need to re-render paper.Rmd after rewrite
