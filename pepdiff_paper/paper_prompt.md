# pepdiff Paper Rewrite Prompt

Use this prompt to start a new Claude Code session for rewriting the paper.

---

## Prompt

```
Rewrite the pepdiff paper (v2).

## Context

pepdiff is an R package for differential abundance analysis of phosphoproteomics data. We have:

- paper.Rmd: First draft (v1) - needs rewriting
- rewrite_plan.md: Detailed critique and restructuring plan
- supplementary_example.Rmd: Working reproducible analysis (keep as-is)
- figure_1.pdf: Generated figure (may need simplification)
- references.bib: Citations ready to use

Read rewrite_plan.md first - it contains the full critique of v1 and the new structure.

## Key Problems with v1

1. **No clear hook** - buries the lead, takes 4 paragraphs to explain what pepdiff does
2. **Code blocks in main text** - belongs in supplementary, not Application Note
3. **Reads like documentation** - lists features instead of telling a story
4. **Weak comparison** - FDR is same for conventional workflow; pooled comparison needs reframing

## The One Message

**Factorial experiments need factorial analysis. Pooled analysis detected 0% of timepoint-specific effects; stratified analysis detected 92%.**

## Key Numbers (from supplementary)

| Approach | Sensitivity | FDR |
|----------|-------------|-----|
| pepdiff (stratified GLM) | 0.92 | 0.02 |
| Conventional (t-test at 24h) | 0.84 | 0.02 |
| Complete cases only | 0.74 | 0.03 |
| Pooled (ignore timepoint) | 0.00 | — |

## New Structure (~1800 words)

### Introduction (~300 words)
- Hook: Factorial designs are standard, but analysis ignores them
- Problem: Pooled = 0% sensitivity for timepoint-specific effects
- Solution: pepdiff automates stratified GLM analysis

### Implementation (~400 words)
- Workflow in PROSE (no code blocks)
- Why Gamma GLM (one paragraph)
- Why stratification matters (one paragraph)
- Diagnostics guide method choice (one paragraph)

### Example Application (~500 words)
- Lead with key result: "0% vs 92%"
- Brief simulation description
- Comparison table
- Figure reference

### Discussion (~200 words)
- One message: respect your experimental design
- Limitations (per-peptide, cross-sectional only)
- peppwR companion

## Comparison Framing

**Critical:** Don't claim other tools CAN'T do stratified analysis. They can (with effort).

Frame as: "We compare against common proteomics workflows, not optimal expert use of existing tools. pepdiff makes appropriate analysis automatic."

Comparisons:
- Conventional workflow = what most users do (fair)
- Complete cases = conservative alternative (fair)
- Pooled pairwise = "what happens if you ignore structure" (fair, but reframe)

## Style

- British English (analyse, summarise)
- NO code blocks in main text
- Scientific but accessible tone
- Lead with results, not methods
- Cite using [@key] format

## Files to Modify

1. **paper.Rmd** - Complete rewrite following new structure
2. **figure_1.pdf** - Consider simplifying (drop workflow panel A?)

Do NOT modify supplementary_example.Rmd - it works and generates correct numbers.

## Verification

1. Render paper.Rmd to PDF
2. Check ~1800 words
3. Confirm no code blocks in main text
4. Verify key message is clear in first paragraph
```

---

## Quick Start

```bash
cd /Users/macleand/Desktop/pepdiff
# Read the rewrite plan first
cat pepdiff_paper/rewrite_plan.md
```

## Key Files

| File | Status | Action |
|------|--------|--------|
| rewrite_plan.md | Reference | Read first |
| paper.Rmd | v1 draft | Rewrite completely |
| supplementary_example.Rmd | Working | Keep as-is |
| figure_1.pdf | Generated | May simplify |
| references.bib | Ready | Use as-is |

## What Changed from v1

The supplementary and simulation are done and working. The task is purely rewriting the paper narrative:

- Remove code blocks
- Restructure with clear hook
- Lead with results (0% vs 92%)
- Reframe comparisons as "common practice"
- Tighten to ~1800 words
