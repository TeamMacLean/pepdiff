# Next Session: Fix GLM/ART Stratification

## Context

The `within` parameter in `compare()` is accepted but silently ignored for GLM and ART methods. The helper function `compare_within_strata()` exists and is complete, but is never called from the main dispatch.

This breaks the glm_analysis vignette, which attempts to show stratified results but currently just returns the main effect twice.

## Your Task

Fix the stratification bug following our **TDD / Ralph Loop** workflow.

## Relevant Files

- `stratification_spec.md` - Full specification with design decisions and test cases
- `R/compare.R` - Main file to modify (dispatch logic ~line 106, `compare_within_strata()` ~line 528)
- `CLAUDE.md` - Project conventions and workflow description

## Workflow: Discuss → TDD → Ralph Loop

### Step 1: Read the spec

```
Read stratification_spec.md
```

Understand the fix and test cases.

### Step 2: TDD - Write failing tests

Create `tests/testthat/test-stratification.R` with the test cases from the spec. Commit the failing tests:

```
git add tests/testthat/test-stratification.R
git commit -m "TDD: Add failing tests for GLM/ART stratification fix"
```

Verify tests fail with:
```
Rscript -e "devtools::test(filter = 'stratification')"
```

### Step 3: Clear context

Start fresh or `/clear` to maximize implementation context.

### Step 4: Ralph Loop - Implement the fix

Use the Ralph Loop to implement:

```
/ralph-loop "Fix GLM/ART stratification in pepdiff.

## Failing Tests
tests/testthat/test-stratification.R

## Spec
stratification_spec.md

## Key Changes
1. R/compare.R ~line 106: Check if within is not NULL, dispatch to compare_within_strata()
2. R/compare.R: Add bf_threshold parameter to compare_within_strata() for bayes_t support

## Verification
Rscript -e 'devtools::test(filter = \"stratification\")'
Must show: OK (all tests pass)

## Success Criteria
- All stratification tests pass
- All existing tests still pass
- devtools::check(vignettes = FALSE) has 0 errors/warnings" --completion-promise "STRATIFICATION-COMPLETE" --max-iterations 15
```

### Step 5: Verify and rebuild vignette

After implementation:

```r
devtools::test()  # All tests pass
devtools::check(vignettes = FALSE)  # 0 errors/warnings

# Rebuild the vignette to verify stratification works
rmarkdown::render('vignettes/glm_analysis.Rmd', output_dir = 'doc')
```

The vignette should now show:
- Group C: FC ≈ 1 at 0h, FC ≈ 4 at 24h (not the diluted ~2)

### Step 6: Commit

```
git add R/compare.R tests/testthat/test-stratification.R
git commit -m "Fix GLM/ART stratification to use compare_within_strata()

When 'within' parameter is specified, analysis is now run separately
for each stratum with FDR correction applied within each stratum.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

## Expected Outcome

After this fix:
- `compare(..., within = "timepoint", method = "glm")` returns separate results for each timepoint
- The glm_analysis vignette correctly shows Group C with FC ≈ 4 at 24h
- All 340+ existing tests still pass

## Notes

- The fix is straightforward - the helper function already exists and works
- Main change is just adding the dispatch check at line ~106
- Don't forget to add `bf_threshold` parameter to `compare_within_strata()` for bayes_t
