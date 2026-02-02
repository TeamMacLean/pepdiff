whl # Semi-Autonomous Feature Development

A workflow for humans and Claude to collaborate on feature development with minimal back-and-forth.

## The Problem

Traditional AI-assisted development has friction:
- Human stays in the loop for every iteration
- Unclear when "done" is actually done
- Context gets consumed by discussion before implementation
- Repeated clarification questions

## The Solution: Discuss → TDD → Ralph Loop

### Phase 1: Discuss

**Goal**: Reach shared understanding of requirements before any code.

- Human describes the feature/bug
- Claude asks clarifying questions
- Agree on scope, edge cases, and what "done" looks like
- Identify relevant files and architecture patterns

**Output**: Clear requirements that can be expressed as a test.

### Phase 2: TDD (Write Failing Test First)

**Goal**: Create an objective, automated definition of "done".

1. Claude writes a test that captures the requirement
2. Run the test to confirm it **fails** (red)
3. The test IS the spec - no ambiguity about completion

**Key principle**: No implementation code until the test exists and fails.

### Phase 3: Ralph Loop (Autonomous Implementation)

**Goal**: Claude iterates autonomously until tests pass.

```
/ralph-loop "Implement [feature description].

## Failing Test
tests/e2e/my-feature.test.ts

## Relevant Files
- src/lib/services/foo.ts
- src/lib/components/Bar.svelte

## Verification
npx playwright test my-feature.test.ts
Must show: X passed

## Success Criteria
All tests pass." --completion-promise "FEATURE-COMPLETE" --max-iterations 20
```

**What happens**:
1. Claude works on implementation
2. Runs verification command
3. If tests fail, continues iterating
4. If tests pass, outputs the promise tag and loop exits
5. Human does final smoke test

---

## Context Management

### Clear Context Before Ralph Loop

Discussion consumes context. Before starting Ralph Loop:

1. **Commit the failing test** (preserves the spec in files)
2. **Start a fresh Claude session** or use `/clear`
3. **Run Ralph Loop** with full context available for implementation

This ensures Claude has maximum context for the actual work, not filled with discussion history.

### What Goes in the Ralph Loop Prompt

Since context is cleared, the prompt must be self-contained:

- **Feature description**: What we're building (brief)
- **Failing test**: Exact file path
- **Relevant files**: Where to look/modify
- **Verification command**: Exact command to run
- **Success criteria**: What "pass" looks like

The prompt references files (test, source) rather than repeating discussion content.

---

## Example Workflow

### 1. Discussion Phase

Human and Claude discuss the bug, identify root cause possibilities, agree on expected behavior.

### 2. TDD Phase

Claude writes a test:

```typescript
// tests/e2e/coordinate-input.test.ts
test('accepts NC_ accession format', async ({ page }) => {
  await page.goto('/');
  const input = page.getByPlaceholder('chr1:1000-2000');
  await input.fill('NC_000913.3:100000-100100');
  await input.press('Enter');
  // Should navigate without error
  await expect(page.locator('.error')).not.toBeVisible();
  await expect(input).toHaveValue(/NC_000913/);
});
```

Run it, confirm failure:
```
npx playwright test coordinate-input.test.ts
# 1 failed
```

### 3. Commit and Clear

```bash
git add tests/e2e/coordinate-input.test.ts
git commit -m "Add failing test for NC_ coordinate format"
```

Then `/clear` or start fresh session.

### 4. Ralph Loop Phase

```
/ralph-loop "Fix coordinate input to accept NC_ accession format.

## Failing Test
tests/e2e/coordinate-input.test.ts

## Relevant Files
- src/lib/types/genome.ts (parseCoordinate function)
- src/lib/components/Header.svelte (input handling)

## Verification
npx playwright test coordinate-input.test.ts

## Success Criteria
Test passes." --completion-promise "BUG-FIXED" --max-iterations 15
```

### 5. Human Smoke Test

Once Ralph Loop exits, human does quick manual verification in browser.

---

## Tips

### Writing Good Tests

- Test the **behavior**, not implementation details
- Include both positive cases (should work) and negative cases (should fail gracefully)
- Make assertions specific enough to catch regressions

### Scoping Ralph Loop Prompts

- List only files likely to need changes (helps Claude focus)
- Use specific test file, not full test suite (faster iterations)
- Set reasonable max-iterations (15-25 for bug fixes, 25-40 for features)

### When Things Go Wrong

If Ralph Loop hits max iterations without success:
1. Read the test output to understand what's failing
2. Check if the test itself is correct
3. May need another discussion phase to refine approach
4. Consider breaking into smaller sub-tasks

---

## Summary

| Phase | Who | Output |
|-------|-----|--------|
| Discuss | Human + Claude | Shared understanding |
| TDD | Claude | Failing test (committed) |
| Clear | Human | Fresh context |
| Ralph Loop | Claude (autonomous) | Passing implementation |
| Smoke Test | Human | Final verification |

The key insight: **tests are the contract**. They let Claude work autonomously with clear, objective success criteria.