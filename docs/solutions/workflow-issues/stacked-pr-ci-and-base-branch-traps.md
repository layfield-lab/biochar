---
title: A stacked PR reports CLEAN while running zero CI
date: 2026-07-17
category: workflow-issues
module: ci
problem_type: workflow_issue
component: development_workflow
severity: high
applies_when:
  - "Opening a PR whose base is another feature branch rather than main"
  - "Merging a PR that other PRs are stacked on"
  - "Reading a green or CLEAN status before merging"
tags:
  - ci
  - github-actions
  - stacked-prs
  - pull-requests
  - silent-failure
---

# A stacked PR reports CLEAN while running zero CI

## Context

`.github/workflows/ci.yml` filters on the **base** branch:

```yaml
on:
  pull_request:
    branches: [main]
```

A stacked PR — one whose base is another feature branch — matches nothing, so no workflow
runs. GitHub then reports:

```
mergeable:  MERGEABLE
mergeState: CLEAN
```

`CLEAN` means *"no required check is failing."* **Zero checks satisfies that trivially.** The
status looks like the one a fully-tested PR shows.

This happened twice in one session. PR #21 was stacked on #20's branch and sat at `CLEAN`
with no checks; PR #23 was stacked on #21's branch and did the same. Both would have merged
untested if the status had been taken at face value.

## Guidance

**1. On a stacked PR, `CLEAN` is not evidence.** Ask what actually ran:

```bash
gh pr checks 23
# no checks reported on the 'fix/ca-s-ca-angle' branch    <- the tell
```

**2. Retargeting to `main` does not start CI.** Changing the base fires
`pull_request: edited`, which is not in the default activity types (`opened`, `synchronize`,
`reopened`), and this workflow has no `workflow_dispatch`. Force a run with **close +
reopen** (fires `reopened`) or by pushing a commit (fires `synchronize`). Retarget *first* —
a reopen while the base is still the old branch still matches nothing.

**3. Retarget dependents BEFORE merging the base.** Merging with `--delete-branch` deletes
the base branch, and GitHub **closes** any PR based on it rather than retargeting. A closed
PR is then wedged — it can be neither retargeted nor reopened:

```
GraphQL: Cannot change the base branch of a closed pull request. (updatePullRequest)
GraphQL: Could not open the pull request. (reopenPullRequest)
```

**4. If you are already wedged, restore the base branch.** Its old tip is still reachable as
the stacked branch's parent, so nothing is lost:

```bash
# quote the refspec: zsh's :r modifier mangles "$SHA:refs/heads/..."
git push <remote-url> "<old-base-sha>:refs/heads/<deleted-base-branch>"
gh pr reopen 23
gh pr edit 23 --base main
gh api -X DELETE repos/<owner>/<repo>/git/refs/heads/<deleted-base-branch>
```

Reopen, then retarget, then delete again. The PR keeps its number, description, and thread.

## Why This Matters

The two failure modes are asymmetric, and the dangerous one is silent:

- **Merging untested code.** Nothing distinguishes "all checks passed" from "no checks ran"
  in `mergeStateStatus`. A stacked PR is exactly where this bites, and stacking is most
  attractive on *large* changes — the ones least safe to merge unverified.
- **Wedging a PR.** Loud, recoverable, but only if you know the base tip is still reachable.

## When to Apply

- Any PR whose base is not `main` in this repo.
- Any repo whose CI filters `pull_request: branches:` — the filter is on the **base**, so
  every stacked PR silently opts out.
- Before merging anything with dependents: check `gh pr list --base <branch-about-to-die>`.

## Examples

The safe sequence, dependents first:

```bash
gh pr edit 23 --base main     # retarget BEFORE the base disappears
gh pr close 23 && gh pr reopen 23   # now CI matches branches:[main] and runs
gh pr merge 21 --merge --delete-branch
# wait for #23 green, then merge it
```

The durable fix is to stop the filter from excluding stacked PRs at all — drop the base
filter so every PR is tested regardless of its base:

```yaml
on:
  push:
    branches: [main]
  pull_request:          # no branches: filter -- runs for any base
```

That trades a little CI time for removing the class. With the filter in place, every future
stacked PR re-earns this trap, and the only defence is remembering to look.

## Related

- PRs [#21](https://github.com/layfield-lab/Biochar-simulator/pull/21),
  [#22](https://github.com/layfield-lab/Biochar-simulator/pull/22),
  [#23](https://github.com/layfield-lab/Biochar-simulator/pull/23) — the stack this surfaced
  on. #23 was closed by merging #21 with `--delete-branch` and recovered as above.
- `.github/workflows/ci.yml` — the `branches: [main]` filter that causes it.
- `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md` — unrelated
  subject, same shape: a check that looks like it passed while never having run. That doc's
  guard also skips in CI, so its green badge is likewise not evidence.
