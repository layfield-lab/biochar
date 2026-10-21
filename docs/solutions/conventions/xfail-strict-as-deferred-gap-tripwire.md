---
title: Mark a deferred known gap with xfail(strict) so its fix must retire the exemption
date: 2026-07-17
category: conventions
module: tests
problem_type: convention
component: testing_framework
severity: medium
applies_when:
  - "Deferring a fix for a gap a test already detects"
  - "Writing a KNOWN_BAD / xfail exemption for something you intend to fix later"
  - "Reconciling branches that each fixed part of the same class of bug"
related_components:
  - testing_framework
tags:
  - pytest
  - xfail
  - test-coverage
  - stale-exemption
  - known-gap
---

# Mark a deferred known gap with xfail(strict) so its fix must retire the exemption

## Context

A check finds a real gap you are not going to fix in this change. The usual ways to defer it
all rot:

- A `# TODO: NC-CA-OH angle missing` comment — invisible to the suite, and stays put long
  after the gap is closed.
- A skipped test (`@pytest.mark.skip` / plain `xfail`) — goes green forever, so when someone
  fixes the underlying gap the test keeps *not running* and the now-false "known bad" note
  lingers.
- A `KNOWN_BAD = {...}` exemption dict with a prose reason — same failure: the fix lands
  elsewhere, the entry keeps passing, and the stale exemption sits there implying a bug that
  no longer exists.

The shared failure mode is **drift between the code's real state and what the test scaffolding
says about it**. It is worst across people and branches: the person who eventually fixes the
gap is often not the person who deferred it, may be on a different branch months later, and
has no reason to know the stale marker exists.

## Guidance

When you defer a gap a test can already detect, mark it **`xfail(strict=True)`** with a
reason string that names the gap and what fixing it requires.

`strict=True` is the whole point. A plain `xfail` passes whether the test fails *or*
unexpectedly passes. `strict=True` inverts the second case: once the underlying gap is
fixed, the test **XPASSes**, and a strict xpass is a **suite failure**. That failure is
impossible to ignore and can only be cleared by deleting the marker. The exemption cannot
outlive the bug it documents — the fix is forced to retire it, by whoever lands the fix,
wherever they land it.

```python
_PYRIDINIC_XFAIL = (
    "A hydroxyl on a ring carbon adjacent to a pyridinic N emits NC-CA-OH, which "
    "stock OPLS does not define. A genuine forcefield gap -- belongs in "
    "SUPPLEMENTARY_ANGLE_PARAMS with a defensible value, which needs its own change."
)

@pytest.mark.parametrize("knob", [
    pytest.param("num_pyridinic", marks=pytest.mark.xfail(strict=True, reason=_PYRIDINIC_XFAIL)),
    "num_pyrrolic",
    "num_graphitic",
])
def test_ring_nitrogen_emits_only_resolvable_terms(self, knob):
    ...
```

The reason string is not decoration. Because the fixer is forced to open this exact spot to
clear the XPASS failure, the reason is where you hand them the diagnosis and the intended
fix — it is read at exactly the moment it is useful.

Two rules keep it honest:

- **Reserve it for gaps that are genuinely deferred, not disguised bugs.** If the "gap" is
  actually a bug you could fix now, fix it. xfail is for "real, out of scope here," not "I'll
  paper over this."
- **The reason must state whether a fix is even the right move.** Two gaps can look identical
  from the test output and need opposite responses (see Examples). Record which.

## Why This Matters

The exemption and the fix are usually separated in time, author, and branch. Any mechanism
that relies on someone *remembering* to remove a stale marker will eventually leave one
behind — and a stale "known bad" note is actively misleading: the next reader trusts it and
either re-does the fix or avoids the now-safe code. `xfail(strict=True)` replaces that
memory dependency with a mechanical one: the suite goes red until the marker is gone. It is
the same principle as the rest of this project's force-field work — prefer a check that fails
loudly over a comment that has to be trusted.

## When to Apply

- Deferring a fix for something a test already catches, in this repo or any pytest project.
- Writing any "known bad" / expected-failure exemption you intend to resolve later.
- Reconciling divergent branches: mark each side's deferred gaps `xfail(strict)` so that
  merging the branch that fixes one *forces* the stale exemption out on the other side,
  rather than leaving two lineages quietly disagreeing about what is broken.

Do **not** use it for permanent, never-to-be-fixed limitations — those are documentation, not
a tripwire. `xfail(strict)` says "someone will fix this, and when they do, make them clean up
here."

## Examples

This session closed three deferred gaps, all held by `xfail(strict)`, and the mechanism
fired every time — each XPASS forced the exemption out at the moment its fix landed, across
three different authors/branches:

1. **`SS -> a carbon` in `KNOWN_BAD_MAPPINGS`.** The thioether-sulfur mapping was exempted
   with a note saying "owned by the S-mapping audit." That audit landed on another branch
   (fixing `SS -> opls_222`); the exemption started passing; the strict marker failed the
   suite until the entry was deleted. The dict is now empty — it could not accumulate stale
   entries.

2. **`carboxyl` emits `O_2-C-OH`.** Marked `xfail(strict)` as a *typing* gap. The pH-protonation
   branch fixed the carboxyl typing independently; on reconciliation the marker XPASSed and had
   to be removed. Crucially its reason string recorded **"the fix is correct typing, not a
   supplement -- supplementing would cement the wrong chemistry"** — telling the fixer what
   *not* to do.

3. **`num_pyridinic` emits `NC-CA-OH`.** Marked `xfail(strict)` as a genuine force-field gap.
   Its reason string said the opposite of carboxyl's: **"belongs in SUPPLEMENTARY_ANGLE_PARAMS
   with a defensible value."** When the value was supplied, the marker XPASSed and was removed.

The pair of #2 and #3 is the point of the second guidance rule: identical-looking test
output (`... emits terms no forcefield entry covers`), opposite correct responses — one a
typing fix, one a data supplement. The reason string is what carried that distinction to the
fixer.

## Related

- `tests/test_opls_type_map.py` — where all three markers lived; now carries none, with a
  comment explaining the tripwire so the next deferred gap follows the pattern.
- `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md` — the check that
  *found* these gaps; this convention is about how you *defer* one it finds.
- PRs [#25](https://github.com/layfield-lab/Biochar-simulator/pull/25) (retired the SS and
  carboxyl exemptions) and [#27](https://github.com/layfield-lab/Biochar-simulator/pull/27)
  (retired the pyridinic one).
