# Feature: Generation Configuration and Strictness <!-- rq-9f0d6f46 -->

`GeneratorConfig` (`src/biochar/pipeline/biochar_generator.py`) is where a caller's request becomes
the numbers the pipeline builds from. It accepts composition targets directly, or a pyrolysis
temperature and feedstock to derive them from, and it decides which requests are refused outright
and which are accepted with an adjustment.

The requirements here are about **what happens to a request that cannot be honoured as written**.
Three outcomes are possible, and the difference between them is the whole subject:

- **Refused at construction.** The request is unbuildable in principle, so `GeneratorConfig(...)`
  raises before any work begins.
- **Adjusted, and the caller told.** The request is buildable only after a change this module
  makes on the caller's behalf, and a warning names the change.
- **Attempted, and the result measured.** The request is plausible, so it is attempted, and what
  was actually achieved is reported on the returned `CompositionResult`.

The failure this area produces is the third silently doing the work of the second: a request
attempted, missed by a wide margin, and reported as success. `strict=True` exists to prevent that,
and the scenarios below pin exactly how far it reaches.

## Feature API <!-- rq-a8498d8f -->

- `GeneratorConfig(target_num_carbons: int = 50, H_C_ratio: Optional[float] = None, O_C_ratio: Optional[float] = None, aromaticity_percent: Optional[float] = None, temperature: Optional[float] = None, feedstock: Optional[str] = None, functional_groups: Optional[Dict[str, int]] = None, pH: Optional[float] = None, charge_method: str = "opls", strict: bool = True, seed: Optional[int] = None, molecule_name: str = "BC", ...)` <!-- rq-e199dc57 -->
  - Resolves composition in `__post_init__`, so a constructed config carries concrete numbers and
    never `None` for `H_C_ratio`, `O_C_ratio`, or `aromaticity_percent`.
  - Raises `ValueError` for a request that cannot be built at all, in preference to failing later:
    an out-of-range `pH`, an `H_C_ratio` at or above the model's asymptotic ceiling, a
    `molecule_name` too long for the `.gro` format, an unknown `feedstock` or `charge_method`.
  - Warns, and proceeds, for a request that is buildable but questionable.

- `BiocharGenerator(config: GeneratorConfig)` <!-- rq-db4a65a9 -->
  - `generate() -> Tuple[Chem.Mol, np.ndarray, CompositionResult]` runs the five-stage pipeline.
  - The returned `CompositionResult` carries `requested_counts` alongside `placed_counts`, so the
    caller can measure the gap between what was asked for and what exists.
  - Raises `ValidationError` when `strict` is set and the structure fails validation.

- `ValidationError` <!-- rq-fda6981b -->
  - Raised in strict mode when the generated structure fails validation. Distinct from the
    `ValueError` raised at construction: a `ValueError` means the request was impossible, a
    `ValidationError` means the attempt fell short.

## Composition Resolution Has Three Levels <!-- rq-eac36cf7 -->

A composition target comes from the first of these that supplies it:

1. **An explicit value.** `H_C_ratio=0.77` is used as given.
2. **The temperature model.** With `temperature` set, any target left as `None` is derived from
   the temperature × feedstock model.
3. **A fixed default.** Anything still unset becomes `H_C_ratio` 0.5, `O_C_ratio` 0.1,
   `aromaticity_percent` 90.0.

Explicit values are not merged with derived ones, and passing `temperature` alongside an explicit
ratio is not a conflict — it derives the targets the caller did not give.

```gherkin
Feature: Resolve composition from the most specific source available

  @rq-90382d04
  Scenario: An explicitly requested ratio wins over the temperature-derived one
    Given a config with both a temperature and an explicit H/C ratio
    When the config is constructed
    Then the H/C ratio is the explicit value
    And the ratios that were not given are derived from the temperature

  @rq-c58a2ebe
  Scenario: A config with no temperature falls back to the fixed defaults
    Given a config with neither a temperature nor explicit ratios
    When the config is constructed
    Then H/C is 0.5, O/C is 0.1, and aromaticity is 90 percent
```

## The Aromatic Floor Is a Floor, Not a Prediction <!-- rq-bae6e6ca -->

The skeleton builder assembles fused aromatic sheets. Below
`MIN_BUILDABLE_AROMATICITY` = 70% it cannot faithfully realise what the temperature model
predicts, because the structure it makes is not that kind of char.

A derived aromaticity below the floor is therefore raised to it, and a warning says so. The value
the caller receives is what will be built, not what the data predicted — those differ, and the
config reports the buildable one because that is what the rest of the pipeline acts on.

This applies to *derived* aromaticity only. An explicitly requested value is the caller's
statement about what they want and is not clamped.

```gherkin
Feature: Clamp a derived aromaticity the builder cannot realise

  @rq-1762f896
  Scenario: A temperature-derived aromatic fraction below the floor is raised to the floor
    Given a temperature and feedstock whose predicted aromaticity is below 70 percent
    When the config is constructed
    Then the aromaticity is 70 percent
    And a warning names the clamp
```

## Methyl Decoration Keeps the Requested Size <!-- rq-0d7ac25a -->

When the requested H/C is above what an aromatic sheet of the requested size can carry, the
skeleton is split into a compact aromatic core plus pendant methyls, and the total carbon count is
the core plus the methyls. The split is a prediction: the methyl count follows from the core's
edge, and the core's edge is only known once it is built.

It was predicted from the wrong sheet. The split came from the ceiling of the elongated first pass,
but the core is built compact, with a lower ceiling, so it needed more methyls than were budgeted.
Wood et al.'s BC400 composition (104 carbons, H/C 0.635) came back with 117-120. The split is now
re-planned from the core actually built, and among the builds tried the one that best meets the
H/C, then the size, is kept: the H/C is what this path exists for.

The promise is one-sided. Each methyl needs an edge site, so at a high enough H/C (around 1.0 at
150 carbons) the methyls run out of sites and the structure falls short of the requested size
instead. Re-planning is not attempted then, since a smaller core only has fewer sites.

```gherkin
Feature: Keep the requested size while reaching a high H/C

  @rq-65ff7e4e
  Scenario: Methyl decoration never takes the size more than ten percent over the request
    Given a target H/C above what an aromatic sheet of the requested size can carry
    When a structure is generated across several seeds
    Then its total carbon count is no more than ten percent over the request
```

## Impossible Requests Are Refused at Construction <!-- rq-31553668 -->

Some requests cannot produce a structure at all. Those raise from `GeneratorConfig(...)`, before
a skeleton is built, so the caller learns at the point of the mistake.

`H_C_ratio` at or above **2.0** is the sharpest case. The aromatic-core-plus-methyls model
approaches H/C = 2 asymptotically and never reaches it — every structure keeps a non-zero aromatic
core — and the methyl-count formula divides by `2.0 - H_C_ratio`, which is zero at exactly 2.0.

`pH` outside its meaningful range is refused for a different reason: past the ends of the range the
Henderson–Hasselbalch fraction has saturated, so a value beyond them is far more likely a unit
error than an intent.

```gherkin
Feature: Refuse a request that cannot produce a structure

  @rq-9aa547f6
  Scenario: An unattainable H/C is refused at construction, not at generation
    Given a config with an H/C ratio of 2.0 or above
    When the config is constructed
    Then construction raises
    And the message explains that no finite structure reaches that ratio

  @rq-a9c68e17
  Scenario: A pH outside the meaningful range is refused
    Given a config with a pH below 0 or above 14
    When the config is constructed
    Then construction raises
```

## A Charge Backend That Would Erase the Charge Is Refused <!-- rq-c5bc6e18 -->

Setting `pH` asks for an ionized structure. The ML refiner constrains its predictions to sum to
zero and was trained on neutral molecules, so it would both erase the net charge the caller asked
for and extrapolate outside its training set.

That combination is refused at construction rather than honoured into a neutral structure. The
`opls` and `qm` backends both carry the formal charge through.

```gherkin
Feature: Refuse a charge backend incompatible with the requested protonation

  @rq-9bace66f
  Scenario: pH combined with the ML charge backend is refused
    Given a config with a pH and charge_method set to ml
    When the config is constructed
    Then construction raises
    And the message names an alternative backend
```

## Strictness Covers Shortfall, Not Only Failure <!-- rq-56ed9407 -->

`strict=True` means the caller wants a structure that meets the request or an exception, not a
structure that quietly does not.

A functional group is placed onto available edge sites, and a request can exceed what a skeleton
of the given size can host. Strictness treats such a shortfall as a failure of the request, not as
a successful build — a request for 40 ether bridges answered with 6 is not a structure the caller
asked for, and returning it as a success makes the shortfall invisible at exactly the point a
sweep records the run as passing.

This covers groups the caller named explicitly. When counts are derived from `O_C_ratio` instead,
the ratio is the target and the count is an implementation detail of reaching it — that case is
judged by composition validation against `O_C_tolerance`, the tolerance the caller actually set.

Requests that are attempted and merely approximate remain the normal case: composition ratios are
targets measured against a tolerance, and `CompositionResult` reports what was achieved.

```gherkin
Feature: Fail strictly on a request that was not met

  @rq-4fc714fe
  Scenario: A partially placed functional group does not satisfy strict validation
    Given a strict config explicitly requesting more functional groups than the skeleton can host
    When the structure is generated
    Then generation raises
    And the message states how many of each group were requested and placed

  @rq-ab986122
  Scenario: A ratio-derived shortfall is judged by the ratio, not the count
    Given a strict config with an O/C ratio and no explicit functional groups
    When the structure is generated
    Then a group count below the ratio-derived estimate does not raise on its own
```

## An Extrapolated Temperature Is Reported <!-- rq-32dc1242 -->

The temperature model is fitted over a finite temperature range, and that range differs by
property. A temperature outside the support of a property this config **derives** is
extrapolation, and the caller is told which property and which range rather than receiving a
silently extrapolated number.

The check is per property and scoped to what is derived here — `H_C_ratio` and `O_C_ratio`.
Properties this config never touches, such as pH and electrical conductivity, are fitted over
narrower ranges; warning about them while building a structure would attach noise to a request
that never consulted them. A caller who wants the full reference table asks the model directly,
and the model reports each property's own support.

```gherkin
Feature: Say when a derived property is being extrapolated

  @rq-a1da340c
  Scenario: A temperature beyond a derived property's support is reported
    Given a temperature outside the fitted range of a property this config derives
    When the config is constructed
    Then a warning names that property and its range

  @rq-0d326d5b
  Scenario: A temperature is not judged against a property this config never derives
    Given a temperature within the range of every derived property
    And outside the narrower range of a property this config does not derive
    When the config is constructed
    Then no extrapolation warning is emitted
```

## Questionable Combinations Warn <!-- rq-d6aaeb81 -->

Some requests are buildable and defensible but rest on parameters this project has not verified.
Ring-substituted nitrogen combined with a functional group other than hydroxyl or phenolic is one:
the resulting atom types are reachable, but the OPLS parameters for that combination are
nearest-analog choices rather than validated values.

These warn and proceed. Refusing would block legitimate exploration; silence would let an
unverified parameter reach a production run unnoticed. The warning is emitted at generation rather
than at construction, because it describes the structure about to be built rather than a request
that cannot be honoured.

```gherkin
Feature: Warn where a combination is buildable but unverified

  @rq-8155aecf
  Scenario: Crossed nitrogen doping and an unverified functional group warns
    Given a config requesting ring nitrogen alongside a functional group other than hydroxyl or phenolic
    When the structure is generated
    Then a warning names the combination
    And a structure is still produced
```

## Cross-references <!-- rq-9a108035 -->

- The functional-group placement that produces a shortfall, and the groups silently substituted
  for others, are specified in `heteroatom-assignment.md`.
- The forcefield status of the nearest-analog nitrogen parameters is specified in
  `opls-typing.md`.
- The temperature × feedstock model itself, including how a property's support range is recorded,
  belongs to `models/temperature_model.py`; see `temperature-model.md` when that document exists.
