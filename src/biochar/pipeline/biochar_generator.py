"""
Biochar Structure Generator

Main API for generating biochar building blocks for GROMACS simulations.
"""

import dataclasses
import logging
import random
from typing import Callable, Optional, Tuple, List, Dict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from rdkit import Chem
import numpy as np

logger = logging.getLogger(__name__)

from .carbon_skeleton import PAHAssembler, SkeletonValidator, _edge_carbon_fraction
from .heteroatom_assignment import (
    OxygenAssigner,
    NitrogenSubstitutor,
    HydrogenAssigner,
    HeteroatomValidator,
    CompositionResult,
    attach_aliphatic_carbons,
    _fix_heteroatom_bond_types,
)
from .multi_unit import (
    LINKER_TYPES,
    MIN_UNIT_CARBONS,
    LinkedSkeleton,
    build_linked_skeleton,
    embed_linked_units,
    plan_units,
    unit_bridging_atoms,
)


class ValidationError(Exception):
    """Raised in strict mode when the generated structure fails validation."""
from .geometry_3d import CoordinateGenerator, GeometryValidator, ClashResolver
from .opls_typing import AtomTyper, ChargeAssigner
from ..export.gromacs_export import GromacsExporter
from .validation import ValidationEngine
from ..constants import MIN_BUILDABLE_AROMATICITY

# Rebuilds of the compact core while the aromatic/methyl split settles (see
# _generate_carbon_skeleton). Each re-plans from the ceiling of the core just
# built; two or three usually agree, and the closest-to-target build is kept.
_CORE_SPLIT_ITERATIONS = 4


@dataclass
class GeneratorConfig:
    """
    Configuration for :class:`BiocharGenerator`.

    All parameters are optional — defaults produce a 50-carbon, mostly-aromatic
    biochar with a moderate hydrogen content and no oxygen.

    Attributes:
        target_num_carbons: Target number of carbon atoms in the skeleton.
            The generator grows the PAH graph until it reaches a count within
            *size_tolerance* of this value.  Minimum 6 (benzene).
        size_tolerance: Fractional tolerance on *target_num_carbons*
            (e.g. 0.10 = ±10 %).
        H_C_ratio: Target hydrogen-to-carbon atomic ratio.
        H_C_tolerance: Fractional tolerance on *H_C_ratio*.
        O_C_ratio: Target oxygen-to-carbon atomic ratio.  Ignored when
            *functional_groups* is not ``None``.
        O_C_tolerance: Fractional tolerance on *O_C_ratio*.
        aromaticity_percent: Target fraction of carbon atoms that are
            aromatic, as a percentage (0–100).
        aromaticity_tolerance: Absolute tolerance on *aromaticity_percent*
            in percentage-point units.
        functional_groups: Explicit dict mapping functional group name to
            exact placement count, e.g. ``{"phenolic": 3, "carboxyl": 1}``.
            Supported groups: ``phenolic``, ``hydroxyl``, ``carboxyl``,
            ``ether``.  If ``None``, total oxygen is derived from *O_C_ratio*
            and placed as phenolic groups.
        periodic_box: If ``True``, include periodic boundary box vectors in
            the exported ``.gro`` file.
        box_size: Explicit box size in nm as a 3-element array.  Used only
            when *periodic_box* is ``True`` and box_size is not ``None``.
        molecule_name: Residue name written to ``.gro`` / ``.itp`` (max 5
            characters — GROMACS hard limit).  Suggested naming: ``BC400``,
            ``BC600``, ``BC800`` (pyrolysis temperature series).
        seed: Integer random seed for reproducibility.  ``None`` = random.
        defect_fraction: Probability [0, 1) that each ring added during
            skeleton growth is a 5-membered pentagon rather than a hexagon.
            0.0 = pure hexagonal PAH.  Values 0.05–0.20 introduce realistic
            topological disorder seen in low-temperature biochar.
        max_ether_span: Maximum C–C shortest-path distance (in bonds) between
            the two ring carbons bridged by each ether oxygen.  Controls the
            ring size of the C–O–C bridge (ring size = max_ether_span + 2).
            Default 3 → 5-membered furan-like ring (always geometrically flat).
            Larger values may fold the aromatic sheet.
    """

    # Size parameters
    target_num_carbons: int = 50
    size_tolerance: float = 0.10

    # Composition parameters.
    # H_C_ratio / O_C_ratio / aromaticity_percent default to None and are resolved
    # in __post_init__:  explicit value  >  (temperature, feedstock)-derived  >
    # the historical hard defaults (0.5 / 0.1 / 90.0).  This lets `temperature`
    # drive composition for every entry point without breaking callers that omit it.
    H_C_ratio: Optional[float] = None
    H_C_tolerance: float = 0.10
    O_C_ratio: Optional[float] = None
    O_C_tolerance: float = 0.10

    # Structural parameters
    aromaticity_percent: Optional[float] = None
    aromaticity_tolerance: float = 5.0

    # Functional groups — dict mapping group name → exact count to place.
    # e.g. {"phenolic": 3, "carboxyl": 1}
    # Valid keys: phenolic, hydroxyl, carboxyl, ether,
    #             carbonyl*, quinone*, lactone*   (* fall back to phenolic)
    # If None, O_C_ratio controls total oxygen using phenolic groups.
    functional_groups: Optional[Dict[str, int]] = None

    # System setup
    periodic_box: bool = False
    box_size: Optional[np.ndarray] = None

    # Naming and identification
    # Residue name (max 5 chars for GROMACS .gro format)
    # Suggested: BC400, BC600, BC800 (temperature), BCH05, BCO10 (composition)
    molecule_name: str = "BC"

    # Random seed for reproducibility
    seed: Optional[int] = None

    # Environmental pH.  None (default) = no protonation stage at all, i.e.
    # every group stays in the neutral form it is built in, reproducing the
    # pre-pH behaviour exactly.
    #
    # When set, each titratable site (carboxyl, phenolic, thiol, aniline N,
    # pyridinic N) is independently ionized with its Henderson-Hasselbalch
    # probability, drawn from `seed`.  The molecule then carries a real net
    # charge; neutralising the system is left to `genion -neutral` at solvation
    # time (see md_setup), not to the molecule definition.
    #
    # A single structure is one sample from the ensemble.  Near a pKa each site
    # is close to a coin flip, so a pH within ~1 unit of a pKa needs replicate
    # seeds to be representative -- `composition.ionized_counts` reports what
    # was actually placed.
    pH: Optional[float] = None

    # Ring defects: probability [0, 1) that each ring addition is a pentagon.
    # 0.0 = pure hexagonal PAH (default).  0.1 ≈ 10% pentagons.
    defect_fraction: float = 0.0

    # Ring defects: probability [0, 1) that each ring addition is a 7-membered
    # heptagon, introducing the negative (saddle) curvature of non-graphitizing
    # carbons (Wood et al. 2024).  0.0 = no heptagons (default).  Wood's ratios
    # are ~1 pentagon : 5 hexagons and ~1 heptagon : 10 hexagons, i.e.
    # defect_fraction ≈ 0.154 with heptagon_fraction ≈ 0.077.
    heptagon_fraction: float = 0.0

    # Maximum graph distance (in bonds along the carbon skeleton) between the
    # two ring carbons that an ether oxygen may bridge.  The ring formed by the
    # C-O-C bridge has (max_ether_span + 2) members.  Minimum enforced = 3
    # (5-membered, furan-like) to avoid strained smaller rings.
    #   3 → 5-membered ring (furan/benzofuran-like) ← default, always flat
    #   4 → 6-membered ring (pyran/chromene-like)   — may cause minor strain
    #   5 → 7-membered ring                         — risk of sheet folding
    # Values > 5 risk cross-sheet bridges that fold the sheet into a nanotube.
    max_ether_span: int = 3

    # Ring-substituting nitrogen doping (counts of ring C atoms replaced by N).
    #   num_pyridinic — edge 6-ring N, no H (pyridine-like)
    #   num_pyrrolic  — 5-ring N-H (requires defect_fraction > 0 for pentagons)
    #   num_graphitic — interior 6-ring N, no H (quaternary / graphitic)
    # Substitution runs after oxygen assignment, before hydrogen saturation.
    # If too few suitable sites exist, as many as possible are placed and a
    # warning is logged (does not raise).
    num_pyridinic: int = 0
    num_pyrrolic: int = 0
    num_graphitic: int = 0

    # Partial charge assignment method.
    # "opls"  — static OPLS-AA lookup table (default, no extra dependencies).
    # "ml"    — environment-aware Gaussian Process model; requires the ``ml``
    #           optional extra (scikit-learn).  Falls back to a model trained on
    #           OPLS reference charges when the bundled artifact is absent.
    # "qm"    — LigParGen-style QM charges: AM1 (via an external MOPAC binary) →
    #           CM1A mapping → ×1.14 scaling.  Requires ``mopac`` on PATH
    #           (``conda install -c conda-forge mopac``).  See
    #           :mod:`biochar.charges.qm_charges`.
    charge_method: str = "opls"

    # Data-driven composition by pyrolysis temperature (°C) and optional feedstock.
    # When `temperature` is set, any of H_C_ratio / O_C_ratio / aromaticity_percent
    # left as None are filled from the UC Davis Biochar Database model
    # (:mod:`biochar.models.temperature_model`).  `feedstock` (one of
    # ``temperature_model.VALID_FEEDSTOCKS``: softwood, hardwood, grass, manure,
    # corn_stover, wood) selects a feedstock-specific H/C·O/C curve where the data
    # supports it, otherwise the pooled curve is used.  A data-derived aromaticity
    # below MIN_BUILDABLE_AROMATICITY is clamped to that floor with a warning.
    temperature: Optional[float] = None
    feedstock: Optional[str] = None

    # Strict validation mode (default True).
    # When True, :meth:`BiocharGenerator.generate` raises :class:`ValidationError`
    # if:
    #   * any explicitly-requested functional group type is completely unplaceable
    #     (zero instances placed despite >0 requested), or
    #   * the final structure fails composition or geometry validation (ratio
    #     drift beyond tolerance, or steric clashes remaining after resolution).
    # Set to False to restore the previous lenient behaviour (print-and-continue).
    strict: bool = True

    # Aliphatic (sp3) carbon decoration.  A purely aromatic flake caps hydrogen
    # only on its perimeter, so its H/C is bounded (~0.5 even fully elongated).
    # When True (default), if the requested H_C_ratio exceeds that aromatic
    # ceiling the generator builds a smaller aromatic core and attaches pendant
    # methyl groups so the total carbon count still matches target_num_carbons
    # while the H/C reaches the target -- matching the aliphatic content of
    # low-temperature biochar.  Set False (or request aromaticity_percent ≥ 99)
    # to force a pure-aromatic structure (H/C then capped, with a warning).
    allow_aliphatic: bool = True

    # Aliphatic oxygen placement.  In O_C_ratio mode, when the aromatic edge
    # sites cannot hold the requested oxygen (typical of low-aromaticity chars
    # whose edges are consumed by aliphatic decoration + H-saturation), the
    # remainder is placed as aliphatic hydroxyls (-CH2-OH) on the sp3 carbons.
    # This is what lets high-O/C low-temperature points reach their O/C target.
    # Set False to keep oxygen on aromatic edges only (the pre-existing
    # behaviour); no effect when there are no aliphatic carbons.
    allow_aliphatic_oxygen: bool = True

    # Linked multi-unit ("island") structures.  Off by default: every request
    # that does not set this builds exactly the one fused sheet it always did.
    # When True, the skeleton is several separate aromatic units joined in a tree
    # by biaryl bonds, diaryl / aryl-alkyl ethers and propyl bridges, the
    # architecture of Wood et al. 2024's low-temperature blocks.  The unit count
    # follows `temperature` (below 500 °C 4-5 units, 500-700 °C 3, 700 °C and
    # above 1 -- the fused sheet) unless `num_units` sets it; `linker_mix` maps
    # linker names (biaryl, diaryl_ether, aryl_alkyl_ether, propyl) to relative
    # weights and otherwise follows the same bands.  See
    # :mod:`biochar.pipeline.multi_unit` and rqm/linked-units.md.
    linked_units: bool = False
    num_units: Optional[int] = None
    linker_mix: Optional[Dict[str, float]] = None

    def __post_init__(self):
        # functional_groups defaults to None → O_C_ratio-driven phenolic placement
        # (no default list here; OxygenAssigner handles the None case)

        # Validate molecule name length
        if len(self.molecule_name) > 5:
            raise ValueError(f"molecule_name must be ≤5 characters (GROMACS .gro format), got '{self.molecule_name}' ({len(self.molecule_name)} chars)")

        if self.charge_method not in ("opls", "ml", "qm"):
            raise ValueError(
                f"charge_method must be 'opls', 'ml', or 'qm', got '{self.charge_method}'"
            )

        self._validate_linked_units()

        if self.max_ether_span < 3:
            raise ValueError(
                f"max_ether_span must be ≥ 3 (minimum for a 5-membered ring), "
                f"got {self.max_ether_span}"
            )

        if self.pH is not None:
            from ..constants import PH_MAX, PH_MIN

            if not PH_MIN <= self.pH <= PH_MAX:
                raise ValueError(
                    f"pH must be within [{PH_MIN}, {PH_MAX}], got {self.pH}"
                )

            # The ML refiner forces its predictions to sum to zero and was
            # trained on neutral molecules, so it would both erase the net
            # charge and extrapolate outside its training set.  Failing here is
            # better than silently handing back a neutral structure the caller
            # explicitly asked to be ionized.  'opls' and 'qm' both honour the
            # formal charge ('qm' passes it to MOPAC as CHARGE=).
            if self.charge_method == "ml":
                raise ValueError(
                    "charge_method='ml' cannot be combined with pH: the ML "
                    "refiner constrains total charge to zero and is trained on "
                    "neutral molecules, so it would erase the very charge pH "
                    "creates. Use charge_method='opls' (default) or 'qm'."
                )

        # --- resolve composition: explicit > (temperature,feedstock)-derived > default ---
        if self.feedstock is not None:
            from ..models.temperature_model import VALID_FEEDSTOCKS
            if self.feedstock not in VALID_FEEDSTOCKS:
                raise ValueError(
                    f"feedstock must be one of {VALID_FEEDSTOCKS} or None, "
                    f"got {self.feedstock!r}"
                )
        if self.temperature is not None:
            from ..models.temperature_model import get_default_model
            model = get_default_model()
            # Checked per property, and only for the properties actually derived
            # here. Properties are fitted over different ranges, so one answer
            # cannot describe them all; and warning about pH or conductivity --
            # which this config never derives -- would be noise attached to a
            # structure request.
            for prop in ("H_C_ratio", "O_C_ratio"):
                valid_range = model.get_valid_range(prop, self.feedstock)
                if valid_range is None:
                    continue
                t_min, t_max = valid_range
                if not (t_min <= self.temperature <= t_max):
                    logger.warning(
                        "temperature=%.0f°C is outside the data range [%.0f–%.0f°C] "
                        "for %s (feedstock=%r). Predictions are extrapolated and "
                        "may be unreliable.",
                        self.temperature, t_min, t_max, prop, self.feedstock,
                    )
            comp = model.composition(self.temperature, self.feedstock)
            if self.H_C_ratio is None:
                self.H_C_ratio = comp["H_C_ratio"]
            if self.O_C_ratio is None:
                self.O_C_ratio = comp["O_C_ratio"]
            if self.aromaticity_percent is None:
                arom = comp["aromaticity_percent"]
                if arom < MIN_BUILDABLE_AROMATICITY:
                    logger.warning(
                        "Predicted aromaticity %.0f%% (T=%s°C, feedstock=%s) is below "
                        "the PAH-buildable floor (%.0f%%); clamping. The aromatic-sheet "
                        "model poorly represents such a low-aromaticity char.",
                        arom, self.temperature, self.feedstock, MIN_BUILDABLE_AROMATICITY,
                    )
                    arom = MIN_BUILDABLE_AROMATICITY
                self.aromaticity_percent = arom
        # Fill any still-unset composition fields with the historical defaults.
        if self.H_C_ratio is None:
            self.H_C_ratio = 0.5
        if self.O_C_ratio is None:
            self.O_C_ratio = 0.1
        if self.aromaticity_percent is None:
            self.aromaticity_percent = 90.0

        # H/C = 2.0 is the asymptotic ceiling of the aromatic-core-plus-methyls
        # model used by aliphatic decoration (_generate_carbon_skeleton): H/C
        # approaches but never reaches 2 as the methyl count grows, since every
        # structure retains a nonzero aromatic core. At exactly 2.0 the methyl-
        # count formula divides by (2.0 - H_C_ratio) == 0; above it the request
        # is unbuildable by this generator regardless.
        if self.H_C_ratio >= 2.0:
            raise ValueError(
                f"H_C_ratio must be < 2.0 (the asymptotic ceiling for this "
                f"generator's aromatic-core-plus-methyls model; no finite "
                f"biochar structure reaches it), got {self.H_C_ratio}"
            )

        if self.max_ether_span > 4:
            logger.warning(
                "max_ether_span=%d may fold the aromatic sheet into a nanotube-like "
                "structure. Values > 4 risk cross-sheet bridges. "
                "Recommended: 3 (5-membered furan-like ring).",
                self.max_ether_span,
            )

    def _validate_linked_units(self) -> None:
        if not self.linked_units:
            # Ignoring these would hand back a fused sheet to a caller who
            # described a linked one.
            for name in ("num_units", "linker_mix"):
                if getattr(self, name) is not None:
                    raise ValueError(
                        f"{name} applies only to linked multi-unit structures; "
                        f"set linked_units=True to use it"
                    )
            return
        if self.num_units is not None:
            if isinstance(self.num_units, bool) or not isinstance(self.num_units, int) \
                    or self.num_units < 1:
                raise ValueError(f"num_units must be an integer >= 1, got {self.num_units!r}")
            if self.target_num_carbons < MIN_UNIT_CARBONS * self.num_units:
                raise ValueError(
                    f"num_units={self.num_units} needs at least "
                    f"{MIN_UNIT_CARBONS * self.num_units} carbons (a benzene per unit); "
                    f"target_num_carbons is {self.target_num_carbons}"
                )
        if self.linker_mix is not None:
            unknown = sorted(set(self.linker_mix) - set(LINKER_TYPES))
            if unknown:
                raise ValueError(
                    f"unknown linker type(s) {unknown}; expected some of {list(LINKER_TYPES)}"
                )
            weights = list(self.linker_mix.values())
            if any(w < 0 for w in weights) or sum(weights) <= 0:
                raise ValueError(
                    f"linker_mix weights must be non-negative with a positive sum, "
                    f"got {self.linker_mix}"
                )

    def to_dict(self) -> dict:
        """Return a JSON-serializable dictionary of this configuration."""
        d = dataclasses.asdict(self)
        if d.get("box_size") is not None:
            d["box_size"] = list(d["box_size"])
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "GeneratorConfig":
        """Reconstruct a GeneratorConfig from a plain dictionary (e.g. loaded from JSON)."""
        d = dict(d)
        if d.get("box_size") is not None:
            d["box_size"] = np.array(d["box_size"])
        return cls(**d)


@dataclass
class BiocharResult:
    """
    Named result returned by :func:`generate_biochar`.

    Fields are accessible by name (``result.mol``, ``result.gro_path``, …).
    The object also supports positional unpacking for backward compatibility::

        mol, coords, gro, top, itp = generate_biochar(...)  # still works

    When ``write_files=False``, the three path fields are ``None``.
    """

    mol: Chem.Mol
    coords: np.ndarray
    composition: CompositionResult
    gro_path: Optional[Path]
    top_path: Optional[Path]
    itp_path: Optional[Path]
    ring_composition: Optional[Dict[str, int]] = None
    """Ring-type breakdown of the carbon skeleton, e.g. ``{"hexagons": 14, "pentagons": 2}``."""
    linked_units: Optional["LinkedSkeleton"] = None
    """Unit sizes and linker census when ``linked_units=True`` built several units; else None."""

    def __iter__(self):
        return iter((self.mol, self.coords, self.gro_path, self.top_path, self.itp_path))


class BiocharGenerator:
    """
    Generate a single biochar molecule and export it to GROMACS files.

    The generator runs a five-step pipeline:

    1. **Carbon skeleton** — grows a PAH graph to the requested carbon count
       using hexagonal ring expansion (or defective with pentagons).
    2. **Heteroatom assignment** — places oxygen-containing functional groups
       then fills remaining valences with hydrogen.
    3. **3D coordinates** — embeds the molecule in 3D; flattens large sheets
       via the hex-lattice path and optimises O–H hydrogen positions.
    4. **OPLS-AA typing** — assigns atom types and partial charges.
    5. **Validation** — checks composition ratios and geometry.

    Use :func:`generate_biochar` for a one-call convenience wrapper.

    Examples::

        config = GeneratorConfig(target_num_carbons=80, H_C_ratio=0.4,
                                 O_C_ratio=0.1, seed=42)
        gen = BiocharGenerator(config)
        mol, coords, composition = gen.generate()
        gro, top, itp = gen.export_gromacs(output_directory="output")
    """

    def __init__(self, config: Optional[GeneratorConfig] = None):
        """
        Initialise the generator.

        Args:
            config: Generator configuration.  Uses :class:`GeneratorConfig`
                defaults if ``None``.
        """
        self.config = config or GeneratorConfig()
        self.mol = None
        self.coords = None
        self.composition = None
        self.atom_types = None
        self.charges = None
        self.validation_report = None
        self.ring_composition: Optional[Dict[str, int]] = None
        # Census of a linked multi-unit build; None for the fused sheet.
        self.linked_units: Optional[LinkedSkeleton] = None

    # Functional groups whose adjacency to a ring-substituted N is covered by the
    # depth-3 bonded-resolution check (tests/test_opls_type_map.py). That check
    # exercises each N-doping mode at the default composition, which places
    # hydroxyl / phenolic O -- so those crossings are verified (NC-CA-OH is
    # supplied in SUPPLEMENTARY_ANGLE_PARAMS). Any other functional group crossed
    # with N-doping is NOT systematically checked; see the warning below.
    _DOPING_VERIFIED_GROUPS = frozenset({"hydroxyl", "phenolic"})

    def _warn_unverified_crossed_doping(self) -> None:
        """Warn when the config crosses ring-N doping with an unverified group.

        Ring-substituted N (pyridinic / pyrrolic / graphitic) combined with an
        explicitly requested functional group other than hydroxyl / phenolic can
        emit a bond or angle that neither stock oplsaa.ff nor
        SUPPLEMENTARY_ANGLE_PARAMS covers -- e.g. an ether O next to a pyridinic
        N. Nothing verifies these crossed combinations today (the depth-3 check
        varies one axis at a time), so grompp may reject the topology at
        simulation time with "No default Angle types" / "No default Bond types".

        This is a heads-up, not an error: the structure is still generated. If
        grompp does reject it, supply the missing parameter in
        SUPPLEMENTARY_ANGLE_PARAMS with a provenanced value, per
        docs/solutions/conventions/verify-opls-types-against-real-forcefield.md.
        Silence it with logging (this logger is 'biochar.pipeline.biochar_generator').
        """
        cfg = self.config
        doped = (cfg.num_pyridinic or 0) + (cfg.num_pyrrolic or 0) + (cfg.num_graphitic or 0)
        if doped <= 0 or not cfg.functional_groups:
            return
        crossed = sorted(
            g for g in cfg.functional_groups
            if g not in self._DOPING_VERIFIED_GROUPS
        )
        if not crossed:
            return
        logger.warning(
            "Crossed N-doping + functional group requested (%s alongside ring N). "
            "These combinations are not verified against the forcefield -- the "
            "topology may emit a bond/angle stock oplsaa.ff lacks, which grompp "
            "rejects at simulation time. If that happens, add the missing term to "
            "SUPPLEMENTARY_ANGLE_PARAMS (see the OPLS verification convention doc).",
            ", ".join(crossed),
        )

    def generate(self) -> Tuple[Chem.Mol, np.ndarray, CompositionResult]:
        """
        Run the full generation pipeline and return the molecular structure.

        Returns:
            Tuple of:

            * **mol** (:class:`rdkit.Chem.Mol`) — molecule with 3-D
              conformer and OPLS-AA atom types assigned.
            * **coords** (:class:`numpy.ndarray`, shape ``(N, 3)``) —
              atomic coordinates in Ångströms.
            * **composition** (:class:`~heteroatom_assignment.CompositionInfo`)
              — atom counts, H/C and O/C ratios, and functional-group census.

        Raises:
            RuntimeError: If carbon skeleton growth fails after retries.
        """
        self._warn_unverified_crossed_doping()

        # Step 1: Generate carbon skeleton
        logger.info("Generating carbon skeleton with %d carbons...", self.config.target_num_carbons)
        skeleton = self._generate_carbon_skeleton()
        self.ring_composition = skeleton.ring_composition

        # Step 2: Assign heteroatoms (O, then ring N substitution, then H)
        logger.info("Assigning heteroatoms...")
        mol, comp_result = self._assign_oxygens(skeleton.mol)
        mol = self._substitute_nitrogens(mol, comp_result)
        # Protonation sits here by necessity: after every heteroatom exists (so
        # there are sites to titrate) and before HydrogenAssigner, which owns
        # acidic-H placement and would otherwise fight the decision.
        mol = self._assign_protonation(mol, comp_result)
        mol = self._assign_hydrogens(mol, comp_result)

        # Step 3: Generate 3D coordinates
        logger.info("Generating 3D coordinates...")
        mol, coords = self._generate_geometry(mol)

        # Geometry re-perceives aromaticity on this molecule with RDKit's
        # default model, which marks ether C-O bonds AROMATIC when the bridge
        # closes a furan-like ring.  Restore correct single-bond types before
        # typing.
        mol = _fix_heteroatom_bond_types(mol)

        # Step 4: Assign OPLS types and charges
        logger.info("Assigning OPLS-AA atom types and charges...")
        self._assign_opls_properties(mol, coords)

        # Step 5: Validate
        logger.info("Validating structure...")
        self._validate(mol, comp_result, coords)

        # Validation reads the molecule and no longer rewrites it -- its
        # sanitisation check runs on a copy -- so this is a guard rather than a
        # repair, kept because the cost is a bond-type pass and the failure it
        # would catch is silent: an ether C-O marked AROMATIC types as OS or OH
        # by chance and the wrong parameter reaches the topology.
        mol = _fix_heteroatom_bond_types(mol)

        # Store results
        self.mol = mol
        self.coords = coords
        self.composition = comp_result

        return mol, coords, comp_result

    def export_gromacs(
        self,
        output_directory: str = ".",
        basename: str = "biochar",
    ) -> Tuple[Path, Path, Path]:
        """
        Write GROMACS structure and topology files.

        Must be called after :meth:`generate`.

        Args:
            output_directory: Directory in which to write output files.
                Created if it does not exist.
            basename: Stem for output filenames
                (e.g. ``"bc400"`` → ``bc400.gro``, ``bc400.top``,
                ``bc400.itp``).

        Returns:
            Tuple of :class:`~pathlib.Path` objects
            ``(gro_path, top_path, itp_path)``.

        Raises:
            RuntimeError: If :meth:`generate` has not been called yet.
        """
        if self.mol is None or self.coords is None:
            raise RuntimeError("Must call generate() before export_gromacs()")

        exporter = GromacsExporter(output_directory)
        gro_path, top_path, itp_path = exporter.export(
            self.mol,
            self.coords,
            self.atom_types,
            self.charges,
            molecule_name=self.config.molecule_name,
            basename=basename,
            include_periodic_box=self.config.periodic_box,
            box_size=self.config.box_size,
            title=self._structure_title(),
            provenance=self._provenance_lines(),
        )

        logger.info("GROMACS files written: %s | %s | %s", gro_path, top_path, itp_path)

        return gro_path, top_path, itp_path

    def _structure_title(self) -> str:
        """The `.gro` title line: what distinguishes this structure.

        Two structures built at pH 3 and pH 11 have the same skeleton, the same
        functional groups, and differ only in their protonation -- which the
        title is the only human-readable record of. A timestamp distinguishes
        them from nothing.

        The pH appears only when one was requested. A structure that was never
        titrated has no pH to report, and stamping a default would be a claim
        about chemistry nobody made.
        """
        stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        parts = [f"Biochar structure generated on {stamp}"]
        if self.config.pH is not None:
            charge = getattr(self.composition, "net_charge", None)
            charge_note = f", net charge {charge:+d}" if charge is not None else ""
            parts.append(f"titrated at pH {self.config.pH:g}{charge_note}")
        return ", ".join(parts)

    def _provenance_lines(self) -> List[str]:
        """Header comments for the `.top` and `.itp`.

        A protonation state is one sample from an ensemble, so the pH and the
        seed it was drawn with are what make a topology reproducible. Neither
        is recoverable from the file: a zero charge is equally consistent with
        pH 2 and with no pH at all.
        """
        lines = []
        if self.linked_units is not None:
            lines.append(f"Linked structure: {self.linked_units.describe()}, seed {self.config.seed}")
        if self.config.pH is None:
            return lines
        lines.append(
            f"Protonation sampled at pH {self.config.pH:g}, seed {self.config.seed}"
        )
        counts = getattr(self.composition, "ionized_counts", None)
        if counts:
            placed = ", ".join(f"{g} {n}" for g, n in sorted(counts.items()))
            lines.append(f"Ionized in this sample: {placed}")
        lines.append(
            "One sample of the ensemble at this pH, not the ensemble mean -- "
            "vary the seed for replicates"
        )
        return lines

    def print_summary(self):
        """Print summary of generated structure."""
        if self.composition is None:
            print("No structure generated yet. Call generate() first.")
            return

        print("\n" + "=" * 60)
        print("BIOCHAR STRUCTURE SUMMARY")
        print("=" * 60)
        print("\nComposition:")
        print(f"  Carbons:     {self.composition.num_carbons}")
        print(f"  Hydrogens:   {self.composition.num_hydrogens}")
        print(f"  Oxygens:     {self.composition.num_oxygens}")
        if self.composition.num_nitrogens:
            print(f"  Nitrogens:   {self.composition.num_nitrogens}")
        if self.composition.num_sulfurs:
            print(f"  Sulfurs:     {self.composition.num_sulfurs}")
        print(f"  Formula:     {self.composition.molecular_formula}")
        print(f"  MW:          {self.composition.molecular_weight:.1f} g/mol")
        if self.config.pH is not None:
            print("\nProtonation:")
            print(f"  pH:          {self.config.pH:.2f}")
            print(f"  Net charge:  {self.composition.net_charge:+d} e")
            for group, n in sorted(self.composition.titratable_counts.items()):
                ionized = self.composition.ionized_counts.get(group, 0)
                print(f"  {group + ':':13s}{ionized}/{n} ionized")
            if self.composition.net_charge != 0:
                print("  (gmx genion -neutral balances this at solvation)")
        print("\nRatios:")
        print(f"  H/C ratio:   {self.composition.H_C_ratio:.3f} (target: {self.config.H_C_ratio:.3f})")
        print(f"  O/C ratio:   {self.composition.O_C_ratio:.3f} (target: {self.config.O_C_ratio:.3f})")
        if self.composition.num_nitrogens:
            print(f"  N/C ratio:   {self.composition.N_C_ratio:.3f}")
        if self.composition.num_sulfurs:
            print(f"  S/C ratio:   {self.composition.S_C_ratio:.3f}")
        if (self.composition.num_pyridinic or self.composition.num_pyrrolic
                or self.composition.num_graphitic):
            print("\nRing Nitrogen:")
            if self.composition.num_pyridinic:
                print(f"  Pyridinic:   {self.composition.num_pyridinic}")
            if self.composition.num_pyrrolic:
                print(f"  Pyrrolic:    {self.composition.num_pyrrolic}")
            if self.composition.num_graphitic:
                print(f"  Graphitic:   {self.composition.num_graphitic}")
        print("\nFunctional Groups:")
        if self.composition.functional_groups:
            for group_name, count in self.composition.functional_groups.items():
                print(f"  {group_name}: {count}")
        else:
            print("  None")

        if self.validation_report:
            print("\nValidation:")
            print(f"  Status: {'VALID' if self.validation_report[0] else 'INVALID'}")
            all_errors = self.validation_report[1]
            valence_errors = [e for e in all_errors if "valence" in e.lower() or "Valence" in e]
            other_errors = [e for e in all_errors if e not in valence_errors]
            if valence_errors:
                print(f"  Valence Issues: {len(valence_errors)}")
                for error in valence_errors[:3]:
                    print(f"    - {error}")
            else:
                print("  Valence Issues: 0")
            if other_errors:
                print(f"  Other Errors: {len(other_errors)}")
                for error in other_errors[:3]:
                    print(f"    - {error}")
            if self.validation_report[2]:
                print(f"  Warnings: {len(self.validation_report[2])}")
                for warning in self.validation_report[2][:3]:
                    print(f"    - {warning}")

        print("=" * 60 + "\n")

    # Private methods

    def _generate_carbon_skeleton(self):
        """Generate the carbon skeleton (aromatic core + optional aliphatic C)."""
        self.linked_units = None
        if self.config.linked_units:
            plan = plan_units(self.config, random.Random(self.config.seed))
            if plan.num_units > 1:
                # One base seed for every linked stage. seed=None means a fresh
                # structure per call, as on the fused path, so draw it here
                # rather than letting each stage fall back to a fixed value.
                self._linked_seed = (
                    random.SystemRandom().randrange(2 ** 31)
                    if self.config.seed is None else int(self.config.seed)
                )
                r = None if self.config.functional_groups else self.config.H_C_ratio
                skeleton, self.linked_units = build_linked_skeleton(
                    self.config, plan, self._aliphatic_enabled(r), self._linked_seed
                )
                logger.info("Linked skeleton: %s", self.linked_units.describe())
                return skeleton
        assembler = PAHAssembler(seed=self.config.seed)
        N = self.config.target_num_carbons
        # Explicit functional groups specify a precise composition and carry
        # their own hydrogen (e.g. amino -NH2), so we do NOT auto-shape the
        # skeleton for H/C in that case -- neither elongation nor aliphatic
        # decoration would otherwise double-count that hydrogen and overshoot.
        # H/C is then whatever the aromatic edges + requested groups produce.
        r = None if self.config.functional_groups else self.config.H_C_ratio

        # First pass: aromatic skeleton at full size (with H/C-aware elongation).
        skeleton = assembler.generate(
            N,
            defect_fraction=self.config.defect_fraction,
            heptagon_fraction=self.config.heptagon_fraction,
            target_h_c=r,
        )

        # If the requested H/C still exceeds what this (already elongated)
        # aromatic flake can carry, rebuild a smaller aromatic core and decorate
        # it with pendant methyls so the total carbon count stays ~N while the
        # H/C reaches the target.  The methyl count is computed from the core
        # that is *actually built* so elongation and methyls never double-count.
        if self._aliphatic_enabled(r) and r > _edge_carbon_fraction(skeleton.mol) + 0.02:
            ceiling = _edge_carbon_fraction(skeleton.mol)
            # Estimate the aromatic/aliphatic split (keeps total C ~ N).  Each
            # methyl adds +1 C and a net +2 H, so H/C = ceiling + n(2-ceiling)/N.
            def split(c):
                n = round(N * (r - c) / (2.0 - c))
                return max(1, min(n, N - 6))

            # Build a COMPACT core (target_h_c=None): the methyls, not
            # elongation, supply the extra hydrogen here, so the aromatic
            # ceiling stays low and the exact methyl count below cannot be
            # driven negative by an already-elongated (overshooting) core.
            #
            # The split above is estimated from the first pass's ceiling, but
            # that pass was elongated and the core is compact, with a lower
            # ceiling -- so the core needs more methyls than were budgeted, and
            # total C overshot N. Re-plan from
            # the ceiling of the core actually built until the total is within
            # size_tolerance, keeping the build that best meets the H/C and then N.
            n_est = split(ceiling)
            best = None
            tried = set()
            for _ in range(_CORE_SPLIT_ITERATIONS):
                core = assembler.generate(
                    N - n_est,
                    defect_fraction=self.config.defect_fraction,
                    heptagon_fraction=self.config.heptagon_fraction,
                    target_h_c=None,
                )
                # Exact methyl count for the core as built.  A saturated
                # aromatic core carries E edge hydrogens; attaching m methyls
                # gives H = E + 2m and C = Ca + m, so H/C = target when
                #   m = (target*Ca - E) / (2 - target).
                Ca, E = self._aromatic_core_stats(core.mol)
                m = round((r * Ca - E) / (2.0 - r))
                m = max(0, min(m, E))
                # Composition first, size second: a smaller core can cap m at
                # its edge-site count and fall short of the H/C, which is the
                # request this path exists for. Misses inside 0.02 count as hits.
                hc = (E + 2 * m) / (Ca + m) if Ca + m else 0.0
                score = (max(0.0, abs(hc - r) - 0.02), abs(Ca + m - N))
                if best is None or score < best[0]:
                    best = (score, core, Ca, E, m, hc)
                # Re-planning is for the overshoot only. A plan within
                # tolerance is kept as it was (a different split changes which
                # structure the comparison below picks), and so is one that
                # falls short: that happens when m is capped at the edge sites,
                # and a smaller core would only fall shorter.
                if Ca + m - N <= self.config.size_tolerance * N:
                    break
                # The core build is deterministic for a given size, so a
                # plan seen before would only rebuild the same core.
                tried.add(n_est)
                n_est = split(E / Ca) if Ca else n_est
                if n_est in tried:
                    break
            _, core, Ca, E, m, cand_hc = best
            if m > 0:
                rng = random.Random(self.config.seed)
                decorated = attach_aliphatic_carbons(core.mol, m, rng)
                cand = PAHAssembler._make_skeleton(decorated)
            else:
                # Core already meets/exceeds the target on its own.
                cand = core
            # Keep the decorated candidate only if it lands closer to the target
            # than the first-pass aromatic skeleton.  A benzene-grown "compact"
            # core is not always minimally condensed, so it can overshoot; in
            # that case the first-pass skeleton is the better structure.
            if abs(cand_hc - r) < abs(ceiling - r):
                skeleton = cand

        # Validate skeleton
        valid, errors = SkeletonValidator.validate(skeleton)
        if not valid:
            logger.warning("Skeleton validation issues: %s", errors)

        return skeleton

    def _aliphatic_enabled(self, target_h_c: Optional[float]) -> bool:
        """Whether aliphatic (sp3) decoration may be used to raise H/C."""
        if (not self.config.allow_aliphatic or target_h_c is None
                or self.config.target_num_carbons < 8):
            return False
        # Honour an explicit request for a (near) fully aromatic structure.
        ap = self.config.aromaticity_percent
        if ap is not None and ap >= 99.0:
            return False
        return True

    @staticmethod
    def _aromatic_core_stats(mol: Chem.Mol) -> Tuple[int, int]:
        """Return (total carbons, aromatic edge carbons with a free valence)."""
        Ca = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        E = sum(
            1 for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6 and a.GetIsAromatic() and a.GetDegree() < 3
        )
        return Ca, E

    def _assign_oxygens(self, mol: Chem.Mol) -> Tuple[Chem.Mol, CompositionResult]:
        """Assign oxygen atoms. Returns (mol, CompositionResult) with placed/requested counts."""
        if self.linked_units is not None:
            mol, comp = self._assign_oxygens_linked(mol)
        else:
            assigner = OxygenAssigner(seed=self.config.seed,
                                      max_ether_span=self.config.max_ether_span)
            mol, comp = assigner.assign_oxygens(
                mol,
                self.config.O_C_ratio,
                O_C_tolerance=self.config.O_C_tolerance,
                functional_group_preference=self.config.functional_groups,
                allow_aliphatic_oxygen=self.config.allow_aliphatic_oxygen,
            )

        # Strict mode fails on a shortfall, not only on total failure. Placing 6
        # of 40 requested ether bridges is not the structure the caller asked
        # for, and returning it as a success hides the gap exactly where a sweep
        # records the row as strict_pass.
        #
        # Scoped to *explicitly requested* groups. When functional_groups is None
        # the counts in requested_counts were derived from O_C_ratio, and there
        # the ratio is the target while the count is an implementation detail of
        # reaching it -- that case is judged by composition validation against
        # O_C_tolerance, which is the tolerance the caller actually set.
        if self.config.strict and self.config.functional_groups and comp.requested_counts:
            shortfall = [
                (g, n, comp.placed_counts.get(g, 0))
                for g, n in comp.requested_counts.items()
                if n > 0 and comp.placed_counts.get(g, 0) < n
            ]
            if shortfall:
                detail = ", ".join(
                    f"{g} {placed}/{requested}" for g, requested, placed in shortfall
                )
                raise ValidationError(
                    f"Strict mode: functional groups were not fully placed "
                    f"({detail}) — the skeleton has too few suitable sites. "
                    f"Reduce the requested counts, raise target_num_carbons, or "
                    f"set strict=False to accept the shortfall."
                )

        return mol, comp

    _LINKED_OXYGEN_ATTEMPTS = 5

    def _assign_oxygens_linked(self, mol: Chem.Mol) -> Tuple[Chem.Mol, CompositionResult]:
        """Oxygen placement on a linked skeleton.

        OxygenAssigner sizes its oxygen count from num_carbons x O_C_ratio and
        ignores oxygen already present, so the linker oxygens are taken off the
        target here; the realised O/C, which counts every O, then lands on it.
        A bridging group placed across two units would merge their ring systems
        and close a cycle in the unit tree, so such a placement is redrawn on a
        derived seed.
        """
        n_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        o_c = max(0.0, self.config.O_C_ratio - self.linked_units.num_linker_oxygens / max(n_c, 1))
        base = self._linked_seed
        for attempt in range(self._LINKED_OXYGEN_ATTEMPTS):
            assigner = OxygenAssigner(seed=base + 7919 * attempt,
                                      max_ether_span=self.config.max_ether_span)
            out, comp = assigner.assign_oxygens(
                mol,
                o_c,
                O_C_tolerance=self.config.O_C_tolerance,
                functional_group_preference=self.config.functional_groups,
                allow_aliphatic_oxygen=self.config.allow_aliphatic_oxygen,
            )
            if not unit_bridging_atoms(out):
                return out, comp
        message = (
            "Linked units: every oxygen placement bridged two units (a bridging "
            "group across a link merges their ring systems); request fewer "
            "bridging groups or fewer units."
        )
        if self.config.strict:
            raise ValidationError("Strict mode: " + message)
        logger.warning(message)
        return out, comp

    def _substitute_nitrogens(
        self, mol: Chem.Mol, comp_result: CompositionResult
    ) -> Chem.Mol:
        """
        Substitute ring carbons with nitrogen (pyridinic/pyrrolic/graphitic).

        Runs after oxygen placement and before hydrogen saturation so that the
        new valence requirements are satisfied by :class:`HydrogenAssigner`.
        Updates *comp_result* ring-N census fields in-place.
        """
        n_py = self.config.num_pyridinic
        n_pr = self.config.num_pyrrolic
        n_gr = self.config.num_graphitic
        if n_py <= 0 and n_pr <= 0 and n_gr <= 0:
            return mol

        substitutor = NitrogenSubstitutor(seed=self.config.seed)
        mol = substitutor.substitute(
            mol,
            n_pyridinic=n_py,
            n_pyrrolic=n_pr,
            n_graphitic=n_gr,
        )

        comp_result.num_pyridinic = substitutor.placed_pyridinic
        comp_result.num_pyrrolic = substitutor.placed_pyrrolic
        comp_result.num_graphitic = substitutor.placed_graphitic
        comp_result.placed_counts["pyridinic"] = substitutor.placed_pyridinic
        comp_result.placed_counts["pyrrolic"] = substitutor.placed_pyrrolic
        comp_result.placed_counts["graphitic"] = substitutor.placed_graphitic
        comp_result.requested_counts["pyridinic"] = n_py
        comp_result.requested_counts["pyrrolic"] = n_pr
        comp_result.requested_counts["graphitic"] = n_gr

        return mol

    def _assign_protonation(
        self, mol: Chem.Mol, comp_result: CompositionResult
    ) -> Chem.Mol:
        """
        Set protonation states from ``config.pH``, updating *comp_result*.

        A no-op when ``pH`` is None, which is what keeps the default path
        byte-for-byte identical to the pre-pH generator.
        """
        if self.config.pH is None:
            return mol

        from .protonation import ProtonationAssigner

        assigner = ProtonationAssigner(seed=self.config.seed)
        mol, _ = assigner.assign(mol, pH=self.config.pH, result=comp_result)

        if comp_result.net_charge != 0:
            logger.info(
                "pH %.2f → net charge %+d e (%s). The topology carries this "
                "charge; `gmx genion -neutral` adds counterions at solvation.",
                self.config.pH,
                comp_result.net_charge,
                ", ".join(
                    f"{n}×{g}" for g, n in sorted(comp_result.ionized_counts.items())
                ) or "no ionized groups",
            )
        return mol

    def _assign_hydrogens(self, mol: Chem.Mol, comp_result: CompositionResult) -> Chem.Mol:
        """Assign hydrogen atoms, updating *comp_result* in-place."""
        assigner = HydrogenAssigner(seed=self.config.seed)
        mol, _ = assigner.assign_hydrogens(
            mol,
            self.config.H_C_ratio,
            H_C_tolerance=self.config.H_C_tolerance,
            result=comp_result,
        )

        # Validate
        valid, errors = HeteroatomValidator.validate_ratios(
            comp_result,
            self.config.H_C_ratio,
            self.config.O_C_ratio,
            self.config.H_C_tolerance,
            self.config.O_C_tolerance,
        )
        if not valid:
            logger.warning("Composition validation issues: %s", errors)

        return mol

    def _generate_geometry(self, mol: Chem.Mol) -> Tuple[Chem.Mol, np.ndarray]:
        """Generate 3D coordinates with clash resolution."""
        if self.linked_units is not None:
            # Flat-sheet embedding would lay the units out coplanar; the linked
            # embedder keeps each unit flat and rotates the links.
            mol, coords = embed_linked_units(mol, seed=self._linked_seed)
            valid, errors = GeometryValidator.validate_geometry(mol, coords)
            if not valid:
                print("Warning: Geometry validation issues:")
                for error in errors[:3]:
                    print(f"  - {error}")
            return mol, coords
        generator = CoordinateGenerator(seed=self.config.seed)
        mol, coords = generator.generate_3d_coordinates(
            mol,
            force_aromatic_planarity=True,
        )

        # Validate and detect clashes
        valid, errors = GeometryValidator.validate_geometry(mol, coords)

        # Resolve clashes if found.
        #
        # When hex-lattice positions are used (generator.used_hex_lattice is
        # True), the generic resolver and FF passes do not run: the resolver
        # displaces atoms pairwise with no notion of a ring, and would shatter
        # the lattice.  What a flat lattice gets wrong is local -- helicene-type
        # edges.  At a cove (edge carbons 2.46 Å apart, four bonds apart) both
        # substituents point into the same vacant site and land 0.1-0.6 Å
        # apart, and the cove carbons sit inside the 2.50 Å C...C floor; a fjord
        # is worse.  relieve_lattice_crowding tilts those substituents out of
        # plane and lets the coves open under a tethered UFF pass, keeping the
        # sheet's bond lengths and overall flatness (rq-8e5f5dc1).
        if generator.used_hex_lattice:
            coords = generator.relieve_lattice_crowding(mol, coords)
            valid, errors = GeometryValidator.validate_geometry(mol, coords)

            # A sheet decorated with off-lattice carbons (methyls, -CH2OH)
            # can carry a substituent on *both* carbons of a cove: two methyl
            # carbons ~0.2 Å apart, which a rigid tilt cannot fully clear.  A
            # free embedding -- every atom movable -- often can, so when the
            # relief leaves a clash, embed freely too and keep whichever has
            # fewer geometry errors.  Neither path is better everywhere:
            # lightly decorated sheets clear on the lattice and not freely
            # (issue #66), heavily decorated ones the other way round.
            # Undecorated sheets never take this branch.
            if generator.offlattice_carbons and any("Steric clash" in e for e in errors):
                free = CoordinateGenerator(seed=self.config.seed)
                free_mol, free_coords = free.generate_3d_coordinates(
                    Chem.Mol(mol), force_aromatic_planarity=True, use_hex_lattice=False,
                )
                free_coords, free_valid, free_errors = self._refine_free_embedding(
                    free, free_mol, free_coords
                )
                # validate_geometry stops at the first NaN with a single
                # error, so a non-finite embedding must not win on count.
                if np.all(np.isfinite(free_coords)) and len(free_errors) < len(errors):
                    mol, coords, valid, errors = free_mol, free_coords, free_valid, free_errors
        else:
            coords, valid, errors = self._refine_free_embedding(generator, mol, coords)

        if not valid:
            print("Warning: Geometry validation issues:")
            for error in errors[:3]:
                print(f"  - {error}")

        # Measure planarity
        planarity, assessment = GeometryValidator.measure_ring_planarity(mol, coords)
        print(f"  Ring planarity: {assessment} (deviation: {planarity:.3f} Å)")

        return mol, coords

    def _refine_free_embedding(
        self, generator: CoordinateGenerator, mol: Chem.Mol, coords: np.ndarray
    ) -> Tuple[np.ndarray, bool, List[str]]:
        """Clash resolution and force-field refinement for an embedding that
        did not come from the hex lattice; returns ``(coords, valid, errors)``."""
        valid, errors = GeometryValidator.validate_geometry(mol, coords)
        steric_clashes = [e for e in errors if "Steric clash" in e]

        # The force-field pass is NOT gated on clashes.  A clean clash report
        # does not mean the embedding is strain-free: ETKDG can leave a
        # compressed aromatic bond (e.g. 1.16 Å where 1.40 Å is expected) that
        # only the FF refinement relaxes.  While every high-oxygen structure
        # carried at least one (false) H-bond "clash", this pass always ran and
        # the coupling was invisible; once H-bonds stopped being reported as
        # clashes it would have silently skipped refinement altogether.
        # Clash *resolution* still runs only when there is a clash to resolve.

        # First pass: iterative clash resolution
        if steric_clashes:
            coords = ClashResolver.resolve_clashes(
                mol, coords, max_iterations=15, displacement_step=0.15, use_vdw_radii=True
            )

        # Second pass: force-field refinement.
        coords, _ = generator.validate_and_relax(mol, coords, max_iterations=200)

        # Third pass: final clash resolution if needed
        valid_check, errors_check = GeometryValidator.validate_geometry(mol, coords)
        steric_clashes_after_ff = [e for e in errors_check if "Steric clash" in e]

        if steric_clashes_after_ff:
            coords = ClashResolver.resolve_clashes(
                mol, coords, max_iterations=10, displacement_step=0.08, use_vdw_radii=True
            )
            coords, _ = generator.validate_and_relax(mol, coords, max_iterations=200)

        # Final validation
        valid, errors = GeometryValidator.validate_geometry(mol, coords)
        if steric_clashes:
            steric_clashes_after = [e for e in errors if "Steric clash" in e]
            clashes_resolved = len(steric_clashes) - len(steric_clashes_after)
            print(f"  Clash resolution: {clashes_resolved}/{len(steric_clashes)} clashes resolved")
        return coords, valid, errors

    def _assign_opls_properties(self, mol: Chem.Mol, coords: np.ndarray):
        """Assign OPLS-AA atom types and charges.

        Atom types are always assigned from the OPLS-AA typer.  Partial charges
        come from the configured ``charge_method``: the static OPLS table
        (``"opls"``, default), the ML refiner (``"ml"``), or LigParGen-style QM
        charges (``"qm"``), which need the 3D ``coords`` for the AM1 calculation.
        """
        typer = AtomTyper()
        self.atom_types = typer.assign_atom_types(mol)

        charger = ChargeAssigner()
        self.charges = charger.assign_charges(mol, self.atom_types)

        if self.config.charge_method == "ml":
            from ..charges.ml_charges import MLChargeRefinement
            refiner = MLChargeRefinement()
            self.charges = refiner.refine(mol, self.atom_types)
            logger.info("ML charge refinement applied.")
        elif self.config.charge_method == "qm":
            from ..charges.qm_charges import QMChargeAssigner
            assigner = QMChargeAssigner()
            self.charges = assigner.assign(mol, coords, self.atom_types)
            logger.info("QM (1.14*CM1A) charge assignment applied.")

        logger.info(
            "Atom types assigned: %d unique types, total charge: %.3f e",
            len(set(self.atom_types.values())),
            sum(self.charges.values()),
        )

    def _validate(
        self, mol: Chem.Mol, composition: CompositionResult, coords: np.ndarray
    ):
        """Validate complete structure; raises ValidationError in strict mode."""
        is_valid, errors, warnings, metrics = ValidationEngine.validate_complete(
            mol,
            composition,
            coords,
            self.config.H_C_ratio,
            self.config.O_C_ratio,
            self.config.aromaticity_percent,
            self.config.H_C_tolerance,
            self.config.O_C_tolerance,
        )

        self.validation_report = (is_valid, errors, warnings, metrics)

        if not is_valid:
            logger.error("Validation FAILED with %d error(s): %s", len(errors), "; ".join(errors))
            if self.config.strict:
                raise ValidationError(
                    f"Strict mode: structure failed validation with "
                    f"{len(errors)} error(s): " + "; ".join(errors)
                )
        else:
            logger.info("Validation PASSED")
            if warnings:
                logger.warning("%d validation warning(s): %s", len(warnings), "; ".join(warnings))


def generate_biochar(
    target_num_carbons: int = 50,
    H_C_ratio: Optional[float] = None,
    O_C_ratio: Optional[float] = None,
    aromaticity_percent: Optional[float] = None,
    functional_groups: Optional[Dict[str, int]] = None,
    pH: Optional[float] = None,
    defect_fraction: float = 0.0,
    heptagon_fraction: float = 0.0,
    max_ether_span: Optional[int] = None,
    num_pyridinic: int = 0,
    num_pyrrolic: int = 0,
    num_graphitic: int = 0,
    charge_method: str = "opls",
    temperature: Optional[float] = None,
    feedstock: Optional[str] = None,
    output_directory: str = ".",
    basename: str = "biochar",
    molecule_name: str = "BC",
    seed: Optional[int] = None,
    write_files: bool = True,
    linked_units: bool = False,
    num_units: Optional[int] = None,
    linker_mix: Optional[Dict[str, float]] = None,
) -> "BiocharResult":
    """
    Convenience function to generate and export biochar in one call.

    Args:
        target_num_carbons: Target number of carbon atoms.
        H_C_ratio: Target hydrogen-to-carbon ratio.
        O_C_ratio: Target oxygen-to-carbon ratio.  Used to determine total
            oxygen when *functional_groups* is None.
        aromaticity_percent: Target aromaticity percentage.
        functional_groups: Dict mapping functional group name → exact count,
            e.g. ``{"phenolic": 3, "carboxyl": 1}``.
        pH: Environmental pH.  ``None`` (default) leaves every group neutral.
            When set, each titratable site is independently ionized with its
            Henderson-Hasselbalch probability and the structure carries a real
            net charge (``result.composition.net_charge``), which
            ``gmx genion -neutral`` balances at solvation time.

            Supported groups:

            * ``phenolic``  — aromatic C–OH             (1 O per group)
            * ``hydroxyl``  — same as phenolic for pure aromatic PAH (1 O)
            * ``carboxyl``  — aromatic C–C(=O)(OH)      (2 O per group)
            * ``ether``     — aromatic C–O–C bridge     (1 O per group)
            * ``carbonyl``  — not supported; substituted with phenolic
            * ``quinone``   — not supported; substituted with phenolic
            * ``lactone``   — not supported; substituted with phenolic

            If ``None`` (default), the total oxygen count is derived from
            *O_C_ratio* and placed as phenolic (–OH) groups.
        defect_fraction: Probability [0, 1) that each ring added during
            skeleton growth is a 5-membered (pentagon) ring rather than a
            hexagon.  0.0 (default) = pure hexagonal PAH.  Values ~0.1–0.2
            introduce realistic topological disorder.
        heptagon_fraction: Probability [0, 1) that each ring added is a
            7-membered (heptagon) ring, adding the saddle curvature of
            non-graphitizing carbons.  0.0 (default) = no heptagons.  Wood
            et al. 2024 ratios ≈ defect_fraction 0.154 + heptagon_fraction 0.077.
        max_ether_span: Maximum number of C–C bonds between the two ring
            carbons bridged by each ether oxygen.  Controls the ring size of
            the C–O–C bridge (ring size = max_ether_span + 2).
            ``None`` (default) uses :attr:`GeneratorConfig.max_ether_span`
            default of 3 (5-membered furan/benzofuran-like ring — always flat).
            Use 4 for pyran/chromene-like (6-membered) or 5 for 7-membered;
            larger values risk cross-sheet bridges that fold the molecule.
        output_directory: Output directory for GROMACS files.
        basename: Base filename for output files.
        molecule_name: Residue name (max 5 chars). Suggested: BC400, BC600,
            BCH05, BCO10.
        seed: Random seed for reproducibility.
        write_files: If ``True`` (default), write ``.gro``, ``.top``, and
            ``.itp`` files to *output_directory*.  Set to ``False`` to skip
            all disk I/O; the path fields on the returned :class:`BiocharResult`
            will be ``None``.
        linked_units: Build several aromatic units joined by linkers instead
            of one fused sheet (see :attr:`GeneratorConfig.linked_units`).
        num_units: Explicit unit count for a linked structure.
        linker_mix: Relative weights of linker types for a linked structure.

    Returns:
        :class:`BiocharResult` with named fields ``mol``, ``coords``,
        ``composition``, ``gro_path``, ``top_path``, ``itp_path``.
        Supports 5-tuple positional unpacking for backward compatibility.
    """
    config_kwargs: Dict = dict(
        target_num_carbons=target_num_carbons,
        H_C_ratio=H_C_ratio,
        O_C_ratio=O_C_ratio,
        aromaticity_percent=aromaticity_percent,
        functional_groups=functional_groups,
        pH=pH,
        defect_fraction=defect_fraction,
        heptagon_fraction=heptagon_fraction,
        num_pyridinic=num_pyridinic,
        num_pyrrolic=num_pyrrolic,
        num_graphitic=num_graphitic,
        charge_method=charge_method,
        temperature=temperature,
        feedstock=feedstock,
        molecule_name=molecule_name,
        seed=seed,
        linked_units=linked_units,
        num_units=num_units,
        linker_mix=linker_mix,
    )
    if max_ether_span is not None:
        config_kwargs["max_ether_span"] = max_ether_span
    config = GeneratorConfig(**config_kwargs)

    generator = BiocharGenerator(config)
    mol, coords, composition = generator.generate()
    generator.print_summary()

    gro_path = top_path = itp_path = None
    if write_files:
        gro_path, top_path, itp_path = generator.export_gromacs(
            output_directory=output_directory,
            basename=basename,
        )

    return BiocharResult(
        mol=mol,
        coords=coords,
        composition=composition,
        gro_path=gro_path,
        top_path=top_path,
        itp_path=itp_path,
        ring_composition=generator.ring_composition,
        linked_units=generator.linked_units,
    )


def generate_biochar_series(
    configurations: List[Dict],
    output_directory: str = ".",
    create_combined_top: bool = True,
    verbose: bool = True,
    progress_callback: Optional[Callable[[int, int, str], None]] = None,
    on_error: str = "raise",
) -> Dict[str, Tuple[Path, Path, Path]]:
    """
    Generate multiple biochar structures for mixed simulations.

    This function is ideal for creating temperature series, composition series,
    or mixed biochar systems for GROMACS simulations.

    Args:
        configurations: List of configuration dictionaries. Each dict should contain:
            - 'molecule_name' (str, required): Residue name (max 5 chars, e.g., 'BC400')
            - 'target_num_carbons' (int, optional): Default 50
            - 'H_C_ratio' (float, optional): Default 0.5
            - 'O_C_ratio' (float, optional): Default 0.1
            - 'aromaticity_percent' (float, optional): Default 90.0
            - 'seed' (int, optional): For reproducibility
            - 'functional_groups' (dict, optional): e.g. {"phenolic": 3, "carboxyl": 1}

        output_directory: Output directory for all files
        create_combined_top: If True, generate a combined topology for all structures
        verbose: If True, print progress information
        progress_callback: Optional callable invoked after each successful
            generation as ``callback(completed, total, molecule_name)``.
            Useful for GUI/web frontends and long overnight runs.
        on_error: What to do when a structure fails generation.
            ``"raise"`` (default) — re-raise the exception immediately.
            ``"skip"`` — skip silently; the molecule is omitted from results.
            ``"warn"`` — log a warning and skip.

    Returns:
        Dictionary mapping molecule_name -> (gro_path, top_path, itp_path)

    Example:
        >>> configs = [
        ...     {'molecule_name': 'BC400', 'H_C_ratio': 0.65, 'O_C_ratio': 0.20},
        ...     {'molecule_name': 'BC600', 'H_C_ratio': 0.55, 'O_C_ratio': 0.12},
        ...     {'molecule_name': 'BC800', 'H_C_ratio': 0.40, 'O_C_ratio': 0.05},
        ... ]
        >>> results = generate_biochar_series(configs, output_directory='output')
    """
    # Create output directory
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)

    results = {}
    itp_files = {}
    molecule_names = []

    if verbose:
        print("\n" + "=" * 70)
        print(f"BATCH BIOCHAR GENERATION - {len(configurations)} structures")
        print("=" * 70 + "\n")

    for i, config in enumerate(configurations, 1):
        # Extract molecule name first (needed for progress/error messages).
        # Missing name always raises immediately — on_error does not apply to
        # misconfigured input.
        molecule_name = config.get("molecule_name")
        if not molecule_name:
            raise ValueError(f"Configuration {i} missing required 'molecule_name'")

        target_carbons = config.get("target_num_carbons", 50)
        h_c = config.get("H_C_ratio", 0.5)
        o_c = config.get("O_C_ratio", 0.1)
        arom = config.get("aromaticity_percent", 90.0)
        seed = config.get("seed", None)
        fg = config.get("functional_groups", None)

        if verbose:
            print(f"[{i}/{len(configurations)}] Generating {molecule_name}...")
            print(f"  Parameters: {target_carbons}C, H/C={h_c:.2f}, O/C={o_c:.2f}")

        try:
            if len(molecule_name) > 5:
                raise ValueError(
                    f"molecule_name '{molecule_name}' exceeds 5 character limit "
                    f"(GROMACS .gro format requirement)"
                )
            result = generate_biochar(
                target_num_carbons=target_carbons,
                H_C_ratio=h_c,
                O_C_ratio=o_c,
                aromaticity_percent=arom,
                functional_groups=fg,
                output_directory=str(output_path),
                basename=molecule_name.lower(),
                molecule_name=molecule_name,
                seed=seed,
            )
            gro_path, top_path, itp_path = result.gro_path, result.top_path, result.itp_path

            results[molecule_name] = (gro_path, top_path, itp_path)
            itp_files[molecule_name] = itp_path
            molecule_names.append(molecule_name)

            if verbose:
                print(f"  ✓ Successfully generated {gro_path.name}\n")

            if progress_callback is not None:
                progress_callback(len(results), len(configurations), molecule_name)

        except Exception as e:
            if verbose:
                print(f"  ✗ Failed: {e}\n")
            if on_error == "raise":
                raise
            elif on_error == "warn":
                logger.warning(
                    "Skipping '%s' after error: %s", molecule_name, e
                )
            # on_error == "skip": silently continue

    # Generate combined topology file if requested
    if create_combined_top and len(results) > 1:
        combined_top_path = output_path / "combined.top"
        _create_combined_topology(
            itp_files, molecule_names, combined_top_path, verbose
        )

    if verbose:
        print("=" * 70)
        print(f"BATCH GENERATION COMPLETE - {len(results)} structures generated")
        print("=" * 70 + "\n")

        if create_combined_top and len(results) > 1:
            print(f"Combined topology: {combined_top_path}")
            print("\nTo use in GROMACS simulations:")
            print("  gmx grompp -f md.mdp -c combined.gro -p combined.top -o topol.tpr")
            print()

    return results


def _create_combined_topology(
    itp_files: Dict[str, Path],
    molecule_names: List[str],
    output_path: Path,
    verbose: bool = True,
) -> Path:
    """
    Create a combined topology file for mixed biochar simulations.

    Args:
        itp_files: Dictionary mapping molecule_name -> itp_path
        molecule_names: List of molecule names in order
        output_path: Path to write combined .top file
        verbose: If True, print status information

    Returns:
        Path to combined topology file
    """
    if verbose:
        print("\nCreating combined topology file...")

    with open(output_path, "w") as f:
        # Header
        f.write("; Combined topology for mixed biochar simulation\n")
        f.write("; Auto-generated by Biochar Generator\n")
        f.write(f"; Contains {len(molecule_names)} biochar structures\n\n")

        # Include forcefield
        f.write('#include "oplsaa.ff/forcefield.itp"\n\n')

        # Include all molecule topologies
        f.write("; Molecule definitions\n")
        for mol_name in molecule_names:
            itp_file = itp_files[mol_name]
            f.write(f'#include "{itp_file.name}"\n')

        f.write("\n")

        # System section
        f.write("[ system ]\n")
        f.write("Mixed Biochar System\n\n")

        # Molecules section
        f.write("[ molecules ]\n")
        f.write("; Name            #molecules\n")
        for mol_name in molecule_names:
            f.write(f"{mol_name:20s} 1\n")

    if verbose:
        print(f"✓ Combined topology written to {output_path.name}")

    return output_path
