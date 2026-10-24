"""
OPLS-AA Force Field Constants

Atom types, charges, and masses for OPLS-AA force field.
Reference: Jorgensen, W. L., et al. JACS 118.45 (1996): 11225-11236.
"""

from dataclasses import dataclass


# OPLS-AA Atom Types for Biochar Systems
# Format: atom_type: (description, mass_amu, default_charge)

OPLS_ATOM_TYPES = {
    # Aromatic carbons and hydrogens
    "CA": ("Aromatic carbon", 12.01, -0.115),
    "HA": ("Aromatic hydrogen", 1.008, 0.115),

    # Aliphatic carbons
    "CT": ("Aliphatic C (sp3)", 12.01, -0.18),
    "HC": ("H on aliphatic C", 1.008, 0.06),

    # Oxygens
    "OC": ("Carbonyl oxygen", 15.999, -0.56),
    "OH": ("Hydroxyl oxygen", 15.999, -0.661),
    "OS": ("Ether oxygen", 15.999, -0.322),
    "OW": ("Water oxygen", 15.999, -0.820),

    # Hydroxyls and carboxylic acid
    "HO": ("H in hydroxyl", 1.008, 0.436),
    "C": ("Carboxylic acid carbonyl C", 12.01, 0.620),
    "O": ("Carboxylic acid carbonyl O", 15.999, -0.540),
    "OH2": ("Hydroxyl on carboxylic acid", 15.999, -0.540),
    "HO2": ("H on carboxylic hydroxyl", 1.008, 0.436),

    # Nitrogen
    "NA": ("Aromatic amine nitrogen (aniline Ar-NH2)", 14.007, -0.60),
    "HNA": ("H on aromatic amine nitrogen", 1.008, 0.30),

    # Sulfur
    "SH_": ("Aromatic thiol sulfur (Ar-SH)", 32.06, -0.39),    # opls_734
    "HSH": ("H on thiol sulfur", 1.008, 0.21),                 # opls_204
    "SS": ("Thioether sulfur bridging two aryl C (Ar-S-Ar)", 32.06, -0.16),  # opls_222
    # Ring-substituting nitrogen (biochar N-doping)
    "NPY": ("Pyridinic N (substituted into 6-ring, no H)", 14.007, -0.36),
    "NPR": ("Pyrrolic N (substituted into 5-ring, with H)", 14.007, -0.52),
    "NGR": ("Graphitic/quaternary N (interior 6-ring, no H)", 14.007, 0.02),
    "HNPR": ("H on pyrrolic nitrogen", 1.008, 0.38),

    # ---- Ionized forms (pH-dependent protonation) -------------------------
    # Charges are the stock OPLS values of the type each one maps to (see
    # GROMACS_OPLS_TYPE_MAP for the mapping and its provenance).  Only
    # carboxylate has a genuine stock OPLS type; the rest are derived from the
    # nearest analog and are flagged as such on every line.
    "CM":  ("Carboxylate carbon (Ar-COO-)", 12.01, 0.700),      # opls_271, exact
    "O2M": ("Carboxylate oxygen (Ar-COO-)", 15.999, -0.800),    # opls_272, exact
    "OM":  ("Phenolate oxygen (Ar-O-)", 15.999, -0.980),        # derived: opls_420 alkoxide
    "SM":  ("Thiophenolate sulfur (Ar-S-)", 32.06, -0.900),     # derived: opls_417 thiolate
    "NPYP": ("Pyridinium N (protonated pyridinic, +1)", 14.007, -0.740),  # derived: opls_379
    "HPYP": ("H on pyridinium N", 1.008, 0.460),                # derived: opls_513 (H on HIP N)
    "NAP": ("Anilinium N (Ar-NH3+)", 14.007, -0.300),           # derived: opls_287 (RNH3+)
    "HNAP": ("H on anilinium N", 1.008, 0.330),                 # derived: opls_290 (H of RNH3+)
}

# Lennard-Jones, bond, and angle parameters intentionally live in oplsaa.ff, not
# here. GROMACS resolves them from the #included forcefield by the opls_XXX name in
# GROMACS_OPLS_TYPE_MAP below, so a hand-copied table here would be dead weight that
# can only drift. A previous table did drift: its values were a mix of AMBER and OPLS
# with no single provenance.
#
# The one exception is SUPPLEMENTARY_ANGLE_PARAMS below, and the rule that keeps it
# from becoming that table again is strict: it may hold ONLY combinations stock
# oplsaa.ff does not define. A value that also exists in the forcefield is
# duplication and will drift; a value that exists nowhere else cannot.
# tests/test_opls_type_map.py enforces both halves -- every emitted combination
# resolves, and nothing here shadows a stock entry.

# Mapping from internal generic type names to GROMACS OPLS-AA opls_XXX names.
# Internal names (CA, HA, etc.) are used throughout the biochar generation pipeline.
# At GROMACS export time, these are translated so the topology is compatible with
# the standard oplsaa.ff forcefield shipped with GROMACS.
#
# When the exported .itp #includes a real oplsaa.ff, only the opls_XXX name reaches
# GROMACS -- mass, charge, LJ and bonded parameters are all resolved from the
# forcefield by that name. A name naming the wrong element therefore silently
# simulates the wrong chemistry. Trailing comments quote the type's description in
# oplsaa.ff/atomtypes.atp verbatim; keep them in sync when editing.
# tests/test_opls_type_map.py checks every entry against an installed oplsaa.ff.
GROMACS_OPLS_TYPE_MAP: dict[str, str] = {
    "CA":  "opls_145",   # aromatic carbon, benzene-type
    "HA":  "opls_146",   # aromatic hydrogen
    "CT":  "opls_135",   # aliphatic sp3 carbon
    "HC":  "opls_140",   # hydrogen on aliphatic carbon
    "OH":  "opls_154",   # phenol / alcohol oxygen
    "HO":  "opls_155",   # phenol / alcohol O-H hydrogen
    "OS":  "opls_467",   # aryl ether oxygen (Ar-O-Ar or Ar-O-R)
    "OC":  "opls_278",   # ketone / aldehyde C=O oxygen
    "C":   "opls_267",   # carboxylic acid carbonyl carbon
    "O":   "opls_269",   # carboxylic acid C=O oxygen
    "OH2": "opls_268",   # carboxylic acid -OH oxygen
    "HO2": "opls_270",   # carboxylic acid -OH hydrogen
    "OW":  "opls_116",   # SPC/E water oxygen
    "NA":  "opls_900",   # "N primary   amines" -- aniline Ar-NH2 nitrogen
    "HNA": "opls_909",   # "H(N)   primary   amines" -- H on that nitrogen
    "SH_": "opls_734",   # "all-atom S: thiophenol (HS is #204)" -- Ar-SH sulfur
    "HSH": "opls_204",   # "all-atom H(S): thiols" -- the HS named by opls_734
    "SS":  "opls_222",   # "S in thioanisoles" -- nearest aryl-S; see note below
    "NPY": "opls_520",   # "N   in pyridine 6-31G*" -- pyridinic ring N
    "NPR": "opls_542",   # "N   in pyrrole" -- pyrrolic ring N
    "HNPR": "opls_545",  # "H1  in pyrrole" -- the pyrrole N-H hydrogen
    "NGR": "opls_379",   # "CytH+ N3 Protonated cytosine." -- a cationic aromatic
                         # ring N. Graphitic N carries a formal +1, so a cationic
                         # aromatic N is a closer analog than the neutral pyridine
                         # N this used to share with NPY. Stock OPLS has no
                         # quaternary aromatic N. See note below.

    # ---- Ionized forms (pH-dependent protonation) -------------------------
    # Only carboxylate has genuine stock OPLS types.  The rest are DERIVED from
    # the nearest available analog -- chemically reasonable, but not validated
    # against QM.  Each line records the analog and why it was chosen.
    "CM":  "opls_271",   # carboxylate C  -- exact ("C in COO- carboxylate")
    "O2M": "opls_272",   # carboxylate O  -- exact ("O in COO- carboxylate")
    "OM":  "opls_420",   # phenolate O    -- DERIVED from "O in CH3O-" (alkoxide).
                         # Stock OPLS has no aryl oxide; alkoxide is the only
                         # deprotonated O available.  Being aliphatic, it likely
                         # overstates the charge on an aryl oxide, whose charge
                         # delocalises into the ring.
    "SM":  "opls_417",   # thiophenolate S -- DERIVED from "S in CH3S-" (thiolate).
                         # Same caveat as OM: stock OPLS has no aryl thiolate.
    "NPYP": "opls_379",  # pyridinium N   -- DERIVED from "CytH+ N3", a protonated
                         # aromatic ring N.  Nearest available cationic aromatic
                         # N; shares the analog with NGR, which is also a
                         # cationic aromatic N.
    "HPYP": "opls_513",  # H on pyridinium N -- DERIVED from "H on N in HIP", the
                         # H on doubly-protonated histidine's cationic ring N.
                         # Structurally the closest stock aromatic N-H+.
    "NAP": "opls_287",   # anilinium N    -- DERIVED from "N (RNH3+)".
                         # Parameterised for alkylammonium, so it does not carry
                         # aniline's ring conjugation.
    "HNAP": "opls_290",  # H on anilinium N -- DERIVED from "H (RNH3+)".
}

# Two entries above are deliberate approximations, not exact matches:
#
#   SS  -> opls_222 is the thioanisole sulfur (Ar-S-CH3). OPLS-AA has no diaryl
#          thioether (Ar-S-Ar) type; opls_222 is the only aryl-attached sulfide S
#          and carries the matching CA-S bond (0.176 nm, "thioanisole" in
#          ffbonded.itp).
#
#          The CA-S bond resolves, but the CA-S-CA angle has no stock angletype
#          -- see SUPPLEMENTARY_ANGLE_PARAMS below, which supplies it.
#
#   NGR -> opls_379 ("CytH+ N3") is a protonated, cationic aromatic ring N.
#          CHOSEN DELIBERATELY, 2026-07-17, and reviewed -- not inherited from
#          whichever branch merged last. Do not "restore" opls_520 on the
#          assumption that this drifted.
#
#          OPLS-AA has no substitutional 3-coordinate aromatic N, so this is an
#          analog either way; a cationic aromatic N is the closer one, because
#          graphitic N carries a formal +1. It previously shared the neutral
#          pyridine N (opls_520) with NPY, which understated exactly that charge
#          -- the reason it moved is the pH work, which makes formal charge real
#          rather than something ChargeAssigner flattened to zero. Element and
#          ring aromaticity are right; the bonded environment is still
#          approximate, and this is not QM-validated. A QM check, or a better
#          analog, would be a legitimate reason to revisit. "It looks like an
#          accident" is not -- it was a decision.

# Angles that stock oplsaa.ff cannot resolve, written inline into [ angles ] so the
# parameters travel with the molecule (the .itp), not with whichever .top includes it.
#
# Keyed by internal atom type, outer two sorted. Values are GROMACS units:
# (theta0_deg, k_kJ/mol/rad^2) -- the same columns as ffbonded.itp [ angletypes ].
#
# Every entry needs a provenance comment naming why the forcefield lacks it and where
# the number came from. Nothing may be added here that oplsaa.ff already defines.
SUPPLEMENTARY_ANGLE_PARAMS: dict[tuple[str, str, str], tuple[float, float]] = {
    # Ar-S-Ar, the diaryl thioether bridge. OPLS-AA has no aryl-S-aryl angle: the
    # closest stock entries are CA-S-CT (thioanisole) and CA-S-CM, both 104.200 /
    # 518.816 and themselves "adjusted from CT-S-CT" per ffbonded.itp. Reusing that
    # value keeps the bridge consistent with the opls_222 sulfur SS already maps to.
    # Approximate, like the SS mapping it accompanies -- not QM-validated.
    ("CA", "SS", "CA"): (104.200, 518.816),
    # Hydroxyl on a ring carbon adjacent to a pyridinic N (3-hydroxypyridine).
    # Resolves to bonded NC-CA-OH, which stock OPLS omits. The value is the phenol
    # angle CA-CA-OH (120.000 / 585.760 in ffbonded.itp) transcribed verbatim: it
    # holds the hydroxyl-on-aromatic-carbon geometry exactly and differs only in the
    # far ring atom (neutral CA vs pyridinic NC). Nearest analog, not QM-validated.
    ("NPY", "CA", "OH"): (120.000, 585.760),
}

# Functional groups definitions
# Each functional group specifies how to add atoms to the carbon skeleton
FUNCTIONAL_GROUPS = {
    "hydroxyl": {
        "description": "Hydroxyl group (-OH)",
        "atoms": [("O", "OH"), ("H", "HO")],  # (atom_type, group_code)
        "connectivity": [(0, "C", 1), (0, 1, 1)],  # (atom_idx, connect_to, bond_type)
        "composition": {"O": 1, "H": 1},
        "O_per_group": 1,
        "H_per_group": 1,
    },
    "carboxyl": {
        "description": "Carboxylic acid group (-COOH)",
        "atoms": [("C", "C"), ("O", "O"), ("O", "OH2"), ("H", "HO2")],
        "connectivity": [
            (0, "C", 2),  # C=O double bond
            (1, "C", 1),  # Single bond C-C
            (2, "C", 1),  # O-H bond
            (3, 2, 1),
        ],
        "composition": {"C": 1, "O": 2, "H": 1},
        "O_per_group": 2,
        "H_per_group": 1,
    },
    "phenolic": {
        "description": "Phenolic group (aromatic -OH)",
        "atoms": [("O", "OH"), ("H", "HO")],
        "connectivity": [(0, "CA", 1), (0, 1, 1)],
        "composition": {"O": 1, "H": 1},
        "O_per_group": 1,
        "H_per_group": 1,
    },
    "aliphatic_hydroxyl": {
        # Primary alcohol on an sp3 carbon: a pendant -CH3 becomes -CH2-OH.
        # This is the dominant O-bearing aliphatic group in low-temperature,
        # cellulose-derived biochar, and unlike phenolic it draws on the sp3
        # (aliphatic) carbons the H/C-shaping stage adds -- which is what lets
        # high-O/C low-aromaticity chars reach their oxygen target at all.
        # Types are identical to phenolic (OH/HO on the O-H); only the carbon
        # it attaches to differs (CT vs CA), so no new OPLS type is needed.
        "description": "Aliphatic hydroxyl (-CH2-OH on sp3 carbon)",
        "atoms": [("O", "OH"), ("H", "HO")],
        "connectivity": [(0, "CT", 1), (0, 1, 1)],
        "composition": {"O": 1, "H": 1},
        "O_per_group": 1,
        "H_per_group": 1,
    },
    "ether": {
        "description": "Ether group (-O-)",
        "atoms": [("O", "OS")],
        "connectivity": [(0, "C", 1), (0, "C", 1)],  # Two C-O bonds
        "composition": {"O": 1},
        "O_per_group": 1,
        "H_per_group": 0,
    },
    "carbonyl": {
        "description": "Carbonyl group (C=O)",
        "atoms": [("O", "OC")],
        "connectivity": [(0, "C", 2)],  # C=O double bond
        "composition": {"O": 1},
        "O_per_group": 1,
        "H_per_group": 0,
    },
    "lactone": {
        "description": "Lactone group (cyclic ester)",
        "atoms": [("C", "C"), ("O", "O"), ("O", "OS")],
        "connectivity": [
            (0, "C", 2),
            (1, "C", 1),
            (2, "C", 1),
        ],
        "composition": {"C": 1, "O": 2},
        "O_per_group": 2,
        "H_per_group": 0,
    },
    "quinone": {
        "description": "Quinone group (C=O with aromatic C)",
        "atoms": [("O", "OC"), ("O", "OC")],
        "connectivity": [(0, "CA", 2), (1, "CA", 2)],
        "composition": {"O": 2},
        "O_per_group": 2,
        "H_per_group": 0,
    },
    "amino": {
        "description": "Amino group (-NH2) on aromatic ring",
        "atoms": [("N", "NA"), ("H", "HNA"), ("H", "HNA")],
        "connectivity": [(0, "CA", 1), (1, 0, 1), (2, 0, 1)],
        "composition": {"N": 1, "H": 2},
        "O_per_group": 0,
        "H_per_group": 2,
    },
    "thiol": {
        "description": "Thiol group (-SH) on aromatic ring",
        "atoms": [("S", "SH_"), ("H", "HSH")],
        "connectivity": [(0, "CA", 1), (1, 0, 1)],
        "composition": {"S": 1, "H": 1},
        "O_per_group": 0,
        "H_per_group": 1,
    },
    "thioether": {
        "description": "Thioether group (-S-) bridging two aromatic carbons",
        "atoms": [("S", "SS")],
        "connectivity": [(0, "C", 1), (0, "C", 1)],  # Two C-S bonds
        "composition": {"S": 1},
        "O_per_group": 0,
        "H_per_group": 0,
    },
}

# ---------------------------------------------------------------------------
# pH-dependent protonation states
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ProtonationState:
    """
    One acid/base equilibrium available to a biochar functional group.

    *pKa* always belongs to the **protonated** member of the pair, following the
    usual convention.  ``kind`` says which side of the equilibrium the charge
    appears on, and is what stops the two blocks being conflated:

        "acidic"  — the neutral form is protonated; it *loses* H+ to become an
                    anion.  Fraction ionized rises with pH.
        "basic"   — the neutral form is deprotonated; it *gains* H+ to become a
                    cation.  Fraction ionized falls with pH.

    ``neutral_type`` / ``ionized_type`` name the OPLS type of the one heteroatom
    whose formal charge changes.  ``h_type`` is the exchangeable hydrogen's OPLS
    type — the H that is removed (acidic) or added (basic).
    """

    pKa: float
    kind: str            # "acidic" | "basic"
    neutral_type: str    # OPLS type of the key heteroatom, neutral form
    ionized_type: str    # OPLS type of the key heteroatom, ionized form
    h_type: str          # OPLS type of the exchangeable hydrogen
    description: str


# Titratable groups, keyed by the group names OxygenAssigner already places
# (plus "pyridinic", which NitrogenSubstitutor substitutes into the ring).
#
# pKa provenance — model-compound values, chosen to sit inside the ranges
# reported for real biochar surfaces by potentiometric and Boehm titration
# (carboxylic 3–6, phenolic 8–11, see refs below):
#
#   carboxyl   4.20  benzoic acid.  Biochar surface carboxyls titrate 3–6;
#                    benzoic acid is the aryl-carboxyl model compound.
#   phenolic   9.50  phenol (9.95) shifted down slightly toward the 9–10 band
#                    reported for biochar phenolic OH, which sits on
#                    electron-poor polyaromatic edges rather than plain benzene.
#   thiol      6.60  thiophenol.  Aryl thiols are far more acidic than alkyl.
#   amino      4.60  anilinium (conjugate acid of aniline).
#   pyridinic  5.20  pyridinium (conjugate acid of pyridine).
#
# Refs: Boehm/potentiometric titration ranges —
#   doi:10.1016/j.scitotenv.2020.142792 (proton uptake vs. pyrolysis temperature)
#   doi:10.1016/j.jcis.2016.01.076      (pKa of graphene-like materials)
#   doi:10.1016/j.carbon.2013.09.048    (limits of the Boehm titration)
#
# NOTE — lactonic groups (pKa 7–9) also titrate in the environmentally
# interesting window but are deliberately absent: "lactone" currently falls back
# to "phenolic" in OxygenAssigner, so there is no lactone to titrate.
#
# NOTE — graphitic N is absent by design.  It is permanently +1 by construction
# (three aromatic ring bonds, pyridinium-like) and does not titrate.
PROTONATION_STATES: dict[str, ProtonationState] = {
    "carboxyl": ProtonationState(
        pKa=4.20,
        kind="acidic",
        neutral_type="OH2",   # -C(=O)OH  hydroxyl oxygen
        ionized_type="O2M",   # -C(=O)O-  carboxylate oxygen
        h_type="HO2",
        description="Ar-COOH <-> Ar-COO- + H+  (benzoic acid, pKa 4.20)",
    ),
    "phenolic": ProtonationState(
        pKa=9.50,
        kind="acidic",
        neutral_type="OH",
        ionized_type="OM",
        h_type="HO",
        description="Ar-OH <-> Ar-O- + H+  (phenol, pKa 9.95 -> 9.50 on PAH edge)",
    ),
    "thiol": ProtonationState(
        pKa=6.60,
        kind="acidic",
        neutral_type="SH_",
        ionized_type="SM",
        h_type="HSH",
        description="Ar-SH <-> Ar-S- + H+  (thiophenol, pKa 6.60)",
    ),
    "amino": ProtonationState(
        pKa=4.60,
        kind="basic",
        neutral_type="NA",
        ionized_type="NAP",
        h_type="HNAP",
        description="Ar-NH2 + H+ <-> Ar-NH3+  (anilinium, pKa 4.60)",
    ),
    "pyridinic": ProtonationState(
        pKa=5.20,
        kind="basic",
        neutral_type="NPY",
        ionized_type="NPYP",
        h_type="HPYP",
        description="pyridinic N + H+ <-> pyridinium NH+  (pyridinium, pKa 5.20)",
    ),
}

# Group kinds, for callers that need to branch on the sign of the transition
# without reaching into the table.
ACIDIC_GROUPS = frozenset(
    g for g, s in PROTONATION_STATES.items() if s.kind == "acidic"
)
BASIC_GROUPS = frozenset(
    g for g, s in PROTONATION_STATES.items() if s.kind == "basic"
)

# Physically meaningful pH bounds.  Outside this range the Henderson-Hasselbalch
# fraction saturates anyway, and the request is much more likely a unit error.
PH_MIN = 0.0
PH_MAX = 14.0


# Common PAH structures (SMILES notation)
# All entries validated: correct carbon count, 100% aromatic, all atoms in 6-membered rings.
# Hex-lattice entries are programmatically generated compact graphene-nanoflake topologies.
PAH_LIBRARY = {
    # --- 6 carbons ---
    "benzene": {
        "smiles": "c1ccccc1",
        "num_atoms": 6,
        "num_aromatic": 6,
        "molecular_formula": "C6H6",
        "references": "Basic aromatic ring",
    },
    # --- 10 carbons ---
    "naphthalene": {
        "smiles": "c1ccc2ccccc2c1",
        "num_atoms": 10,
        "num_aromatic": 10,
        "molecular_formula": "C10H8",
        "references": "Two fused rings",
    },
    # --- 14 carbons ---
    "anthracene": {
        "smiles": "c1ccc2cc3ccccc3cc2c1",
        "num_atoms": 14,
        "num_aromatic": 14,
        "molecular_formula": "C14H10",
        "references": "Three fused rings (linear acene)",
    },
    "phenanthrene": {
        "smiles": "c1ccc2ccc3ccccc3c2c1",
        "num_atoms": 14,
        "num_aromatic": 14,
        "molecular_formula": "C14H10",
        "references": "Three fused rings (angular)",
    },
    # --- 16 carbons ---
    "pyrene": {
        "smiles": "c1cc2ccc3cccc4ccc(c1)c2c34",
        "num_atoms": 16,
        "num_aromatic": 16,
        "molecular_formula": "C16H10",
        "references": "Four fused rings (pericondensed)",
    },
    # --- 18 carbons ---
    "chrysene": {
        "smiles": "c1ccc2c(c1)cc1ccc3ccccc3c1c2",
        "num_atoms": 18,
        "num_aromatic": 18,
        "molecular_formula": "C18H12",
        "references": "Four fused rings (chrysene topology)",
    },
    "tetracene": {
        "smiles": "c1ccc2cc3cc4ccccc4cc3cc2c1",
        "num_atoms": 18,
        "num_aromatic": 18,
        "molecular_formula": "C18H12",
        "references": "Four fused rings (linear acene, naphthacene)",
    },
    "triphenylene": {
        "smiles": "c1ccc2c(c1)c1ccccc1c1ccccc21",
        "num_atoms": 18,
        "num_aromatic": 18,
        "molecular_formula": "C18H12",
        "references": "Four fused rings (triphenylene topology)",
    },
    # --- 22 carbons ---
    "pentacene": {
        "smiles": "c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1",
        "num_atoms": 22,
        "num_aromatic": 22,
        "molecular_formula": "C22H14",
        "references": "Five fused rings (linear acene)",
    },
    "picene": {
        "smiles": "c1ccc2cc3ccc4ccc5ccccc5c4c3cc2c1",
        "num_atoms": 22,
        "num_aromatic": 22,
        "molecular_formula": "C22H14",
        "references": "Five fused rings (picene/[5]helicene topology)",
    },
    "hex_lattice_22": {
        "smiles": "c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56",
        "num_atoms": 22,
        "num_aromatic": 22,
        "molecular_formula": "C22H12",
        "references": "Six-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 24 carbons ---
    "coronene": {
        "smiles": "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",
        "num_atoms": 24,
        "num_aromatic": 24,
        "molecular_formula": "C24H12",
        "references": "Seven fused rings (central + 6 surrounding)",
    },
    # --- 26 carbons ---
    "hexacene": {
        "smiles": "c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1",
        "num_atoms": 26,
        "num_aromatic": 26,
        "molecular_formula": "C26H16",
        "references": "Six fused rings (linear acene)",
    },
    "dibenzo_bc_ef_coronene": {
        "smiles": "c1ccc2cc3cc4ccc5ccc6ccccc6c5c4cc3cc2c1",
        "num_atoms": 26,
        "num_aromatic": 26,
        "molecular_formula": "C26H14",
        "references": "Six fused rings (angular pericondensed topology)",
    },
    # --- 28 carbons ---
    "hex_lattice_28": {
        "smiles": "c1ccc2c(c1)c1ccc3ccc4ccc5ccc6ccc2c2c6c5c4c3c12",
        "num_atoms": 28,
        "num_aromatic": 28,
        "molecular_formula": "C28H14",
        "references": "Eight-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 30 carbons ---
    "hex_lattice_30": {
        "smiles": "c1cc2ccc3cc4ccc5ccc6ccc7ccc8c(c1)c2c3c1c4c5c6c7c81",
        "num_atoms": 30,
        "num_aromatic": 30,
        "molecular_formula": "C30H14",
        "references": "Nine-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 38 carbons ---
    "hex_lattice_38": {
        "smiles": "c1cc2cc3ccc4cc5cccc6c7ccc8ccc9ccc%10c(c1)c2c1c3c4c(c56)c2c7c8c9c%10c12",
        "num_atoms": 38,
        "num_aromatic": 38,
        "molecular_formula": "C38H16",
        "references": "Twelve-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 40 carbons ---
    "hex_lattice_40": {
        "smiles": "c1cc2cc3ccc4cc5ccc6ccc7cc8ccc9ccc%10c(c1)c2c1c3c4c2c5c6c7c3c8c9c%10c1c32",
        "num_atoms": 40,
        "num_aromatic": 40,
        "molecular_formula": "C40H16",
        "references": "Thirteen-ring compact graphene nanoflake (hex-lattice seed)",
    },
}

# Atomic masses (in amu)
ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.01,
    "N": 14.007,
    "O": 15.999,
    "S": 32.065,
    "P": 30.974,
    "Cl": 35.45,
    "Br": 79.904,
}

# Van der Waals radii (in Angstrom) for steric clash detection
VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
}

# Van der Waals diameter of graphitic carbon (graphite interlayer spacing).
# Used as effective sheet thickness when computing slit-pore geometry.
CARBON_VDW_DIAMETER = 3.4  # Angstroms

# --- Hydrogen-bond-aware clash detection ----------------------------------
# A polar hydrogen sitting close to an O/N acceptor is a hydrogen bond, not a
# steric clash.  The generic 0.75 × vdW-sum floor is 2.04 Å for an O/H pair,
# which lands squarely inside the physical H...A range (~1.6-2.2 Å), so every
# intramolecular H-bond between adjacent -OH groups is otherwise reported as a
# clash.  On high-oxygen chars (O/C >= ~0.2, e.g. 400 °C softwood) such pairs
# are unavoidable and strict validation fails on every seed.
#
# Donor-acceptor pairs that satisfy the angle criterion below are held to this
# reduced floor instead.  It still catches genuine overlap: an H...O contact
# shorter than this is too close even for a low-barrier hydrogen bond.
HBOND_MIN_H_ACCEPTOR_DISTANCE = 1.5  # Angstroms

# Minimum D-H...A angle (degrees) for a contact to count as a hydrogen bond.
# Real H-bonds are near-linear (typically > 120°); 90° is a deliberately
# permissive gate that only requires the H to point *toward* the acceptor
# rather than the acceptor being jammed into the side of the D-H bond.
HBOND_MIN_DHA_ANGLE_DEG = 90.0

# Elements that act as hydrogen-bond donors (when carrying an H) and acceptors.
HBOND_DONOR_ACCEPTOR_ELEMENTS = (7, 8)  # N, O (atomic numbers)

# Slack (Å) below a clash floor before a contact is reported as a clash.  A
# contact is a clash only when it is deeper than this below the floor
# (0.75 × vdW-sum, or HBOND_MIN_H_ACCEPTOR_DISTANCE for a hydrogen bond).
#
# Neither floor is meaningful to 0.01 Å: the ETKDG embedding and MMFF relaxation
# leave contacts scattered by more than that, so a hard `distance < floor`
# knife-edge reports contacts sitting a hundredth of an ångström inside the
# floor as clashes.  On the densest chars (e.g. 400 °C manure, O/C 0.24) an
# H-bond 0.01 Å shorter than 1.5 Å, or a peri H...C 0.005 Å inside the vdW
# floor, then fails strict validation on seed after seed even though GROMACS
# energy minimisation would move the atom in the first step.
#
# 0.05 Å is well below a genuine overlap: the real clashes this system catches
# (two functional groups embedded on top of each other) sit 0.2-0.5 Å inside
# the floor.  Anything within 0.05 Å of the floor is embedding noise, not a
# clash.
CLASH_SEVERITY_TOLERANCE = 0.05  # Angstroms

# --- Hex-lattice crowding relief -------------------------------------------
# The hex-lattice path places ring atoms on an exact flat lattice and every
# substituent radially outward from its parent.  That is right almost
# everywhere and wrong at helicene-type edges: at a cove (two edge carbons
# 2.46 Å apart, four bonds apart -- a [4]helicene motif) both substituents point
# into the same vacant lattice site and land 0.1-0.6 Å apart, and the cove
# carbons themselves sit at 2.46 Å, inside the 2.50 Å C...C floor.  Real
# molecules relieve this by twisting out of plane; GROMACS energy minimisation
# of an exported sheet opens such coves to ~2.95 Å and buckles the ring carbons
# up to ~1.0 Å out of plane.  The relief pass reproduces that locally.
#
# X-O-H angle (degrees) at which a hydroxyl hydrogen is placed.  Radial
# placement put it collinear with the C-O bond (180°) on ~30% of hydroxyls --
# the ones no clash drew attention to -- which is both wrong geometry and a
# singular torsion for any force field handed the structure.
HYDROXYL_H_PLACEMENT_ANGLE_DEG = 108.5

# Largest out-of-plane tilt (degrees) applied to a crowded substituent before
# relaxation, and the step the search takes toward it.  75° stops short of
# standing a substituent perpendicular to the sheet, where its tilt axis
# degenerates.
LATTICE_TILT_MAX_DEG = 75.0
LATTICE_TILT_STEP_DEG = 7.5

# Restrained relaxation.  Each ring atom is held by a flat-bottomed tether: free
# within LATTICE_TETHER_RADIUS of its lattice position, penalised beyond it, so
# a cove can open and buckle locally while the sheet as a whole stays where the
# lattice put it.  Aromatic bonds are held to LATTICE_BOND_RANGE -- the
# force-field copy types them as sp2 single bonds, whose UFF rest length
# (~1.47 Å) would otherwise stretch every bond in the sheet.  The bond term is
# stiff, so the constraint must be stiffer: at 500 the bonds settled at ~1.46 Å,
# and a relaxed sheet sits at the constraint's upper edge, which is why the
# range tops out at the lattice value rather than straddling it (bonds land at
# 1.41-1.43 Å, median 1.42).  Much stiffer (50000) left the minimiser too
# little room and relieved fewer structures.  Values were chosen on a
# 36-structure grid (100-200 C, H/C 0.35, O/C 0-0.2): a 0.2 Å tether relieved
# fewer structures, a 0.5 Å one let the whole sheet drift.
LATTICE_TETHER_RADIUS = 0.3  # Angstroms
LATTICE_TETHER_FORCE_CONSTANT = 10.0
LATTICE_BOND_RANGE = (1.40, 1.42)  # Angstroms
LATTICE_BOND_FORCE_CONSTANT = 20000.0
# The minimiser stops early once converged, so this is a ceiling, not a cost.
# At 400 it was routinely cut off mid-descent: 150 C, H/C 0.35, O/C 0.3, seed 4
# kept a fjord contact at 400 and cleared it at 1000 from every start tried.
LATTICE_RELAX_MAX_ITERATIONS = 2000

# Largest RMS deviation (Å) of the ring atoms from their best-fit plane that a
# relaxed sheet may show and still count as flat.  Local cove buckling of up to
# ~1 Å at a few sites puts a relieved 100-200 C sheet at 0.25-0.55 Å RMS
# (GROMACS EM of the same sheets: ~0.33 Å); a sheet the relaxation folded is
# well beyond this.  A result over it is discarded in favour of the unrelaxed
# (tilted-only) coordinates.
LATTICE_MAX_SHEET_RMS = 0.75  # Angstroms

# --- Bond-length validation ------------------------------------------------
# COVALENT_RADII are *single-bond* radii, so their sum only predicts a single
# bond.  Scale by bond order to get the expected length: an aromatic C-C is
# 1.40 Å, not the 1.52 Å the radii sum implies, and a C=O is 1.23 Å, not 1.42 Å.
# Reporting the unscaled sum made every such message wrong about what it
# expected, even when the bond really was out of range.
#
# Factors are the ratio of the observed length to the single-bond radii sum:
# aromatic C-C 1.40/1.52 = 0.92; C=O 1.23/1.42 = 0.87; C#C 1.20/1.52 = 0.79.
BOND_ORDER_LENGTH_FACTORS = {
    "SINGLE": 1.00,
    "AROMATIC": 0.92,
    "DOUBLE": 0.87,
    "TRIPLE": 0.79,
}

# Fractional tolerance on the expected bond length.  These are tighter than the
# old 0.8/1.5 band on purpose: correcting `expected` downward for aromatic and
# multiple bonds would otherwise *lower* the absolute floor and let a genuinely
# compressed bond through.  At these factors an aromatic C-C is accepted over
# 1.19-1.96 Å, versus 1.22-2.28 Å before -- comparable at the low end, and no
# longer absurdly permissive at the high end.
BOND_LENGTH_MIN_FACTOR = 0.85
BOND_LENGTH_MAX_FACTOR = 1.40

# ---------------------------------------------------------------------------
# Experimental-data model provenance & tunables
# ---------------------------------------------------------------------------
# Primary characterization dataset behind the temperature/feedstock composition
# model (see :mod:`biochar.models.temperature_model`):
UC_DAVIS_DB_URL = "https://biochar.ucdavis.edu/"   # UC Davis Biochar Database
# Methodological parent: Wood, Mašek & Erastova, Cell Reports Physical Science
# 5(7), 2024, DOI 10.1016/j.xcrp.2024.102036.
#
# Minimum aromaticity (%) the PAH skeleton builder can faithfully realise.  When
# the data-derived aromaticity (from temperature/feedstock) falls below this it
# is clamped to this floor and a warning is emitted.
MIN_BUILDABLE_AROMATICITY = 70.0

# Ring-curvature ratios of non-graphitizing carbons (Wood et al. 2024): the
# published island building blocks average roughly 10 hexagons : 2 pentagons :
# 1 heptagon.  As per-ring-addition probabilities over the 13-ring total this is
# 2/13 pentagons and 1/13 heptagons (hexagons are the ~10/13 remainder).  Pass
# these as ``defect_fraction`` / ``heptagon_fraction`` to reproduce that mix.
WOOD_PENTAGON_FRACTION = 2.0 / 13.0   # ≈ 0.154
WOOD_HEPTAGON_FRACTION = 1.0 / 13.0   # ≈ 0.077

# Covalent radii (in Angstrom) for bond length validation
COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "S": 1.05,
    "P": 1.07,
    "Cl": 1.02,
    "Br": 1.20,
}


def get_atom_mass(atom_type: str) -> float:
    """Get atomic mass from OPLS atom type."""
    if atom_type in OPLS_ATOM_TYPES:
        return OPLS_ATOM_TYPES[atom_type][1]
    # Extract element from atom type (first character or two)
    element = atom_type[0] if atom_type[0].isupper() else atom_type[:2]
    return ATOMIC_MASSES.get(element, 12.01)


def get_vdw_radius(element: str) -> float:
    """Get Van der Waals radius for element."""
    return VDW_RADII.get(element, 1.70)


def get_covalent_radius(element: str) -> float:
    """Get covalent radius for element."""
    return COVALENT_RADII.get(element, 0.76)
