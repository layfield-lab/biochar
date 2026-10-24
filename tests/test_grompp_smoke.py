"""
Depth 4 of rqm/opls-typing.md: run grompp on an exported structure and read
back the binary topology it produces.

Depths 1-3 read oplsaa.ff the way grompp would and check that every term the
topology emits can be resolved. That cannot catch a term the topology never
emits, or a preprocessor macro that fails to expand -- a missing interaction
shows up only as a count of zero in the .tpr. Both of those were real: exported
topologies carried no [ pairs ] section and no improper terms, and every depth-3
check passed throughout.

These need the `gmx` binary, not just the forcefield files. CI installs the full
gromacs package and asserts the binary is present, so a skip here is a local
environment gap rather than the normal state.
"""

import re
import shutil
import subprocess

import pytest

from biochar.export.gromacs_export import GromacsExporter
from biochar.pipeline.biochar_generator import BiocharGenerator, GeneratorConfig

GMX = shutil.which("gmx") or shutil.which("gmx_mpi")
requires_gmx = pytest.mark.skipif(
    GMX is None,
    reason="no gmx binary on PATH; the conda env ships one -- "
           "export PATH=\"$CONDA_PREFIX/bin:$PATH\"",
)

MINIMAL_MDP = """integrator    = steep
nsteps        = 10
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
pbc           = xyz
"""


def _export(tmp_path, basename, **cfg):
    kwargs = {"target_num_carbons": 20, "strict": False, "seed": 1}
    kwargs.update(cfg)
    gen = BiocharGenerator(GeneratorConfig(**kwargs))
    mol, coords, _ = gen.generate()
    GromacsExporter(str(tmp_path)).export(
        mol, coords, gen.atom_types, gen.charges,
        molecule_name="BCX", basename=basename, include_periodic_box=True,
    )
    (tmp_path / "min.mdp").write_text(MINIMAL_MDP)
    return mol


def _grompp(tmp_path, basename, *extra):
    return subprocess.run(
        [GMX, "grompp", "-f", str(tmp_path / "min.mdp"),
         "-c", str(tmp_path / f"{basename}.gro"),
         "-p", str(tmp_path / f"{basename}.top"),
         "-o", str(tmp_path / f"{basename}.tpr"), *extra],
        cwd=tmp_path, capture_output=True, text=True,
    )


def _warnings(stderr: str) -> list:
    """The numbered warning blocks grompp emitted.

    Counting the bare word overcounts: grompp repeats it in the "Too many
    warnings (1)" banner when it refuses, so a single warning reads as two.
    Each real warning opens with `WARNING <n> [file ...]`.
    """
    blocks = re.split(r"^WARNING \d+ \[", stderr, flags=re.MULTILINE)[1:]
    return [" ".join(b.split())[:200] for b in blocks]


def _interaction_counts(tmp_path, basename):
    """Entry counts per interaction class, read back from the .tpr.

    `gmx dump` prints `nr:` as the number of *fields*, not interactions: each
    entry is the parameter index plus its atoms, so a bond is 3 and a dihedral
    is 5. The ratios matter here, not the absolute numbers.
    """
    dump = subprocess.run(
        [GMX, "dump", "-s", str(tmp_path / f"{basename}.tpr")],
        cwd=tmp_path, capture_output=True, text=True,
    )
    counts, current = {}, None
    for raw in dump.stdout.splitlines():
        line = raw.strip()
        if line.endswith(":") and not line.startswith("nr"):
            current = line[:-1]
        elif line.startswith("nr:") and current:
            counts.setdefault(current, int(line.split(":")[1]))
            current = None
    return counts


@requires_gmx
class TestGromppAcceptsTheTopology:
    # rq-be5907a0
    def test_a_neutral_topology_is_accepted_without_warnings(self, tmp_path):
        _export(tmp_path, "neutral")
        result = _grompp(tmp_path, "neutral")

        assert result.returncode == 0, (
            f"grompp refused the exported topology:\n{result.stderr[-2000:]}"
        )
        assert not _warnings(result.stderr), (
            f"grompp warned on a neutral structure: {_warnings(result.stderr)}\n"
            f"{result.stderr[-2000:]}"
        )
        assert (tmp_path / "neutral.tpr").exists()

    # rq-85bcfa59
    def test_every_emitted_term_reaches_the_binary_topology(self, tmp_path):
        _export(tmp_path, "neutral")
        assert _grompp(tmp_path, "neutral").returncode == 0

        counts = _interaction_counts(tmp_path, "neutral")

        # The 1-4 pair list. Zero here is what a missing [ pairs ] section looks
        # like, and it is invisible to every depth-3 check.
        assert counts.get("LJ-14", 0) > 0, (
            f"no 1-4 pairs in the binary topology; nrexcl=3 excluded them and "
            f"nothing restored them. counts={counts}"
        )
        # Proper torsions, written as Ryckaert-Bellemans (funct 3).
        assert counts.get("Ryckaert-Bell.", 0) > 0, f"no proper torsions: {counts}"
        # Impropers. OPLS implements them in the periodic *proper* form (funct 1),
        # so they land under "Proper Dih." and "Improper Dih." is legitimately 0 --
        # checking the wrong one of those two would pass on an empty topology.
        assert counts.get("Proper Dih.", 0) > 0, (
            f"no improper terms; aromatic centres are free to pyramidalise. "
            f"counts={counts}"
        )

    # rq-256e31b8
    def test_a_charged_structure_needs_exactly_one_warning_allowed(self, tmp_path):
        mol = _export(tmp_path, "charged", O_C_ratio=0.2, pH=11.0)
        assert sum(a.GetFormalCharge() for a in mol.GetAtoms()) != 0, (
            "fixture is void: raise the pH until the structure carries a charge"
        )

        refused = _grompp(tmp_path, "charged")
        assert refused.returncode != 0, "grompp accepted a net-charged system under PME"

        warnings = _warnings(refused.stderr)
        assert len(warnings) == 1, (
            f"expected exactly one warning to account for, got {warnings}\n"
            f"{refused.stderr[-2000:]}"
        )
        assert "net charge" in warnings[0], (
            f"the one warning is not the net-charge one: {warnings[0]}"
        )

        # This is the budget md_setup spends on pre-neutralisation stages.
        allowed = _grompp(tmp_path, "charged", "-maxwarn", "1")
        assert allowed.returncode == 0, (
            f"-maxwarn 1 was not enough:\n{allowed.stderr[-2000:]}"
        )


# rq-e72878a2
@requires_gmx
@pytest.mark.parametrize("linker", ["biaryl", "diaryl_ether", "aryl_alkyl_ether", "propyl"])
def test_grompp_accepts_each_linked_unit_linker(tmp_path, linker):
    """Each linker of a linked multi-unit structure (rqm/linked-units.md) emits
    bonds, angles and torsions stock oplsaa.ff resolves: CA-OS-CA, CT-OS-CA,
    CA-CT-CT and the biaryl CA-CA across units."""
    _export(tmp_path, linker, target_num_carbons=40, linked_units=True, num_units=3,
            linker_mix={linker: 1.0}, H_C_ratio=0.6, O_C_ratio=0.05)
    result = _grompp(tmp_path, linker)
    assert result.returncode == 0, (
        f"grompp refused a {linker}-linked topology:\n{result.stderr[-2000:]}"
    )
    assert not _warnings(result.stderr), _warnings(result.stderr)


def test_the_grompp_skip_is_visible():
    """A depth that quietly did not run is worse than one never claimed.

    Mirrors TestForcefieldAbsentSkipIsVisible in test_opls_type_map.py: assert
    the gate is a real skip with a reason naming what is missing, rather than a
    silent pass.
    """
    if GMX is not None:
        pytest.skip("gmx is present; the absent-binary path is what needs proving")

    marker = requires_gmx.markname
    assert marker == "skipif"
    reason = requires_gmx.kwargs["reason"]
    assert "gmx" in reason and "PATH" in reason, (
        f"the skip reason does not say what is missing or how to fix it: {reason}"
    )
