"""
Validation Engine

Comprehensive validation for biochar structures including valence checking.
"""

from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass

from rdkit import Chem
import numpy as np

from .heteroatom_assignment import CompositionInfo
from .geometry_3d import GeometryValidator
from .valence import ValenceValidator


@dataclass
class ValidationReport:
    """Result of validation check."""

    is_valid: bool
    errors: List[str]
    warnings: List[str]
    metrics: Dict[str, float]


class CompositionValidator:
    """Validate molecular composition against targets."""

    @staticmethod
    def validate(
        composition: CompositionInfo,
        target_H_C: float,
        target_O_C: float,
        target_aromaticity: Optional[float] = None,
        H_C_tolerance: float = 0.10,
        O_C_tolerance: float = 0.10,
        aromaticity_tolerance: float = 5.0,
    ) -> ValidationReport:
        """
        Validate composition against targets.

        Args:
            composition: CompositionInfo object
            target_H_C: Target H/C ratio
            target_O_C: Target O/C ratio
            target_aromaticity: Target aromaticity % (optional)
            H_C_tolerance: Tolerance for H/C ratio
            O_C_tolerance: Tolerance for O/C ratio
            aromaticity_tolerance: Tolerance for aromaticity %

        Returns:
            ValidationReport
        """
        errors = []
        warnings = []
        metrics = {
            "actual_H_C": composition.H_C_ratio,
            "target_H_C": target_H_C,
            "H_C_error_pct": abs(composition.H_C_ratio - target_H_C) / max(target_H_C, 0.01) * 100,
            "actual_O_C": composition.O_C_ratio,
            "target_O_C": target_O_C,
            "O_C_error_pct": abs(composition.O_C_ratio - target_O_C) / max(target_O_C, 0.01) * 100,
        }

        # Check H/C ratio
        H_C_error_pct = metrics["H_C_error_pct"]
        if H_C_error_pct > H_C_tolerance * 100:
            # The "H/C ratio" prefix is matched by
            # BiocharGenerator._h_c_floor_reason to add the unreachable-floor hint.
            errors.append(
                f"H/C ratio {composition.H_C_ratio:.3f} exceeds tolerance "
                f"(target: {target_H_C:.3f}, error: {H_C_error_pct:.1f}%)"
            )
        elif H_C_error_pct > H_C_tolerance * 100 * 0.5:
            warnings.append(
                f"H/C ratio {composition.H_C_ratio:.3f} is at edge of tolerance "
                f"(target: {target_H_C:.3f}, error: {H_C_error_pct:.1f}%)"
            )

        # Check O/C ratio
        if target_O_C > 0:
            O_C_error_pct = metrics["O_C_error_pct"]
            if O_C_error_pct > O_C_tolerance * 100:
                errors.append(
                    f"O/C ratio {composition.O_C_ratio:.3f} exceeds tolerance "
                    f"(target: {target_O_C:.3f}, error: {O_C_error_pct:.1f}%)"
                )
            elif O_C_error_pct > O_C_tolerance * 100 * 0.5:
                warnings.append(
                    f"O/C ratio {composition.O_C_ratio:.3f} is at edge of tolerance "
                    f"(target: {target_O_C:.3f}, error: {O_C_error_pct:.1f}%)"
                )
        else:
            if composition.O_C_ratio > 0.01:
                warnings.append(
                    f"Target O/C is zero but molecule has oxygens (O/C={composition.O_C_ratio:.3f})"
                )

        # Check atom counts are reasonable
        if composition.num_carbons < 4:
            errors.append("Too few carbons for biochar structure")
        if composition.num_hydrogens == 0:
            errors.append("No hydrogens in structure")

        # Check functional groups if specified
        if composition.functional_groups:
            metrics["num_functional_groups"] = sum(composition.functional_groups.values())

        return ValidationReport(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            metrics=metrics,
        )


class ChemicalFeasibilityValidator:
    """Validate chemical feasibility of structures."""

    @staticmethod
    def validate(mol: Chem.Mol) -> ValidationReport:
        """
        Validate chemical feasibility.

        Args:
            mol: RDKit molecule

        Returns:
            ValidationReport
        """
        errors = []
        warnings = []
        metrics = {}

        # Check for valid RDKit molecule
        if mol is None:
            errors.append("Molecule object is None")
            return ValidationReport(is_valid=False, errors=errors, warnings=warnings, metrics=metrics)

        if mol.GetNumAtoms() == 0:
            errors.append("Molecule has no atoms")
            return ValidationReport(is_valid=False, errors=errors, warnings=warnings, metrics=metrics)

        # Check valence constraints
        valence_errors = ChemicalFeasibilityValidator._check_valences(mol)
        errors.extend(valence_errors)

        # Check for unusual atoms
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num > 18 and atomic_num not in [16, 17, 35]:  # Allow S, Cl, Br
                warnings.append(f"Unusual atom: {atom.GetSymbol()} (atomic number {atomic_num})")

        # Check aromaticity
        num_aromatic = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        metrics["num_aromatic_atoms"] = num_aromatic
        metrics["aromaticity_percent"] = num_aromatic / mol.GetNumAtoms() * 100

        # Check connectivity, on a copy.
        #
        # Sanitisation is being used here as a question -- does RDKit accept
        # this molecule -- and its only output is the warning below. But it is
        # not a read: it kekulises, and kekulisation of a large fused-ring
        # sheet fails, having already rewritten aromatic bonds to SINGLE on its
        # way to raising. Asking the question of the caller's own molecule
        # therefore answered it and destroyed the subject in the same call,
        # leaving a sheet whose atoms are still flagged aromatic and whose
        # bonds are all single. Nothing downstream re-perceives it, so that is
        # what was written to the topology.
        try:
            Chem.SanitizeMol(Chem.Mol(mol))
        except Exception as e:
            warnings.append(f"Molecule sanitization raised warning: {str(e)}")

        # Check for charged atoms
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        metrics["total_formal_charge"] = total_charge

        if abs(total_charge) > 2:
            warnings.append(f"High total formal charge: {total_charge}")

        return ValidationReport(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            metrics=metrics,
        )

    @staticmethod
    def _check_valences(mol: Chem.Mol) -> List[str]:
        """Check for valence constraint violations using comprehensive validator."""
        # Use the dedicated valence validator
        is_valid, errors = ValenceValidator.validate_molecule(mol)
        return errors


class StructureValidator:
    """Validate overall structure quality."""

    @staticmethod
    def validate(
        mol: Chem.Mol,
        coords: Optional[np.ndarray] = None,
    ) -> ValidationReport:
        """
        Validate overall structure quality.

        Args:
            mol: RDKit molecule
            coords: Optional 3D coordinates

        Returns:
            ValidationReport
        """
        errors = []
        warnings = []
        metrics = {}

        if mol is None:
            errors.append("Molecule is None")
            return ValidationReport(is_valid=False, errors=errors, warnings=warnings, metrics=metrics)

        # Check molecule size
        num_atoms = mol.GetNumAtoms()
        metrics["num_atoms"] = num_atoms

        if num_atoms < 6:
            warnings.append("Molecule is very small (< 6 atoms)")
        if num_atoms > 5000:
            warnings.append("Molecule is very large (> 5000 atoms)")

        # Check atom composition
        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

        metrics.update({f"num_{symbol}": count for symbol, count in atom_counts.items()})

        if atom_counts.get("C", 0) == 0:
            errors.append("No carbon atoms in structure")
        if atom_counts.get("H", 0) == 0:
            errors.append("No hydrogen atoms in structure")

        # Check connectivity
        if not StructureValidator._is_connected(mol):
            errors.append("Molecule is not fully connected")

        # Check 3D geometry if coordinates provided
        if coords is not None:
            geom_valid, geom_errors = GeometryValidator.validate_geometry(mol, coords)
            if not geom_valid:
                # Every error, not a sample. Truncating here biases the report
                # toward whichever check runs first: a sheet with one out-of-range
                # bond and a long tail of clash contacts shows the contacts and
                # drops the bond, which is how the bond-length defect stayed hidden
                # (docs/solutions/bugs/physical-features-misread-as-geometry-errors.md).
                # Callers that want a summary truncate for display -- print_summary
                # already shows three and states the full count. A caller handed a
                # truncated list cannot recover what was dropped.
                errors.extend(geom_errors)

            # Check planarity of aromatic systems
            planarity, assessment = GeometryValidator.measure_ring_planarity(mol, coords)
            metrics["avg_ring_planarity_deviation_A"] = planarity
            if planarity > 0.5:
                warnings.append(f"Poor aromatic ring planarity: {assessment}")

        return ValidationReport(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            metrics=metrics,
        )

    @staticmethod
    def _is_connected(mol: Chem.Mol) -> bool:
        """Check if molecule is fully connected."""
        if mol.GetNumAtoms() == 0:
            return True

        visited = set()
        stack = [0]  # Start from atom 0

        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)

            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor.GetIdx())

        return len(visited) == mol.GetNumAtoms()


class ValidationEngine:
    """Master validation engine combining all validators."""

    @staticmethod
    def validate_complete(
        mol: Chem.Mol,
        composition: CompositionInfo,
        coords: Optional[np.ndarray] = None,
        target_H_C: float = 0.5,
        target_O_C: float = 0.1,
        target_aromaticity: Optional[float] = None,
        H_C_tolerance: float = 0.10,
        O_C_tolerance: float = 0.10,
    ) -> Tuple[bool, List[str], List[str], Dict]:
        """
        Run comprehensive validation on structure.

        Args:
            mol: RDKit molecule
            composition: CompositionInfo object
            coords: Optional 3D coordinates
            target_H_C: Target H/C ratio
            target_O_C: Target O/C ratio
            target_aromaticity: Target aromaticity %
            H_C_tolerance: H/C tolerance
            O_C_tolerance: O/C tolerance

        Returns:
            (is_valid, all_errors, all_warnings, metrics_dict)
        """
        all_errors = []
        all_warnings = []
        all_metrics = {}

        # Composition validation
        comp_report = CompositionValidator.validate(
            composition,
            target_H_C,
            target_O_C,
            target_aromaticity,
            H_C_tolerance,
            O_C_tolerance,
        )
        all_errors.extend(comp_report.errors)
        all_warnings.extend(comp_report.warnings)
        all_metrics.update(comp_report.metrics)

        # Chemical feasibility validation
        chem_report = ChemicalFeasibilityValidator.validate(mol)
        all_errors.extend(chem_report.errors)
        all_warnings.extend(chem_report.warnings)
        all_metrics.update(chem_report.metrics)

        # Structure validation
        struct_report = StructureValidator.validate(mol, coords)
        all_errors.extend(struct_report.errors)
        all_warnings.extend(struct_report.warnings)
        all_metrics.update(struct_report.metrics)

        return len(all_errors) == 0, all_errors, all_warnings, all_metrics
