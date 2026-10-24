"""
Tests for the data-driven temperature × feedstock property model.
"""

import pytest

from biochar.models.temperature_model import (
    TemperatureModel,
    properties,
    get_default_model,
    compare_models,
    VALID_FEEDSTOCKS,
    _classify_feedstock,
    _MODEL_PATH,
)
from biochar.pipeline.biochar_generator import GeneratorConfig, BiocharGenerator
from biochar.constants import MIN_BUILDABLE_AROMATICITY, UC_DAVIS_DB_URL

M = get_default_model()


# ---------------------------------------------------------------------------
# Pooled predictions
# ---------------------------------------------------------------------------

class TestPooledPredictions:
    def test_feedstocks_constant(self):
        assert VALID_FEEDSTOCKS == (
            "corn_stover", "grass", "hardwood", "manure", "softwood", "wood"
        )

    def test_hc_decreases_with_temperature(self):
        assert (M.predict(300, "H_C_ratio")
                > M.predict(600, "H_C_ratio")
                > M.predict(800, "H_C_ratio"))

    def test_oc_decreases_with_temperature(self):
        assert M.predict(300, "O_C_ratio") > M.predict(800, "O_C_ratio")

    def test_aromaticity_increases_with_temperature(self):
        assert (M.composition(300)["aromaticity_percent"]
                < M.composition(800)["aromaticity_percent"])

    def test_hc_plausible_at_600(self):
        assert 0.2 <= M.predict(600, "H_C_ratio") <= 0.6

    def test_oc_plausible_at_600(self):
        assert 0.03 <= M.predict(600, "O_C_ratio") <= 0.2

    def test_ratios_nonnegative(self):
        for T in (100, 300, 600, 1000):
            assert M.predict(T, "H_C_ratio") >= 0
            assert M.predict(T, "O_C_ratio") >= 0

    def test_clamp_below_grid(self):
        assert M.predict(50, "H_C_ratio") == M.predict(100, "H_C_ratio")

    def test_clamp_above_grid(self):
        assert M.predict(1500, "H_C_ratio") == M.predict(1000, "H_C_ratio")

    def test_unknown_property_raises(self):
        with pytest.raises(KeyError):
            M.predict(600, "not_a_property")

    def test_pH_increases_with_temperature(self):
        assert M.predict(800, "pH") > M.predict(300, "pH")


# ---------------------------------------------------------------------------
# Feedstock behaviour
# ---------------------------------------------------------------------------

class TestFeedstock:
    def test_overrides_differ_between_groups(self):
        assert (M.predict(450, "H_C_ratio", "hardwood")
                != M.predict(450, "H_C_ratio", "grass"))

    def test_feedstock_differs_from_pooled(self):
        assert (M.predict(450, "H_C_ratio", "grass")
                != M.predict(450, "H_C_ratio"))

    def test_wood_falls_back_to_pooled(self):
        # wood is its own category but below the gate -> no override -> pooled.
        assert (M.predict(450, "H_C_ratio", "wood")
                == M.predict(450, "H_C_ratio"))

    def test_out_of_range_temperature_falls_back_to_pooled(self):
        # corn_stover support starts ~300 C; at 100 C it must use the pooled curve.
        assert (M.predict(100, "H_C_ratio", "corn_stover")
                == M.predict(100, "H_C_ratio"))

    def test_classify_maps_known_labels(self):
        assert _classify_feedstock("Soft Wood") == "softwood"
        assert _classify_feedstock("Hard Wood") == "hardwood"
        assert _classify_feedstock("Grass") == "grass"
        assert _classify_feedstock("Corn stover") == "corn_stover"
        assert _classify_feedstock("Manure") == "manure"

    def test_classify_wood_is_own_category(self):
        assert _classify_feedstock("Wood") == "wood"
        assert _classify_feedstock("wood") == "wood"
        # must NOT collapse soft/hard into plain wood
        assert _classify_feedstock("Soft Wood") == "softwood"

    def test_classify_unmapped_returns_none(self):
        assert _classify_feedstock("Sludge") is None
        assert _classify_feedstock("Algae") is None
        assert _classify_feedstock("Other/ Mixed") is None


# ---------------------------------------------------------------------------
# properties() query + aromaticity calibration
# ---------------------------------------------------------------------------

class TestPropertiesQuery:
    def test_returns_reference_properties(self):
        p = properties(600)
        for k in ("surface_area_m2g", "pH", "ash_pct", "C_pct", "VM_pct"):
            assert k in p

    def test_includes_composition(self):
        p = properties(600)
        assert "H_C_ratio" in p and "O_C_ratio" in p and "aromaticity_percent" in p

    def test_invalid_feedstock_raises(self):
        with pytest.raises(ValueError):
            properties(600, feedstock="zzz")

    def test_feedstock_query_differs(self):
        assert (properties(450, "hardwood")["H_C_ratio"]
                != properties(450, "grass")["H_C_ratio"])

    def test_aromaticity_from_hc_clamped_high(self):
        # Very high H/C -> negative raw -> clamp to 0. Both of these sit outside
        # the H/C range the relation was fitted over, so the clamp is reported as
        # extrapolation rather than returned as a confident value
        # (rqm/temperature-model.md rq-27f04167).
        with pytest.warns(UserWarning, match="extrapolat"):
            assert M.aromaticity_from_hc(5.0) == 0.0

    def test_aromaticity_from_hc_clamped_low(self):
        # near-zero H/C -> >100 raw -> clamp to 100
        with pytest.warns(UserWarning, match="extrapolat"):
            assert M.aromaticity_from_hc(0.0) == 100.0

    def test_aromaticity_monotonic_in_hc(self):
        assert M.aromaticity_from_hc(0.3) > M.aromaticity_from_hc(1.0)


class TestProvenance:
    def test_provenance_url(self):
        assert M.provenance["primary_url"] == UC_DAVIS_DB_URL == "https://biochar.ucdavis.edu/"


# ---------------------------------------------------------------------------
# GeneratorConfig integration
# ---------------------------------------------------------------------------

class TestGeneratorConfigIntegration:
    def test_bare_defaults_unchanged(self):
        c = GeneratorConfig()
        assert (c.H_C_ratio, c.O_C_ratio, c.aromaticity_percent) == (0.5, 0.1, 90.0)

    def test_temperature_fills_composition(self):
        c = GeneratorConfig(temperature=600)
        assert 0.2 <= c.H_C_ratio <= 0.6
        assert c.O_C_ratio > 0
        assert c.aromaticity_percent >= MIN_BUILDABLE_AROMATICITY

    def test_explicit_value_overrides_temperature(self):
        c = GeneratorConfig(temperature=600, H_C_ratio=0.3)
        assert c.H_C_ratio == 0.3

    def test_feedstock_changes_result(self):
        a = GeneratorConfig(temperature=450, feedstock="hardwood")
        b = GeneratorConfig(temperature=450, feedstock="grass")
        assert a.H_C_ratio != b.H_C_ratio

    def test_low_temperature_clamps_aromaticity(self):
        c = GeneratorConfig(temperature=300)
        assert c.aromaticity_percent == MIN_BUILDABLE_AROMATICITY

    def test_invalid_feedstock_raises(self):
        with pytest.raises(ValueError):
            GeneratorConfig(temperature=600, feedstock="zzz")

    def test_feedstock_without_temperature_is_inert(self):
        c = GeneratorConfig(feedstock="softwood")
        assert c.H_C_ratio == 0.5  # nothing to derive from; historical default

    def test_to_dict_from_dict_roundtrip(self):
        c = GeneratorConfig(temperature=600, feedstock="softwood")
        d = c.to_dict()
        assert d["temperature"] == 600 and d["feedstock"] == "softwood"
        c2 = GeneratorConfig.from_dict(d)
        assert c2.temperature == 600 and c2.feedstock == "softwood"


# ---------------------------------------------------------------------------
# End-to-end + build helpers
# ---------------------------------------------------------------------------

class TestEndToEnd:
    def test_generate_with_temperature(self):
        c = GeneratorConfig(target_num_carbons=24, temperature=600, seed=1, strict=False)
        mol, coords, comp = BiocharGenerator(c).generate()
        assert mol.GetNumAtoms() > 0

    def test_generate_with_feedstock(self):
        c = GeneratorConfig(target_num_carbons=24, temperature=500,
                            feedstock="softwood", seed=1, strict=False)
        mol, coords, comp = BiocharGenerator(c).generate()
        assert mol.GetNumAtoms() > 0


class TestStrictTemperatureMode:
    """Issue #63: strict temperature-mode builds failed on H/C with no reason given.

    Both end-to-end tests above pass strict=False, which is how the failure went
    unnoticed.
    """

    # rq-ef032605
    @pytest.mark.slow
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_a_reachable_target_passes_strict(self, seed):
        # 500 °C softwood asks H/C 0.414, above a compact 250-carbon sheet's
        # ceiling (~0.22) but within elongated growth's reach. It failed on every
        # seed: elongation overshot to 0.512.
        c = GeneratorConfig(target_num_carbons=250, temperature=500,
                            feedstock="softwood", seed=seed)
        mol, coords, comp = BiocharGenerator(c).generate()
        assert abs(comp.H_C_ratio - c.H_C_ratio) <= c.H_C_tolerance * c.H_C_ratio

    # rq-907cd890
    def test_an_unreachable_target_names_the_floor_and_a_size(self):
        from biochar.pipeline.biochar_generator import ValidationError

        # 800 °C softwood asks H/C 0.168; a 100-carbon sheet has ~0.3 of its
        # carbons on the edge, each carrying a hydrogen or a group.
        c = GeneratorConfig(target_num_carbons=100, temperature=800,
                            feedstock="softwood", seed=1)
        with pytest.raises(ValidationError) as exc:
            BiocharGenerator(c).generate()
        message = str(exc.value)
        assert "below the lowest H/C" in message, message
        assert "target_num_carbons=" in message, message


class TestBuildHelpers:
    def test_model_artifact_loads(self):
        m = TemperatureModel()
        assert "H_C_ratio" in m._m["properties"]
        assert m._m["aromaticity_fit"]["r2"] > 0.8

    def test_compare_models_self_is_zero(self):
        d = compare_models(_MODEL_PATH, _MODEL_PATH)
        assert d  # non-empty
        assert all(v["max_abs_delta"] == 0 for v in d.values())
