"""Unit tests for the shared ensemble summary module, on synthetic inputs."""

from __future__ import annotations

import numpy as np
import pandas as pd

from climate_risk_dc.climate.ensemble_summary import (
    compute_period_deltas,
    ensemble_raster_stats,
    ensemble_summary_stats,
)


def test_ensemble_summary_stats_statewide() -> None:
    df = pd.DataFrame(
        {
            "archetype": ["a", "a", "a", "a", "b", "b"],
            "period": ["historical"] * 4 + ["historical"] * 2,
            "value": [1.0, 2.0, 3.0, 4.0, 10.0, 20.0],
        }
    )
    out = ensemble_summary_stats(df, "value", ["archetype", "period"])
    row_a = out[out["archetype"] == "a"].iloc[0]
    assert row_a["mean"] == 2.5
    assert row_a["median"] == 2.5
    assert row_a["n"] == 4
    row_b = out[out["archetype"] == "b"].iloc[0]
    assert row_b["mean"] == 15.0
    assert row_b["n"] == 2


def test_ensemble_summary_stats_by_site() -> None:
    df = pd.DataFrame(
        {
            "facility_id": [1, 1, 2, 2],
            "period": ["historical"] * 4,
            "value": [1.0, 3.0, 5.0, 15.0],
        }
    )
    out = ensemble_summary_stats(df, "value", ["facility_id", "period"])
    assert out.set_index("facility_id").loc[1, "mean"] == 2.0
    assert out.set_index("facility_id").loc[2, "mean"] == 10.0


def test_compute_period_deltas() -> None:
    df = pd.DataFrame(
        {
            "facility_id": [1, 1, 1, 2, 2, 2],
            "gcm": ["m1"] * 3 + ["m1"] * 3,
            "period": ["historical", "midcentury", "endcentury"] * 2,
            "value": [1.0, 1.5, 2.0, 10.0, 9.0, 8.0],
        }
    )
    out = compute_period_deltas(df, "value", ["facility_id", "gcm"])
    row1 = out[out["facility_id"] == 1].iloc[0]
    assert row1["value_delta_midcentury"] == 0.5
    assert row1["value_delta_endcentury"] == 1.0
    row2 = out[out["facility_id"] == 2].iloc[0]
    assert row2["value_delta_midcentury"] == -1.0
    assert row2["value_delta_endcentury"] == -2.0


def test_compute_period_deltas_missing_period_raises() -> None:
    df = pd.DataFrame(
        {"facility_id": [1, 1], "gcm": ["m1", "m1"], "period": ["historical", "midcentury"], "value": [1.0, 2.0]}
    )
    try:
        compute_period_deltas(df, "value", ["facility_id", "gcm"])
        raised = False
    except KeyError:
        raised = True
    assert raised


def test_ensemble_raster_stats() -> None:
    stack = np.array(
        [
            [[1.0, np.nan], [3.0, np.nan]],
            [[2.0, 5.0], [4.0, np.nan]],
            [[3.0, 6.0], [5.0, np.nan]],
        ]
    )
    out = ensemble_raster_stats(stack)
    np.testing.assert_allclose(out["mean"][0, 0], 2.0)
    np.testing.assert_allclose(out["median"][0, 0], 2.0)
    # column (0, 1) has only 2 finite values across the stack (5.0, 6.0)
    np.testing.assert_allclose(out["mean"][0, 1], 5.5)
    # cell (1, 1) is NaN in every GCM -> stays NaN, not fabricated
    assert np.isnan(out["mean"][1, 1])
