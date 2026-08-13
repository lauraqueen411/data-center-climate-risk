from pathlib import Path

import geopandas as gpd

from climate_risk_dc.exposure.water import load_public_water_sources, plot_public_water_sources


def test_load_and_plot_public_water_sources(tmp_path: Path) -> None:
    gdf = load_public_water_sources()
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert len(gdf) > 0

    output_path = tmp_path / "water_sources.png"
    plot_public_water_sources(gdf, output_path=output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0
