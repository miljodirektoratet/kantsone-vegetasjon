"""Region resolution: one name, three national layers it can come from."""

import geopandas as gpd
import pytest
from shapely.geometry import box

from kantsoner.config import Config, RegionPaths, locate_region

CRS = "EPSG:25833"


def _write(path, rows, geoms):
    gpd.GeoDataFrame(rows, geometry=geoms, crs=CRS).to_file(path, driver="GPKG")


@pytest.fixture
def layers(tmp_path):
    vass = tmp_path / "vassdragomr.gpkg"
    _write(
        vass,
        {"vassomrnr": ["001", "002"], "vassomr": ["Foo/Bar", "Baz"]},
        [box(0, 0, 1000, 1000), box(2000, 0, 3000, 1000)],
    )
    vann = tmp_path / "vannomraader.gpkg"
    _write(
        vann,
        {"vannomradeId": ["5101-06"], "navn": ["Glomma sør for Øyeren"]},
        [box(0, 0, 2000, 2000)],
    )
    regine = tmp_path / "regine.gpkg"
    _write(
        regine,
        {
            "elvNavnHierarki": [
                "Alna",
                "Østensjøbekken/Alna",
                "elv fra Aurevatnet/Østensjøbekken/Alna",
                "Storelva",
                "Storelva",
                None,
            ],
            "nedbfeltHavVassdragNr": [
                "006.2Z",
                "006.2Z",
                "006.2Z",
                "010.Z",
                "020.Z",
                "030.Z",
            ],
        },
        [
            box(0, 0, 100, 100),
            box(100, 0, 200, 100),
            box(200, 0, 300, 100),
            box(0, 500, 100, 600),
            box(0, 700, 100, 800),
            box(0, 900, 100, 1000),
        ],
    )
    return {"vassdragomr": vass, "vannomraader": vann, "regineenhet": regine}


def _config(tmp_path, layers, regions, **drop):
    paths = {
        "base_dir": tmp_path / "out",
        "chm": tmp_path / "chm.tif",
        "dtm": tmp_path / "dtm.tif",
        "grunnkart_dir": tmp_path,
        "elvenett": tmp_path / "elv.gpkg",
        "fkb": dict.fromkeys(
            ("bygning_dir", "bygnanlegg_dir", "vann_dir"), tmp_path
        ),
        **{k: v for k, v in layers.items() if k not in drop},
    }
    return Config.model_validate({"paths": paths, "regions": regions})


def test_bare_strings_are_vassdragsomraader_and_mappings_name_a_source(
    tmp_path, layers
):
    cfg = _config(
        tmp_path,
        layers,
        ["Foo/Bar", {"name": "Alna", "source": "regine"}],
    )
    assert cfg.region_names == ["Foo/Bar", "Alna"]
    assert cfg.region_spec("Foo/Bar").source == "vassdragomr"
    assert cfg.region_spec("Alna").source == "regine"
    with pytest.raises(ValueError, match="repeats"):
        _config(tmp_path, layers, ["Foo/Bar", {"name": "Foo/Bar"}])
    with pytest.raises(ValueError, match="empty"):
        _config(tmp_path, layers, [])


def test_resolve_writes_one_aoi_per_source(tmp_path, layers):
    cfg = _config(
        tmp_path,
        layers,
        [
            "Foo/Bar",
            {"name": "Glomma sør for Øyeren", "source": "vannomraader"},
            {"name": "Alna", "source": "regine"},
        ],
    )
    p = RegionPaths.resolve(cfg, "Foo/Bar")
    assert p.folder == "001_Foo_Bar" and p.aoi.exists()
    assert p.area_km2 == pytest.approx(1.0)
    aoi = gpd.read_file(p.aoi)
    assert aoi.iloc[0][["region", "region_nr", "source"]].tolist() == [
        "Foo/Bar",
        "001",
        "vassdragomr",
    ]

    p = RegionPaths.resolve(cfg, "Glomma sør for Øyeren")
    assert p.folder == "5101-06_Glomma_sør_for_Øyeren"
    assert p.area_km2 == pytest.approx(4.0)

    # Every unit with the river anywhere in its hierarchy, dissolved.
    p = RegionPaths.resolve(cfg, "Alna")
    assert p.folder == "006.2Z_Alna"
    aoi = gpd.read_file(p.aoi)
    assert len(aoi) == 1 and aoi.geometry.iloc[0].equals(box(0, 0, 300, 100))
    assert p.area_km2 == pytest.approx(0.03)


def test_regine_name_in_two_catchments_is_refused(tmp_path, layers):
    cfg = _config(tmp_path, layers, [{"name": "Storelva", "source": "regine"}])
    with pytest.raises(ValueError, match="ambiguous"):
        locate_region(cfg, cfg.region_spec("Storelva"))
    cfg = _config(tmp_path, layers, [{"name": "Nope", "source": "regine"}])
    with pytest.raises(ValueError, match="no REGINE unit"):
        locate_region(cfg, cfg.region_spec("Nope"))


def test_source_layer_must_be_configured(tmp_path, layers):
    cfg = _config(
        tmp_path, layers, [{"name": "Alna", "source": "regine"}], regineenhet=1
    )
    with pytest.raises(ValueError, match="regineenhet is not configured"):
        locate_region(cfg, cfg.region_spec("Alna"))
