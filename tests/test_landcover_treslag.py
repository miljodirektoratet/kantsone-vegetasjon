"""Species of tall canopy comes from arealdekkeniva2."""

from kantsoner.spec import constants as C


def test_every_wooded_class_has_a_species():
    table = C.TRESLAG_FROM_AREALDEKKE
    assert len(table) == 15
    for ground in ("skog", "myrTresatt", "jordbrukInnmarksbeiteTresatt"):
        assert table[f"{ground}Lauv"] == "Lauv"
        assert table[f"{ground}Blanding"] == "Blandingsskog"
        assert table[f"{ground}Barblanding"] == "Barskog"
        assert table[f"{ground}Furu"] == "Barskog"
        assert table[f"{ground}Gran"] == "Barskog"
