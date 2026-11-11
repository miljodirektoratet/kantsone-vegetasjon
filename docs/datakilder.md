# Datakilder

Alle inndatasett som kantsonene er beregnet fra. Stiene står i
`config/config.yaml` (`paths.*`). Unntaket er de fem siste radene, som bare
skriptene i `scripts/` leser. Alle data er i EPSG:25833. Kilder i et annet
koordinatsystem (FKB er i EPSG:3044) transformeres når de leses inn.
`kantsoner config check` sjekker at alle kildene kan åpnes.

| Datasett | Leverandør | Versjon / uttrekk | Format | Hvordan det hentes | Brukes til |
|---|---|---|---|---|---|
| Vassdragsområder (`vassdragomr`) | NVE | 2026-02-27 | GeoPackage | Geonorge, NVE Nedbørfelt (`nedborfelt_nve_25833_20260227.gpkg`), basert på laget Vassdragsomr; 262 områder, felt omdøpt til vassomrnr/vassomr | avgrensning og nummer for hvert område |
| Vegetasjonshøyde, CHM 1 m (`chm`) | Kartverket / Høydedata | 2026-09-11 | GeoTIFF | hoydedata.no | trehøyde, tredekning, enkelttrær |
| Terrengmodell, DTM 1 m (`dtm`) | Kartverket / Høydedata | 2026-09-11, hull fylt | GeoTIFF | hoydedata.no | helning, rutingen av elveløpet, bakkehøyde ved trærne |
| Grunnkart for arealanalyse (`grunnkart_dir`) | SSB / NIBIO | årsversjon 2025 | FileGDB, nasjonal | Geonorge | elve- og innsjøflater, arealdekke, bebygd areal |
| Elvenett, ELVIS supplert med Vann-Nett (`elvenett`) | NVE, bearbeidet i prosjektet | 2026-11-06 | GeoPackage | lages med `scripts/patch_elvis.py`, se under | elvegeometrien som 10 m-segmentene kuttes fra |
| FKB-Bygning, FKB-BygnAnlegg (`fkb.bygning_dir`, `fkb.bygnanlegg_dir`) | Kartverket / Geovekst | nedlastet 2026-04-28 (data oppdatert til 2025-07) | GeoPackage per kommune | Geonorge, SOSI konvertert til GeoPackage | bygninger: tas ut av tredeteksjonen, gir Bygg-dekning |
| FKB-Vann (`fkb.vann_dir`) | Kartverket / Geovekst | nedlastet 2026-04-28 (data oppdatert til 2025-07) | GeoPackage per kommune | Geonorge, SOSI konvertert til GeoPackage | bekker som linjer, vannflater som elveløpet legges gjennom |
| Sikringstiltak: flom- og erosjonssikring (`nve.*sikring_*`) | NVE | uttrekk 2026-02-05 | GeoPackage, linjer og punkter | NVE (sikringstiltak, Geonorge-fil 2026-02-27) | Flomvern, Forbygning |
| Damlinjer og dampunkter (`nve.dam_*`) | NVE | uttrekk 2026-02-05 | GeoPackage | NVE (Vassdragsanlegg) | vandringshindre |
| Plastring/erosjonssikring (`vegvesen`) | Statens vegvesen | eksport 2026-04-29 | GeoPackage | NVDB, eksport fra Vegkart | Forbygning, sammen med NVE |
| Grøfter (`grofter`, valgfri) | Field Geospatial | 2025-11-28 | GeoPackage | modellprediksjon, se under | grøfteattributtene |
| Fremmede arter (`fremmedart`, valgfri) | Artsdatabanken, samlet i prosjektet | v01, uttrekk 2026-08-26 | GeoPackage | Artskart (artskart.artsdatabanken.no), ett uttrekk per art, slått sammen | fremmedart_risiko, K_S_d |
| Elvevannforekomster (`elvevannforekomster`) | Miljødirektoratet, Vann-Nett | 2026-10-26 | GeoPackage | Geonorge, konvertert med `scripts/vannnett_from_geonorge.py` | vannforekomstid per segment |
| Innsjøvannforekomster (`vannforekomster`) | Miljødirektoratet, Vann-Nett | 2026-10-26 | GeoPackage | samme nedlasting | hvilke innsjøer som leveres som innsjø |
| NHM_ProsjektDekning (`nhm_prosjektdekning`) | Kartverket | 2026-09-11 | Shapefile | hoydedata.no, metadata levert med høydemodellene | skanneår (K_S_Laser_aar) |
| Administrative enheter, riksgrensen (`landegrense`) | Kartverket | 2026-02-27 | GeoPackage | Geonorge, Administrative enheter (fylker), slått sammen til én flate for Norge | klipping ved riksgrensen, Utenfor_Norge |
| ELVIS elvenett, original | NVE | uttrekk 2026-02-05 | GeoPackage | NVE (nedlasting.nve.no) | inndata til `scripts/patch_elvis.py` |
| REGINE-enheter | NVE | 2026-02-27 | GeoPackage | Geonorge, NVE Nedbørfelt (`nedborfelt_nve_25833_20260227.gpkg`), laget RegineEnhet | elvenavn i `scripts/patch_elvis.py` |
| Vannområder | Miljødirektoratet | 2026-02-27 | GeoPackage | Geonorge | `merge_delivery.py`: vannområde per segment (eUVannomradeId) |
| Vannforekomster, alle lag | Miljødirektoratet, Vann-Nett | 2026-10-26 | GeoPackage | Geonorge, samme nedlasting som over, FileGDB konvertert til GeoPackage med ogr2ogr | `merge_delivery.py`: vannforekomstid, elvenavn og innsjø-id i leveransen |
| CHM 1 m | Kartverket / Høydedata | 2026-09-11 | GeoTIFF | som over | `merge_delivery.py`: middelhøyde per trekrone i single_trees |

## Merknader

**CHM og DTM.** Disse er de nasjonale høydemodellene med 1 m oppløsning fra
laserskanning, uttrekk 2026-09-11. Hull i DTM er fylt før bruk (filen heter
`_filled`). Vegetasjonshøyden leses direkte fra CHM-filen.

**Elvenett (`elvenett_vn_patched_20261106.gpkg`).** Dette er NVEs ELVIS-elvenett
med linjene fra Vann-Nett lagt til der ELVIS mangler dem.
`scripts/vannnett_from_geonorge.py` gjør Vann-Nett-nedlastingen fra Geonorge
om til de to filene pipelinen leser. Deretter lager `scripts/patch_elvis.py`
det supplerte nettet av ELVIS, elvevannforekomstene og REGINE-enhetene (som
gir elvenavnene).

**Grøfter.** Grøftene er en modellprediksjon som Field Geospatial har laget fra
flybilder (RGB) og DTM, og som er gjort om til midtlinjer. Det er foreløpig ikke
et offentlig register. Laget er valgfritt. Mangler det, regnes det som ingen
grøfter (`Grofte_tetthet` = 0, grøfteretningen står tom).

**Fremmede arter.** Funn av parkslirekne, kjempespringfrø og vasspest fra
Artsdatabanken, samlet i én fil. Laget er valgfritt.

**Vann-Nett.** Kjøringen og sammenstillingen må bruke samme nedlasting
(2026-10-26). Ellers stemmer ikke vannforekomst-id-ene, fordi
Miljødirektoratet deler og nummererer om vannforekomster mellom utgavene.
