# Leveransen

Leveransen er et sett nasjonale GeoPackage-filer, én fil og ett lag per
produkt, satt sammen av `scripts/merge_delivery.py` fra regionkjøringene.
Alle filer er i EPSG:25833. Den fullstendige feltlisten står i
[attributter.html](attributter.html).

## Filene

| filnavn | innhold | geometri | nøkkel/ID-felt | viktigste felt |
| --- | --- | --- | --- | --- |
| `official_centerline_10m.gpkg` | Vann-Nett-linjen delt i 10 m-segmenter, med alle målte og beregnede verdier | linje | `segment_id_vann-nett`, `river_id_vann-nett` | `E_S_Kantvegetasjon_A`, `E_S_Kantvegetasjon_grad`, `E_S_Score_total`, `K_S_prioritet`, `K_L/R_B0-10_Vtot` |
| `official_centerline_100m.gpkg` | Ca. 100 m sammenhengende elv | multilinje | `aggr_segment_id_vann-nett_100m` + `region_folder` | `100m_Kantvegetasjon(_final)`, `100m_datadekning`, `100m_Score_total(_final)` |
| `official_centerline_500m.gpkg` | Ca. 500 m sammenhengende elv | multilinje | `aggr_segment_id_vann-nett_500m` + `region_folder` | som 100 m, med `500m_` |
| `official_centerline_1000m.gpkg` | Ca. 1000 m sammenhengende elv | multilinje | `aggr_segment_id_vann-nett_1000m` + `region_folder` | som 100 m, med `1000m_` |
| `official_centerline_river.gpkg` | Én rad per elvevannforekomst | multilinje | `river_id_vann-nett` | `river_*` |
| `official_centerline_vannomraade.gpkg` | Én rad per vannområde | multilinje | `eUVannomradeId`, `vannomradeId` | `vannomraade_*`, `n_rivers` |
| `official_centerline_region.gpkg` | Én rad per beregningsregion (nedbørfelt) | multilinje | `region_id` | `region_*`, `n_rivers` |
| `lake_shore.gpkg` | Innsjøstrand i 10 m-biter, med verdiene for 0–20 m-sonen | flate | `segment_id`, `lake_id` | som 10 m, uten feltene for elvebredde og høyre bredd |
| `lake_polygon.gpkg` | Én rad per innsjø | flate (Vann-Netts innsjøpolygon) | `lake_id` | `innsjo_*`, `innsjo_strandlinje_m` |
| `buffer10.gpkg` | 0–10 m-sonen langs elvene, én flate per segment og bredd | flate | `segment_id_vann-nett`, `side` | `K_B0-10_Vtot`, `E_Kantvegetasjon_A`, `E_S_Kantvegetasjon_grad` |
| `lake_buffer10.gpkg` | 0–10 m-sonen rundt innsjøene, én flate per strandsegment | flate | `segment_id_vann-nett`, `lake_id` | som `buffer10` |
| `single_trees.gpkg` | Trekroner funnet i laserdataene | flate | `tree_id` + `region_folder` | `max_height`, `mean_height`, `crown_area_m2`, `height_class` |

Felles for alle lag:

- `vannforekomstid`, `river_id_vann-nett` og `vannforekomstnavn` er hentet fra
  Vann-Nett-uttrekket fra 2026-10-26, det samme som kjøringen bygde på.
- `eUVannomradeId` er vannområdet som har størst del av objektet.
- `region_folder` er regionen raden kommer fra. ID-ene på 100/500/1000 m og
  trær begynner på nytt i hver region og er bare unike sammen med den.
- Elv som krysser en regiongrense, ble kjørt i begge regionene. Den er levert
  én gang, fra regionen som eier vannforekomsten.

I tillegg ligger det kontrollfiler: `README.md` (lag, antall rader og felt),
`vannforekomst_lines_check.csv` og `vannforekomst_lines_missing.gpkg` (hvilke
deler av Vann-Nett-linjene som ikke er levert, utenom innsjøer),
`lakes_not_delivered.csv` og `lake_shore_double_lines.csv`.

## Rapporteringsnivåene

Alle nivåer over 10 m regnes fra 10 m-segmentene. Bare summer (lengder,
antall, arealer) legges sammen; hver andel og hvert gjennomsnitt regnes på nytt
fra summene på hvert nivå. Andeler er lengdevektet. Klassene
`*_Kantvegetasjon_final` og `*_Score_total_final` står tomme når `datadekning`
er under 50 %; de andre `_final`-klassene fylles alltid.

**10 m.** Grunnenheten: én rad per segment av Vann-Nett-linjen, med alle
målinger og klasser. `E_S_Score_total` og `K_S_prioritet` står tomme på
segmenter som ikke ble beregnet (`E_S_Kantvegetasjon_A` = Uklassifisert, eller
uten kantsone, også lukkede bekker som ikke ble beregnet).

**100 m, 500 m, 1000 m.** Blokker av sammenhengende segmenter innenfor samme
vannforekomst: 10 segmenter blir 100 m, 5 slike blir 500 m, 2 slike blir
1000 m. Blokkene er omtrent så lange, fordi segmentlengden varierer litt.

**Elv (vannforekomst).** Én rad per elvevannforekomst, slått sammen fra alle
dens 10 m-segmenter på tvers av regioner. Innsjøstrand er ikke med.

**Vannområde.** Slås sammen på samme måte, gruppert på `eUVannomradeId`.

**Region.** Hele beregningsregionen (nedbørfeltet) i én rad.

**Innsjøstrand og innsjø.** `lake_shore` er strandsegmentene, målt som en
elvebredd med én side; flaten er 0–20 m-sonen verdiene gjelder. `lake_polygon`
slår sammen strandsegmentene per innsjø med de samme reglene som elvene, så en
innsjøs andeler betyr det samme som en elvs. Doble strandlinjer gjennom åpent
vann (der to vannflater møtes) er fjernet.

**buffer10.** Selve 0–10 m-sonen per bredd (`side` = L eller R), som
målingene ble gjort på. Hver flate har bredden sine verdier og segmentets
`E_S_Kantvegetasjon_fjernet` og `E_S_Kantvegetasjon_grad`, så sonene kan tegnes
etter tilstand uten kobling.

**lake_buffer10.** Det samme rundt innsjøene, én flate per strandsegment,
med verdiene fra `lake_shore`.

**single_trees.** Kronepolygon fra tredeteksjonen. `max_height` er
høydemodellen (CHM) i toppunktet, `mean_height` snittet over kronen og
`height_class` 1–5, 5–10 eller 10pluss m. En krone som to regioner fant, er
levert én gang. Omrisset er forenklet; areal og høyde er målt før det.

### Hvordan verdiene slås sammen

- `Kantvegetasjon`: lengdevektet snitt av `E_S_Kantvegetasjon_fjernet` (andel
  av forventet kantvegetasjon som mangler) over vurderte segmenter.
- `datadekning`: andel av lengden, utenom lukkede bekker, som ble beregnet.
- `Flomvern`, `Forbygning`, `Filter_score_final`, `Struktur_score`,
  `Kontinuitet`, `Score_total`: lengdevektet snitt over beregnede segmenter.
- `Barrier`: antall. `Fragmentering`: den dårligste klassen i enheten.
- `Score_total_final`: 1 (Svært dårlig, under 0,20) til 5 (Svært god, fra
  0,80).

Noen kolonner beregnes, men tas ut rett før filene leveres: andelene etter den
gamle ja/nei-regelen (`Kantvegetasjon_janei`, `Kantskog_fjernet`,
`Nedbygd_kantsone`, `ikke_mulig_Kantvegetasjon`, `Potensiell_Kantvegetasjon`),
dekningstall som gjentar `datadekning` (`frac_matched`, `Uklassifisert`,
`length_uklassifisert`), konstante kolonner (`K_S_Kilde`,
`K_S_beregningsbredde`), tomme kolonner og `segment_id`, som gjentar
`segment_id_vann-nett`. Listen står i `DELIVERY_DROP` i
`scripts/merge_delivery.py`. Innsjøstranden (`lake_shore`) beholder alle sine
kolonner.

## Mellomtrinnet: regionmappene

Hver region kjøres for seg til `<base_dir>/<vassomrnr>_<region>/`:

- `geometry/`: segmentene, rutene gjennom elveløpet, kantsonene per bredd og
  0–3/0–10/0–20 m-sonene.
- `attributes/`: `attributes.gpkg` (målingene) og `derived.gpkg` (alle
  beregnede vurderinger) per 10 m-segment.
- `aggregation/`: alle nivåene med alle mellomregninger.
- `delivery/`: de samme nivåene, bare med leverte felt. Det er disse
  `merge_delivery.py` slår sammen.
