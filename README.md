# Kantsone vegetasjon 2026

- **Prosjekteier:** Miljødirektoratet
- **Oppdragstaker:** Field Geospatial AS
- **Leveringsdato:** 2026-11-11
- **Nedlastingslenke:** <https://dataleveranser.miljodirektoratet.no/nedlasting/089a5626-c799-4465-871f-66da5bd6afaf/>

## Kort beskrivelse

**English.** This repository contains the code that produced the national
riparian zone (*kantsone*) dataset for Norwegian rivers and lakes. It measures
vegetation, terrain, land cover and physical interventions in the 0–20 m zone
along every Vann-Nett river stretch, from airborne laser data, maps and
registers, and scores them with two methods. The results are delivered as
GeoPackages at 10 m, 100 m, 500 m, 1000 m, water body, water area and catchment
level.

**Norsk.** Dette repositoriet inneholder koden som laget det nasjonale
datasettet for kantsoner langs elver og innsjøer i Norge. Den måler
vegetasjon, terreng, arealdekke og inngrep i sonen 0–20 m fra vannkanten langs
hver strekning i Vann-Nett, fra laserdata, kart og registre, og vurderer dem
etter to metoder. Resultatene leveres som GeoPackage-filer på 10 m, 100 m,
500 m, 1000 m, vannforekomst, vannområde og vassdragsområde.

## 1. Bidragsytere

| Navn | Organisasjon | E-post |
| --- | --- | --- |
| Luis Barreiro | Field Geospatial AS | luis.barreiro@fieldgeo.no |
| Martijn Vermeer | Field Geospatial AS | martijn.vermeer@fieldgeo.no |
| Ourania Tsilomitrou | Field Geospatial AS | ourania.tsilomitrou@fieldgeo.no |

Eksterne samarbeidspartnere: Sállir Natur AS og NIVA.

## 2. Formål og bakgrunn

Miljødirektoratet trenger et landsdekkende og etterprøvbart bilde av tilstanden
i kantsonene langs vassdragene. Prosjektet beregner dette automatisk for hele
landet, knyttet til Vann-Netts vannforekomster, slik at resultatene kan brukes
direkte i vannforvaltningen.

Metoder:

- **System A** – Vannportalens klassifisering av hydromorfologi:
  kantvegetasjon, flom- og erosjonssikring og vandringshindre, som *grad av
  endring* på skalaen 1–5 (5 er best).
- **System B** – funksjonell metode: filtrering, temperatur, næringstilførsel
  og kontinuitet, samlet til én score og en restaureringsprioritet.
- Begge beregnes per 10 m-segment og summeres opp til de høyere nivåene.

Referanser:

- [Vannportalen, klassifiseringsveileder for hydromorfologi i elv](https://www.vannportalen.no/veiledere/klassifiseringsveileder/klassifisering-av-okologisk-tilstand/klassifisering-av-okologisk-tilstand--i-elv/hydromorfologiske-kvalitetselement-i-elv/)
- Produktspesifikasjonen *Kantsoner Morfologi* (regneark). Den forvaltes av
  prosjektet og ligger ikke i repositoriet.
- [docs/metode.md](docs/metode.md): hvordan beregningen gjøres, valg og
  begrensninger.

## 3. Forutsetninger

**Maskinvare** (maskinen leveransen ble beregnet på):

| | |
| --- | --- |
| Operativsystem | Ubuntu 24.04.4 LTS (Linux-kjerne 7.0) |
| CPU | Intel Core i9-14900K, 24 kjerner / 32 tråder |
| RAM | 128 GB |
| GPU | NVIDIA GeForce RTX 4090 24 GB (driver 580.173), ikke brukt av beregningen |
| Lagring | lokal NVMe SSD 3,6 TB for resultatene, nettverksdisk for de nasjonale høydemodellene |

**Ytelse og kjøretid.** Regionene (NVEs vassdragsområder) er uavhengige og
kjøres flere samtidig, én prosess per region (`run.workers` eller
`--workers`). En region kjøres som én jobb, så den største bestemmer total
tid. Glomma (ca. 105 000 strekninger) ble kjørt først, alene med 4 samtidige
prosesser: ca. 53 timer og opptil 84 GB minne. De øvrige 261 regionene ble
kjørt med 12 samtidig innenfor en minneramme på 100 GB. Total kjøretid for
hele landet var ca. 11 døgn (2026-10-29 til 2026-11-09), medregnet en
omstart av Glomma etter at minnet tok slutt.

Mer om dette i [docs/metode.md, avsnitt 6](docs/metode.md#6-beregningstid).

**Kjøremiljø.** Lokal PC med Linux, uten container. Ingen
systembiblioteker trengs utover det Python-pakkene har med seg: rasterio og
pyogrio inneholder GDAL 3.12 og PROJ.

**Programvare.** QGIS (en nyere versjon) for å se på GeoPackage-filene.
[uv](https://docs.astral.sh/uv/) ≥ 0.13 for Python-miljøet.
**Programmeringsspråk:** Python 3.12 (`.python-version`, `requires-python >=3.12`).

**Viktigste pakker** (låste versjoner i `uv.lock`): geopandas 1.1.4,
shapely 2.1.2, rasterio 1.5.1, pyogrio 0.13.0, scikit-image 0.26.0,
scipy 1.18.0, pandas 3.0.5, numpy 2.5.2.

**Miljødefinisjon.** Definert av `pyproject.toml` og `uv.lock`; installeres
nøyaktig med `uv sync --locked`.

**KI-verktøy.**

- Claude Code (Anthropic) er brukt både til kode (generering og
  gjennomgang) og til dokumentasjon. Modellene som står i git-historikken er Claude Opus 5.5,
  Opus 5, Sonnet 5.5, Sonnet 5, Sonnet 4.6, Sonnet 4.5 og Fable 5.1.
- Det er ikke brukt egne instruksjonsfiler (CLAUDE.md/AGENTS.md).

## 4. Struktur

```
config/config.yaml      stier til inndata og kjøreinnstillinger (plassholdere)
config/testomraade.yaml samme for testområdet i data/
src/kantsoner/          selve beregningsprogrammet
  cli.py, stages.py       kommandoen `kantsoner` og de fire stegene
  io.py                   all lesing og skriving av data
  geometry/               10 m-segmenter, elveløp, buffere, innsjøer
  primary/                målinger fra laserdata, kart og registre
  derived/                formlene for System A og System B
  network/                elvenettet: hindre, konnektivitet, prioritet
  delivery/               aggregering til 100 m, 500 m, 1000 m, vannforekomst, vannområde og region
  spec/                   attributtlisten og metodekonstantene
scripts/                verktøy før og etter kjøringen (se avsnitt 6)
docs/                   dokumentasjon (se avsnitt 5)
tests/                  automatiske tester, blant dem hele kjeden på testområdet
data/testomraade/       lite utsnitt av alle inndata (ca. 20 MB), se README der
```

Resultatene leveres separat via nedlastingslenken øverst. Kjøringer av
testområdet skriver til `outputs/`, som ikke er en del av repoet.

## 5. Input og output

- **Input:** nasjonale høydemodeller fra laser (CHM og DTM, 1 m), Grunnkart
  for arealanalyse, FKB, NVEs elvenett (ELVIS) supplert med Vann-Nett, NVEs
  sikringstiltak og dammer, Vann-Netts vannforekomster m.m. Full liste med
  versjoner: [docs/datakilder.md](docs/datakilder.md).
- **Output:** én GeoPackage per produkt i EPSG:25833: elvesegmenter på 10 m,
  blokker på 100/500/1000 m, vannforekomst, vannområde, region, innsjøstrand,
  innsjø, 0–10 m-sonene og enkelttrær. Se
  [docs/leveranse.md](docs/leveranse.md).
- **Felt:** alle attributter med forklaring i
  [docs/attributter.html](docs/attributter.html).
- **Dataflyt:** fra kilde til leveranse i
  [docs/dataflyt.html](docs/dataflyt.html).

## 6. Installasjon og kjøring

### Kontroller at alt kjører (testområdet)

Uten de nasjonale dataene kan hele kjeden kjøres på det lille testområdet i
`data/testomraade/`. FKB er ikke åpne data, så FKB-filene der er tomme (samme lag og felt, ingen objekter), se `data/testomraade/README.md`. Det tar under ett minutt:

```bash
uv sync --locked --group test
uv run --group test pytest -q tests/test_integrasjon.py
```

Testen kjører de samme kommandoene som under, på testområdet, og sjekker at
alle leveransefilene blir laget og har beregnede verdier. Vil du se på
resultatet i QGIS, kjør stegene selv:

```bash
KANTSONER_EXPORT_SINGLETREE=1 uv run kantsoner run -c config/testomraade.yaml
uv run python scripts/merge_delivery.py -s outputs/testomraade -o outputs/leveranse \
    --work outputs/merge -w 1 --min-free-gb 0 --stop-below-gb 0 \
    --vannomrader data/testomraade/vannomraader/vannomraader_miljoedir_25833_20260227.gpkg \
    --vannforekomster data/testomraade/vannnett/vannforekomster_20261026.gpkg
```

### Hele landet

1. Installer miljøet:

   ```bash
   uv sync --locked
   ```

2. Klargjør elvenettet fra Vann-Nett-nedlastingen og NVEs ELVIS og REGINE:

   ```bash
   uv run python scripts/vannnett_from_geonorge.py \
       vannforekomster_0000_norge_25833_FILEGDB.gdb <DATA>/vannnett 20261026
   uv run python scripts/patch_elvis.py \
       --elvis <DATA>/elvenett/elvenett_nve.gpkg \
       --vannnett <DATA>/vannnett/elvevannforekomster_20261026.gpkg \
       --regine <DATA>/nve_felt/regineenhet_nve.gpkg \
       --out <DATA>/elvenett/elvenett_vn_patched_20261106.gpkg
   ```

   Kommandoen supplerer og kompletterer i én kjøring. Det leverte nettet ble
   laget i to omganger (supplert 2026-10-27, så kompletert med `--complete`),
   så en ny kjøring kan gi en litt annen fil.

3. Bytt ut plassholderne `<DATA>` og `<OUT>` i `config/config.yaml` med
   lokale mapper. Kjøreinnstillingene (`run:`) er de som ble brukt.

4. Sjekk at alle kilder kan åpnes (tar sekunder, ikke timer):

   ```bash
   uv run kantsoner config check
   ```

5. Kjør alle regioner. Med `KANTSONER_EXPORT_SINGLETREE=1` lagres
   trekronene, og hver region får `delivery/single_trees.gpkg`, som trengs
   for `single_trees`:

   ```bash
   export KANTSONER_EXPORT_SINGLETREE=1
   uv run kantsoner run -w 4 -r "Glommavassdraget/Hvaler og Singlefjorden"
   uv run kantsoner run -c config/uten_glomma.yaml
   ```

   `uten_glomma.yaml` er en kopi av `config/config.yaml` der Glomma er tatt
   ut av `regions:`.

   Enkeltsteg og enkeltregioner: `uv run kantsoner --help`.
6. Sett sammen leveransen fra regionene, deretter enkelttrærne:

   ```bash
   uv run python scripts/merge_delivery.py -s <OUT> -o <LEVERANSE> --work <TMP> \
       --vannomrader <DATA>/geonorge/natur/vannomraader_miljoedir_25833_20260227.gpkg \
       --vannforekomster <DATA>/vannnett/vannforekomster_20261026.gpkg
   uv run python scripts/merge_delivery.py <samme argumenter som over> --products single_trees \
       --chm <DATA>/hoydedata/chm_norway_1m_epsg25833_20260911.tif
   ```

En ny kjøring gir samme metode, men tallene kan avvike litt fra de leverte
dataene: der to strekninger er like nære samme vannflate, kan fordelingen av
vannet og rekkefølgen på radene bli en annen, og det kan oppstå små
restflater.

Alle tester: `uv run --group test pytest -q`.

## 7. Lisens

MIT. Se [LICENSE](LICENSE).
