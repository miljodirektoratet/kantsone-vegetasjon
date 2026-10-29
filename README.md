# [Prosjektnavn]-[Årstall]

_Vennligst fyll ut denne README.md-malen i henhold til instruksjonene. Gjør prosjektspesifikke tilpasninger der det er nødvendig. Fjern alle instruksjoner og kommentarer når du er ferdig._

- **Norsk tittel**: [Prosjektnavn]-[Årstall]
- **Prosjekteier**: Miljødirektoratet
- **Oppdragstaker**: [Oppdragstaker]
- **Leveringsdato**: [Leveringsdato]
- **Nedlastingslenke:** [dataleveranser.miljodirektoratet.no/nedlasting/leveranse-id](https://dataleveranser.miljodirektoratet.no/nedlasting/<leveranse-id>) _Fyll inn URL til nedlastingssiden for leveransen i Dataleveranse-applikasjonen (kun for offentlige leveranser)._

Kort beskrivelse av formål og innhold i repoet.

- **English**:
- **Norwegian**:

Eksempel: "Beregninger av bestandsindekser og flerartsindekser i Norsk hekkefuglovervåking 2026"

## 1. Bidragsytere / Kontakt

- **Navn**: [Navn på bidragsyter(e)]
- **Organisasjon**: [Navn]
- **Epost**: [e-post]
- **Eventuelle eksterne samarbeidspartnere:**

## 2. Formål / Bakgrunn

- Kort om overvåkingsprogrammet
- Hvilke metoder som brukes
- Relevante referanser

## 3. Forutsetninger / Dependencies

_Dokumenter alle forutsetninger som kreves for å kjøre analysen og gjenskape resultatene._

- **Hardware**:
  - **OS:** [name, version/build]
  - **Memory (RAM):** [total GB, speed]
  - **CPU:** [model, total cores/threads]
  - **GPU:** [model, VRAM, driver/CUDA]
  - **Storage:** [type, size]
  - **Ytelse:** _Gjerne spesifiser forventet kjøretid eller benchmarktester for analysene._
- **Kjøringsmiljø:**
  - **Container:** _Dockerfile, docker compose, andre relevante konfigfiler._
  - **Lokal PC:** _OS og version/build, viktige systemavhengigheter (f.eks. GDAL) og eventuelle andre avhengigheter som kreves av Python-/R-pakker._
  - **Databricks:** _Databricks Runtime-versjon, eventuelle cluster init-scripts eller policykrav._
  - **Andre kjøringsmiljø:** _Beskriv relevante konfigurasjoner og systemkrav._
- **Programvare:**
  - [Programvare (versjon ≥ x.x)]
  - _Dokumenter alle nødvendige programvareavhengigheter og versjoner (QGIS, ArcGIS Pro, E-cognition, FME, osv.)_
- **KI-verktøy:**
  - [KI-verktøy (modell + versjon ≥ x.x)]
  - _Dokumenter alle KI-verktøy og versjoner som brukes (ChatGPT, Claude, osv.)_
  - _Alle relevante instruksjonsfiler (f.eks `AGENTS.md` eller `CLAUDE.md`) legges ved i repoet._
  - _Beskriv kort hvordan KI-verktøyene er brukt (f.esk. dokumentasjon, analyse, generering av kode)._
- **Programmeringsspråk:**
  - [Programmeringsspråk (versjon ≥ x.x)]
  - _Dokumenter alle programmeringsspråk og versjoner som brukes (Python, R, SQL, osv.)_
- **Definisjon av analysemiljø:**
  - Alle nødvendige pakker, avhengigheter og versjoner skal spesifiseres i en miljøfil som er lagret i repoet.
  - Miljøfilen skal være tilstrekkelig til at analysemiljøet kan gjenskapes uten manuell installasjon eller spesifikasjon av pakker og versjoner.
  - `renv.lock` (_anbefalt for R_), `pyproject.toml` (_anbefalt for Python_), `environment.yml`, `pixi.toml` eller `requirements.txt`.

  **R-miljø:**
  - Vi anbefaler å bruke `renv` for å håndtere pakker og avhengigheter.
  - Eventuelle systemavhengigheter som kreves for å bygge eller kjøre pakker (f.eks. GDAL, PROJ og GEOS for geospatiale pakker), skal dokumenteres.
  - Dersom installasjon av systemavhengigheter krever egne kommandoer, anbefaler vi å legge ved et skript (f.eks. `.sh`) som dokumenterer installasjonstrinnene.
  - R versjon skal være >4.4

  **Python-miljø:**
  - Vi anbefaler å bruke `uv` som pakkebehandler og definere miljøet i `pyproject.toml`.
  - Dersom conda-støtte er nødvendig, anbefaler vi å bruke `pixi`.
  - Se [uv-demo](https://github.com/miljodirektoratet/uv-demo)-repoet for anbefalt praksis for Python-utvikling og reproducerbare miljøer.
  - Python versjon skal være >3.11

  **ArcPy-miljø:**
  - For analyser som bruker ArcPy, anbefales ArcGIS Pro sitt conda-miljø.
  - Miljøet bør eksporteres til en `environment.yml`-fil som lagres i repoet.
  - Se [gis-utils-public](https://github.com/miljodirektoratet/gis-utils-public) for anbefalt praksis for ArcPy-utvikling og reproducerbare miljøer.

## 4. Struktur

Dette repoet følger en standard mappestruktur for å gjøre det enkelt å navigere

- _/config_ - Konfigurasjonsfiler, parameterfiler.
- _/notebooks_ - Python eller Rmd Notebooks
- _/src_, _/R_, _/sql_ - Alle kodefiler
- _/data_ - Eksempeldata eller små testfiler (Dersom hensiktsmessig), samt metadata
- _/docs_ - Metodikk, rapporter, dokumentasjon
- _/outputs_ - Resultater (figurer, tabeller, rapporter)

Oppdragstaker kan opprette andre nødvendige filer (som feks. .gitignore etc.) selv.

## 5. Input / Output

- **Input**: Beskriv hvilke filer som trengs, filformat og eksempeldata. Dersom data hentes fra eksterne kilder, skal versjon, uttrekkstidspunkt og framgangsmåte for å hente ut de samme dataene dokumenteres.
- **Output**: Beskriv hva skriptene genererer, f.eks. figurer, tabeller, filer eller rapporter

## 6. Installasjon / Kjøring

- Kort forklaring på hvordan og i hvilke rekkefølge script(ene) kjøres og miljøet settes opp.

## 7. Lisens / Rettigheter

Se [**LICENSE**-filen](LICENSE) for mer informasjon.
​‌
