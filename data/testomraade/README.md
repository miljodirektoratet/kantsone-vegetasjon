# Testområde Korsrudtjernet

Et lite utsnitt (ca. 1,1 × 0,75 km, Lunner og Gran kommuner, vassdragsområde
012) av alle inndatasettene, for å kontrollere at hele kjeden kjører. Området
har 7 elvestrekninger (blant annet Sløvikelva og en innsjømidtlinje),
innsjøen Korsrudtjernet, jordbruk, skog og bebyggelse.

Kjøres med `config/testomraade.yaml`, se README.md avsnitt 6. Det tar under ett
minutt og bruker under 1 GB minne.

## Slik er utsnittet laget

Hvert datasett er klippet fra den samme filen som den nasjonale kjøringen
brukte (filnavn og datoer som i `config/config.yaml` og
[docs/datakilder.md](../../docs/datakilder.md)):

- **Testområdet** (`vassdragomr/`) er ett rektangel rundt Korsrudtjernet
  (+250 m), lagt inn som et eget vassdragsområde med nummeret 012.
- **Høydemodellene** (`hoydedata/`, CHM og DTM 1 m) dekker området og alle
  strekninger som berører det, pluss 200 m.
- **Vektordata** er hentet for samme utstrekning pluss 1000 m. Flater er
  klippet, linjer og punkter er hele. Grunnkart er lagret som GeoPackage i
  stedet for FileGDB, med samme lag og felt.
- **Kilder uten objekter i området** (NVE-sikringstiltak og dammer, Statens
  vegvesen, fremmede arter) har de tre nærmeste objektene utenfor området,
  så filene og lagene finnes. De påvirker ikke resultatet.
- **FKB** (Bygning, BygnAnlegg og Vann) er ikke åpne data og er derfor ikke med. Filene i `fkb/` har samme navn, lag og felt som originalen (én fil per kommune, EPSG:3044), men ingen objekter. Testområdet kjører derfor uten bygningsfotavtrykk og uten bekkelinjer fra FKB, og resultatet blir litt annerledes enn med de ekte dataene. Den nasjonale kjøringen brukte FKB, se `docs/datakilder.md`.

`omraade.json` har koordinatene til testområdet og utsnittene (EPSG:25833).

Dataene tilhører de respektive leverandørene og følger deres lisenser.
