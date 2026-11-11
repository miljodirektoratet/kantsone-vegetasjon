# Metode

Kort beskrivelse av hvordan kantsonene er beregnet. Tallene og reglene står i
koden; formlenes konstanter ligger i `src/kantsoner/spec/constants.py`, og hver
egenskap er listet i `src/kantsoner/spec/attributes.py`.

## 1. Hva beregningen gjør

Landet deles i regioner (NVEs vassdragsområder). Hver region beregnes for seg,
i fire steg.

1. **Geometri** (`src/kantsoner/geometry/`). Utgangspunktet er NVEs elvenett
   (ELVIS) supplert med Vann-Nett; bare segmenter med vannforekomst-id
   leveres. Hver strekning deles i segmenter på
   10 m. Langs strekningen finnes det fysiske elveløpet i vannflatene fra
   Grunnkart og FKB: en linje legges der det er billigst å gå gjennom vann og
   dyrest å gå over land eller oppover i terrenget. Fra vannkanten lages
   buffere på 0–3, 0–10 og 0–20 m. Bufferen deles mellom 10 m-segmentene og i
   venstre og høyre bredd. Hver kantsone vet dermed fra starten hvilket
   segment og hvilken vannforekomst den hører til.
2. **Måling** (`src/kantsoner/primary/`). For hver bredd og buffer måles
   vegetasjon (fra laserdata), helning (fra terrengmodellen), arealtype og
   treslag (fra Grunnkart), bygninger, flom- og erosjonssikring, dammer og
   grøfter. Regionen måles i puljer av nærliggende strekninger.
3. **Formler** (`src/kantsoner/derived/`, `src/kantsoner/network/`). System A
   og System B beregnes per 10 m-segment (avsnitt 2). Her kobles også
   segmentene til et elvenett. Strømretningen er rekkefølgen linjene er
   digitalisert i (elvenettet er digitalisert nedstrøms), ikke
   terrengmodellen. Nettverket brukes til fragmentering, konnektivitet,
   fremmedarter og restaureringsprioritet.
4. **Aggregering** (`src/kantsoner/delivery/aggregate.py`). 10 m-segmentene
   summeres til blokker på 100 m, 500 m og 1000 m langs elva, og til
   vannforekomst, vannområde og region. Blokker inneholder bare
   sammenhengende vann. Reglene står i
   [leveranse.md](leveranse.md#rapporteringsnivåene).

Hvert steg skriver filer til disk, så alle mellomresultater kan åpnes i GIS.

## 2. System A og System B

**System A** følger Vannportalens klassifiseringsveileder for hydromorfologi.
Resultatet oppgis som *grad av endring* fra 5 (svært god, ingen endring) til 1
(svært dårlig).

- **Kantvegetasjon.** For hver bredd beregnes *forventet dekning*: hvor tett
  vegetasjonen ville vært i naturtilstanden, ut fra arealtype og skogbonitet
  i Grunnkart (f.eks. 95 % på skog med svært høy bonitet, 75 % på lav
  bonitet). Andelen fjernet er `(forventet − målt) / forventet`. Segmentet
  får gjennomsnittet av de to breddene. En bredd vurderes bare når minst
  halvparten av arealet kan ha vegetasjon (ikke berg, myr, vann eller snø),
  eller når vegetasjonen alt er så tett som forventet. Klassen *Kantskog fjernet* /
  *Nedbygd kantsone* settes etter hvilken arealbruk som dominerer der
  vegetasjonen mangler (se `derived/system_a.py` og `derived/kantveg_grad.py`).
- **Flomvern og forbygning.** Finnes det sikring i bredden (flomvern: NVE;
  forbygning: NVE og Statens vegvesen)? Hver bredd teller halvparten.
- **Barrierer.** NVE-dammer i drift som krysser elveløpet. Telles per enhet.
- **Fragmentering.** Avstand oppstrøms til nærmeste barriere: under 2 km
  svært dårlig, under 10 km dårlig, under 50 km moderat, ellers god; ingen
  barriere oppstrøms gir svært god.

Andeler (kantvegetasjon fjernet, flomvern, forbygning) omregnes til grad slik:
under 5 % → 5, 5–20 % → 4, 20–50 % → 3, 50–80 % → 2, minst 80 % → 1.
Barrierer: ingen → 5, 1–2 → 3, 3 → 2, 4 eller flere → 1.

**System B** er den funksjonelle metoden (Sállir-rapporten). Fire funksjoner
gir hver en verdi 0–1 per segment.

| Del | Hva den måler | Viktigste inndata |
| --- | --- | --- |
| Filtrering | Hvor mye sediment, fosfor, nitrogen og plantevernmidler kantsonen holder tilbake | Effektiv vegetasjonsbredde, helning, vegetasjonstype (høyde, treslag), arealbruk, grøfter |
| Temperatur | Skygge over vannet | Trehøyde (median og 95-persentil), elvebredde, elvas himmelretning, lauv/bar |
| Næring | Tilførsel av løv til vannet | Vegetasjonsdekning, treslag, elvebredde |
| Kontinuitet | Sammenhengende tre- og buskdekke langs bredden | Dekning av trær (over 5 m) og busker (0,5–5 m) |

De settes sammen slik (`derived/scores.py`):

- Filterscore = vektet snitt av sediment 25 %, fosfor 35 %, nitrogen 25 %,
  plantevern 15 %; justeres for næring og reduseres der grøfter krysser
  bredden og leder vann forbi vegetasjonen.
- Strukturscore = 0,5 × temperatur + 0,5 × næring.
- Kontinuitet = 0,7 × tre + 0,3 × busk.
- **Samlet score** = 0,5 × filter + 0,3 × struktur + 0,2 × kontinuitet.

På 100 m og høyere nivåer gis samlet score også som grad 1–5 i steg på 0,2: under 0,20 → 1 (svært
dårlig), fra 0,80 → 5 (svært god).

**Restaureringsprioritet** (`network/priority.py`) sammenligner bredden
kantsonen har med bredden den burde ha. Anbefalt bredde er 10 m + 1,5 ×
elvebredde (mellom 10 og 50 m), økt for bratt terreng, intensiv arealbruk
(dyrket mark, bebyggelse) og svak konnektivitet. Prioriteten er størst der
manglende bredde er stor, samlet score lav, arealbruken intensiv og
fremmedarter er registrert i eller oppstrøms segmentet.

**Konnektivitet** (`network/connectivity.py`) er snittet av tre ledd: andel
naturlig mark i 0–20 m, samlet score oppstrøms (vektet ned med avstanden,
inntil 10 km) og samlet score innenfor 5 km uansett strømretning.

## 3. Innsjøer

En innsjø som er egen innsjøvannforekomst i Vann-Nett, måles fra strandlinjen,
ikke fra midtlinja elvenettet trekker gjennom den (`geometry/lakes.py`).

- Vann-Netts innsjøpolygon brukes bare til å velge og navngi innsjøen.
  Strandlinjen er Grunnkart/FKBs, fordi Vann-Netts omriss er generalisert.
- Strandlinjen deles i 10 m-segmenter, også rundt øyer. Bufferen er land
  innenfor 20 m av strandlinjen, utenom alt vann.
- Formlene er de samme som for elv, med to unntak: temperatur og næring
  avhenger av elvebredde og beregnes ikke. Samlet score for innsjø er
  0,7 × filter + 0,3 × kontinuitet. Anbefalt bredde er 10 m.
- Barrierer og fragmentering hentes fra innsjøens midtlinje i elvenettet; den
  dårligste verdien gjelder hele innsjøen.
- Leveres som strandsegmenter (`lake_shore`) og én rad per innsjø
  (`lake_polygon`). Midtlinjene i slike innsjøer leveres ikke som elv. En
  innsjø som ikke er egen vannforekomst, regnes som del av elva gjennom den.

## 4. Automatisk gjenkjenning og modellerte data

Noen inndata er ikke målt direkte, men tolket fra laserdata eller hentet fra
datasett som selv er modellberegnet.

- **Enkelttrær** (`primary/singletree.py`). Trekroner finnes i
  vegetasjonshøydemodellen (CHM): høydene glattes, hver topp blir en krone,
  og kronen avgrenses der høyden faller mot nabokronen. Bygninger maskeres
  bort. Kronearealet gir vegetasjonsdekningen og høydeklassene. Trekronene
  leveres som eget lag (`single_trees`).
- **Trehøyde** (`primary/canopy.py`). Median og 95-persentil av CHM-pikslene
  over 1 m i bufferen.
- **Grøfter** (`primary/grofter.py`). Grøftekartet er et landsdekkende,
  maskinelt framstilt kart over sannsynlige grøfter, ikke innmålte grøfter.
  En grøft telles når den krysser bredden (minst 45° vinkel mot elva).
- **Treslag** (`primary/landcover.py`). Lauv, blandingsskog eller barskog fra
  arealdekke-klassene i Grunnkart. Treslag mangler på mange bredder.
- **Elveløpet** (`geometry/elvis.py`). Linjen gjennom vannflatene er beregnet,
  ikke kartlagt. Segmenter der den ligger langt fra Vann-Nett-linja, får
  `Usikkerhet` = Ja.
- **Mulig rørlagt bekk** (`primary/rorlagt.py`). En regel, ikke en måling: en
  rekke segmenter uten kartlagt vann i nærheten og med jordbruk på begge
  sider. Kantvegetasjon vurderes ikke for slike segmenter («Ikke vurdert
  (bekkelukking)»).

## 5. Begrensninger og tolkninger

Valg som er tatt der spesifikasjonen var uklar eller selvmotsigende:

- Skalaen 1–5 har 5 som beste klasse, for alle parametere.
- Kantvegetasjon aggregeres som andel *fjernet*, slik at 0 % er
  naturtilstanden, som for flomvern og forbygning.
- Ingen grad settes for kantvegetasjon eller samlet score når under 50 % av
  enhetens lengde er målt (`datadekning`; rørlagte bekker regnes ikke med).
- Formler som i arbeidsboka blandet venstre og høyre bredd, beregnes likt for
  begge bredder (filterscore, kontinuitet, lateral konnektivitet).
- Grøftevinkelen på 45° gjelder retningen i grader, ikke
  orienteringsfaktoren.
- Fragmentering regnes oppstrøms: segmentet nedenfor en dam er det
  fragmenterte.
- Oppstuving leveres ikke: elvebredde mangler på de fleste innsjømidtlinjer,
  så dammer ville stått uten oppstuving.
- All beregning bruker 0–10 m-bufferen; 0–20 m brukes til effektiv bredde,
  arealandeler og konnektivitet.
- En bredd uten buffer (ingen areal i 0–20 m) regnes som ikke målt og holdes
  utenfor snittene. Et segment uten målt bredd får `Beregnet` = no og klassen
  «Uklassifisert».
- Mangler arealbruk, brukes dyrket mark; mangler treslag under høye trær,
  brukes lauv. Begge markeres i `E_S_usikkerhet_flagg` (bare på 10 m-nivå).
- Mangler helning, brukes 0 %; mangler grøftekart, regnes det som ingen
  grøfter. Ingen av delene markeres.
- Ved riksgrensen måles bare den norske bredden (`Grenseelv`, `Utenfor_Norge`).
- Høydeklassene for vegetasjon er 0,5–1, 1–5, 5–10 og over 10 m.

## 6. Beregningstid

Regionene er uavhengige av hverandre, så `kantsoner run` kjører flere samtidig,
hver i sin egen prosess (`run.workers` i konfigurasjonen, eller `--workers`).
Innenfor en region går de fire stegene etter hverandre, og målingen tar én
pulje (`run.batch_segments` segmenter) om gangen.

En region kjøres som én jobb. Bare geometritrinnet deler den på flere
kjerner: maskinens 32 tråder delt på `run.workers`. De andre stegene bruker
én kjerne. Derfor setter den største regionen en nedre
grense for hvor lang en landsdekkende kjøring blir, uansett hvor mange
regioner som kjøres samtidig.

- Regionene varierer fra 62 til 43 000 km². Tid og minne følger antall
  strekninger mer enn arealet.
- Glomma er størst: rundt 105 000 strekninger og 5,9 millioner
  10 m-segmenter. Den ble kjørt først, alene med `-w 4`: ca. 53 timer og
  opptil 84 GB minne. Med `-w 1` fikk geometritrinnet 32 prosesser og gikk
  tom for minne.
- De øvrige 261 regionene ble kjørt med 12 samtidig (`run.workers: 12`)
  innenfor en minneramme på 100 GB (`run.mem_budget_gb`).
- Hele landet tok ca. 11 døgn (2026-10-29 til 2026-11-09), medregnet en
  omstart av Glomma etter at minnet tok slutt.
- De største regionene startes først, så de ikke blir gående alene til slutt
  (ordnet etter areal, eller etter antall strekninger når
  `run.mem_budget_gb` er satt). Store og små blandes, og en ny region startes
  bare når minst `run.min_free_gb` minne er ledig, så de største ikke tar
  hele maskinen samtidig.

## 7. Ordliste

| Ord | Betydning |
| --- | --- |
| Kantsone | Vegetasjonsbeltet langs elv eller innsjø, her 0–20 m fra vannkanten |
| Vannforekomst | Vann-Netts rapporteringsenhet: en elvestrekning eller innsjø med egen id |
| Vann-Nett | Miljødirektoratets register over vannforekomster |
| ELVIS | NVEs elvenett; Vann-Netts elvelinjer bygger på det |
| Strekning | Én ELVIS-linje; deles i 10 m-segmenter |
| Segment | En 10 m lang bit av en strekning eller strandlinje; minste enhet |
| Buffer | Landstripen målt fra vannkanten, f.eks. 0–10 m |
| Bredd (L/R) | Venstre og høyre side, sett i linjens retning |
| FKB | Kartverkets detaljerte kartdata (vann, bygninger) |
| Grunnkart | Nasjonalt arealdekkekart (arealtype, skogbonitet, treslag) |
| DTM | Terrengmodell: bakkens høyde, fra laser |
| CHM | Vegetasjonshøydemodell: høyden av vegetasjon over bakken, fra laser |
| Grad av endring | Vannportalens skala 1–5 for avvik fra naturtilstand; 5 er best |
| Forventet dekning | Vegetasjonsdekningen bredden ville hatt i naturtilstanden |
| Datadekning | Andel av en enhets lengde som faktisk er målt |
| Beregnet | Om segmentet har minst én målt bredd |
| Region | NVE-vassdragsområde; enheten beregningen kjøres per |
| Pulje | En gruppe nærliggende strekninger som måles sammen, for å spare minne |
| System A / System B | Vannportalens gradering (avsnitt 2) / den funksjonelle scoren 0–1 |
| `_final` | Feltnavn for en grad 1–5 på 100 m og høyere nivåer |
| Lengdevektet | Snitt der hvert segment teller etter sin lengde |
| GeoPackage | Åpent filformat for kartdata (.gpkg), leses av QGIS |
| EPSG:25833 | ETRS89 / UTM sone 33N, Norges vanlige kartprojeksjon |
