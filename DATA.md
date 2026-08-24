# Data in this repository

All files in `data/` are **small derived extracts prepared for teaching**, not
copies of the original databases. Total size ~4 MB. If you use any of them
beyond the summer school, go back to the original source and cite it.

| file | contents | source | licence / terms |
|---|---|---|---|
| `pat_green_inventors_IT.csv.gz` | inventor × patent records for EPO applications with ≥1 Italian inventor and ≥1 green CPC code (Y02/Y04S), priority years 2010–2019 | OECD REGPAT, May 2025 | derived extract, **pseudonymised** (see below); original data available from OECD on request |
| `pat_all_inventors_ITgreen.csv.gz` | as above, plus the non-green EPO applications of the same inventors (used to show how the network boundary changes results) | OECD REGPAT, May 2025 | idem |
| `pat_green_tech_IT.csv.gz` | patent × CPC 4-digit class for the same patents | OECD REGPAT / CPC | derived extract |
| `green_region_tech_EU.csv.gz` | NUTS region × CPC4 × period, fractional patent counts, EU + CH/NO/UK | OECD REGPAT | aggregate counts |
| `green_tech_cooc_EU.csv.gz`, `green_tech_npat_EU.csv` | co-classification counts between CPC4 classes, EU green patents 2015–2019 | OECD REGPAT / CPC | aggregate counts |
| `green_coinvention_country_benchmark.csv` | co-invention network statistics for 17 countries, green patents 2015–2019 | OECD REGPAT | aggregate statistics |
| `cpc4_def.csv` | CPC 4-digit class labels | EPO/USPTO Cooperative Patent Classification | public |
| `cordis_he_projects.csv`, `cordis_he_participants.csv.gz`, `cordis_he_scivoc.csv.gz` | Horizon Europe projects with a 100% climate policy marker or a sustainability euroSciVoc term, and their participants | CORDIS, release 2026-08-06, <https://cordis.europa.eu/dataset> | **CC BY 4.0** |
| `openalex_ce_authorships.csv.gz` | 2,000 works on the circular economy (2020–2024) flattened to work × author × institution; offline fallback for the API block | OpenAlex, <https://openalex.org> | **CC0** |

## Pseudonymisation of the patent files

OECD REGPAT is licensed data and its inventor files carry inventor names and
REGPAT `person_id`s. Because the class needs the *network structure* and not the
identities, in the two inventor files:

* `person_id` was replaced by sequential ids (`I00001`, `I00002`, …), stable
  across both files;
* `inv_name` was replaced by matching labels (`INV_00001`, …).

Everything else — patent id, priority year, NUTS region, country, green flag —
is unchanged, so every degree, betweenness, component and regional indicator
computed in the session is exactly the one you would obtain from the raw data.
The crosswalk back to REGPAT ids is not distributed.

## Citing the sources

* OECD (2025), *OECD REGPAT database, May 2025*, OECD, Paris.
* European Commission, *CORDIS — EU research projects under Horizon Europe*,
  Publications Office of the EU, dataset release 6 August 2026 (CC BY 4.0).
* Priem, J., Piwowar, H., Orr, R. (2022), *OpenAlex: a fully-open index of
  scholarly works, authors, venues, institutions, and concepts*, arXiv:2205.01833.
