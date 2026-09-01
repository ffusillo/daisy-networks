# Network analysis: indicators, mapping and empirics with innovation data
### DAISY International Summer School 2026 — hands-on session
Fabrizio Fusillo (University of Turin) · Taranto, 7–11 September 2026

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ffusillo/daisy-networks/blob/main/colab/DAISY_networks_colab.ipynb)

Networks for innovation data such as patent data, R&D project and
publication data share one feature: the relation between actors is
not observed directly, it is *inferred* from co-participation to an event.
This session goes from those raw tables to networks, to positional measures,
and finally to the network-based indicators that end up in network analysis and econometric models.

---

## Contents

| file | what it does |
|---|---|
| `00_setup.R` | packages, data access, and the helper functions reused everywhere |
| `01_patents.R` | OECD REGPAT green patents → co-invention network → inventor positions → regional indicators |
| `02_cordis.R` | CORDIS Horizon Europe → organisation & regional collaboration networks |
| `03_publications.R` | OpenAlex API → co-authorship and institution networks |
| `04_indicators.R` | networks as *measurement devices*: relatedness, knowledge space, variety, coherence, complexity, diversification |
| `05_trade.R` | OECD TiVA + CEPII BACI: observed, directed, **valued** flows → filtering as method, trade blocs, brokerage between regions, gross vs value added, the product space |
| `06_brokerage_communities.R` | brokerage (Gould–Fernandez roles) and community detection choices (algorithms, resolution, stability, consensus, reporting) |
| `99_exercises.R` | home exercises + the questions checklist |
| `CODEBOOK.md` | one table per data source: variables, what you can build, pitfalls |
| `DATA.md` | provenance and licence of every file in `data/` |
| `colab/DAISY_networks_colab.ipynb` | the same code as a Colab notebook (R runtime), one result per cell |
| `data/` | pre-processed extracts used in class |

## Running the session

### Option A — VS Code (or RStudio), locally
```r
# Download scripts and the data folder and put all into a "lesson" folder
# set the folder "lesson" as your working directory, then
source("00_setup.R")     # installs what is missing
# and run 01_patents.R ... 06_brokerage_communities.R line by line
```
Requirements: R ≥ 4.2 with `data.table`, `igraph`, `Matrix`, `ggplot2`, `ggraph`,
`jsonlite` (installed automatically by `00_setup.R`). In VS Code use the
**REditorSupport R extension** + `radian`, and keep the interactive terminal open
(`Ctrl/Cmd+Enter` sends the current line to R).

### Option B — Google Colab
1. Open `colab/DAISY_networks_colab.ipynb` in Colab.
2. `Runtime → Change runtime type → R`.
3. Run the **first code cell** (packages + helpers + data access) and wait: it
   installs the packages, which takes a few minutes even from the binary
   repository. Everything after it is split one result per cell, so you can read
   the notebook by scrolling.



## Data sources used

| source | what it is | access | licence |
|---|---|---|---|
| **OECD REGPAT** (May 2025) | EPO/PCT applications with disambiguated inventors and applicants, geocoded to NUTS-3 / TL3 | on request to OECD (free, academic) | OECD terms, redistribution of raw data not allowed |
| **CPC / Y02 tagging** | "green" technology identification (Y02A–Y02W, Y04S) | in PATSTAT/REGPAT | — |
| **CORDIS** Horizon Europe | all EU-funded projects, participants, contributions, euroSciVoc topics, policy markers | bulk CSV, `https://cordis.europa.eu/dataset` | CC-BY 4.0 |
| **OpenAlex** | ~250M works, authors, institutions, topics, citations | REST API + snapshots, no key | CC0 |
| **OECD TiVA** | value added of country *i* embodied in the final demand of country *j*; directed and valued | SDMX API | OECD terms, reuse with attribution |
| **CEPII BACI** | bilateral gross exports by HS6 product | direct download from CEPII | free, citation required |

The extracts in `data/` are small derived aggregates: see [DATA.md](DATA.md) for
the provenance of each file and [CODEBOOK.md](CODEBOOK.md) for its variables.
**Raw REGPAT files are not redistributed**, and the two inventor files are
pseudonymised: the class needs the network structure, not the identities.

Code is MIT; the teaching text is CC BY 4.0; the data keep the licences of their
original sources. If you reuse any of the data, cite that source rather than this
repository.

## Other sources worth knowing
* **PATSTAT** (EPO) — the reference patent database; citations (incl. non-patent
  literature: the science→technology link), families, legal status.
* **PatentsView / USPTO** — disambiguated US inventors, assignees, locations; API + bulk.
* **Google Patents Public Data** (BigQuery) — full text, CPC, citations.
* **SDC Platinum / Zephyr–Orbis** — strategic alliances, JVs, M&A.
* **Community Innovation Survey (CIS)** — declared cooperation partners (direct ties).
* **OECD TiVA, WIOD, COMTRADE** — trade / GVC networks; combined with ANBERD R&D
  data they give embodied-R&D flow networks.
* **Toy datasets** for teaching: <https://toreopsahl.com/datasets/>
