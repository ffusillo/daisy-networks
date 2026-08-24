# Network analysis with innovation and sustainability data

**DAISY International Summer School 2026 — hands-on session (90 min)**
Fabrizio Fusillo · University of Turin · Taranto, 7–11 September 2026

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ffusillo/daisy-networks/blob/main/colab/DAISY_networks_colab.ipynb)

Patent data, EU-funded project data and publication data share one feature:
**the relation between two actors is never observed directly**. It is inferred
from co-participation in a document — a patent, a project, a paper. This session
goes from those raw tables to networks, from networks to positional measures, and
from measures to the **network-based indicators** that end up in econometric
models of innovation and sustainability.

Everything runs in R with `igraph`, `data.table` and `ggraph`, on ~4 MB of
prepared extracts included in `data/`.

---

## Start here

### In the browser (Google Colab)
1. Click the badge above.
2. `Runtime → Change runtime type → **R**`.
3. Run the first cell (installs packages, downloads the data). ~1 minute.

### On your machine (VS Code, RStudio, or plain R)
```bash
git clone https://github.com/ffusillo/daisy-networks.git
cd daisy-networks
```
```r
source("00_setup.R")   # installs anything missing, defines the helpers
```
then work through `01_patents.R` → `04_indicators.R` line by line. The scripts
read from the local `data/` folder, so they also work with no internet.

In VS Code: install the **R extension** (REditorSupport) and, ideally,
[`radian`](https://github.com/randy3k/radian) as the R console;
`Ctrl/Cmd+Enter` sends the current line to R.

---

## What is in here

| file | what you build | ~time |
|---|---|---|
| `00_setup.R` | packages, `daisy_data()`, and the three helpers reused everywhere: `proj_two_mode()`, `make_net()`, `giant()` | — |
| `01_patents.R` | OECD REGPAT green patents → **co-invention network** → components, small-worldness, centralities, structural holes → regional indicators | 22 min |
| `02_cordis.R` | CORDIS Horizon Europe climate projects → **organisation and country collaboration networks** → brokers, communities, cross-border ties, tie persistence | 20 min |
| `03_publications.R` | OpenAlex API → **co-authorship, institution and country networks**; aggregate queries as free indicators | 12 min |
| `04_indicators.R` | region × technology matrix → RTA → **relatedness, knowledge space, variety, coherence, complexity, relatedness density** → does relatedness predict entry into new technologies? | 22 min |
| `99_exercises.R` | seven exercises, and the five questions to ask before trusting a network result | — |
| `colab/` | the same code as a Colab notebook (generated from the `.R` files) | — |
| `data/` | the prepared extracts — see [DATA.md](DATA.md) | — |

The four blocks deliberately use **one pipeline**:

```r
long table  →  proj_two_mode()  →  igraph object  →  measures / aggregation
```

`proj_two_mode()` is the whole inferential step: it builds the sparse
event × actor incidence matrix `B` and returns the projection `A = Bᵗ B`, where
`A[i,j]` counts the documents shared by actors *i* and *j*. Change the two column
names and the same function gives you co-inventors, project partners, co-authors,
co-classified technologies, or collaborating countries.

## Requirements

R ≥ 4.2 and `data.table`, `igraph`, `Matrix`, `ggplot2`, `ggraph`, `jsonlite` —
`00_setup.R` installs whatever is missing (and on Linux/Colab points R at the
Posit binary repository, which turns a ten-minute compile into about one minute).
Optional, for two clearly marked chunks: `sna`, `intergraph`, `graphlayouts`.

## Data and licences

Extracts from **OECD REGPAT** (pseudonymised — see [DATA.md](DATA.md)),
**CORDIS** Horizon Europe (CC BY 4.0) and **OpenAlex** (CC0). Code is MIT, the
teaching text is CC BY 4.0. If you reuse any of the data, cite the original
source, not this repository.

## Other sources worth knowing

**PATSTAT** (EPO) for patent citations, families and non-patent-literature
references — the science → technology link; **PatentsView** for disambiguated US
inventors and assignees; **SDC Platinum / Zephyr–Orbis** for alliances, JVs and
M&A; the **Community Innovation Survey** for *declared* cooperation (a rare case
of directly observed ties); **WIOD / TiVA / COMTRADE** with **OECD ANBERD** for
trade and embodied-R&D networks. Toy datasets for practising network mechanics:
<https://toreopsahl.com/datasets/>.
