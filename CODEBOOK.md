# Codebook
### Network analysis: indicators, mapping and empirics with innovation data
DAISY International Summer School 2026

Five sources, ten networks. Every file in `data/` is a small teaching extract:
the tables below list what is inside, what you can build with it, and what will
bite you. Sizes and counts refer to the extracts shipped with the session, not
to the full databases.

**Read the "what you can build" rows as a menu**: each source supports several
levels of analysis (individual, organisation, region, country, category), and
each level is a different research question — not a robustness check.

Provenance and licences: [DATA.md](DATA.md). Loading is always through
`daisy_data("<file>")`, defined in `00_setup.R`.

---

## 1. Patents — OECD REGPAT (May 2025), EPO applications

Unit: **one inventor on one patent application**. Green = at least one CPC code
in Y02/Y04S (climate-change mitigation and adaptation). The extract keeps every
application with at least one Italian inventor, priority years 2010–2019, plus
the full inventor team of those applications (so foreign co-inventors are in).

### `pat_green_inventors_IT.csv.gz` — 12,149 rows · 4,590 patents · 10,069 inventors

| variable | type | description | notes |
|---|---|---|---|
| `appln_id` | integer | patent application id (the **event**) | unique within REGPAT; not a publication number |
| `prio_year` | integer | priority year, 2010–2019 | closest date to the act of invention; use it, not the publication year |
| `person_id` | character | disambiguated inventor id (the **actor**) | pseudonymised in the public extract (`I00001`…) |
| `inv_name` | character | inventor label | pseudonymised (`INV_00001`…); the instructor copy holds real names |
| `reg_code` | character | NUTS-3 / TL3 region of residence | `substr(reg_code, 1, 4)` gives NUTS-2 |
| `ctry_code` | character | country of residence, ISO-2 | 53 countries: the team, not only Italy |

**What you can build**
- co-invention network of **inventors** (projection on `appln_id`) → degree, strength, betweenness, eigenvector, Burt constraint, brokerage;
- co-invention network of **regions** (projection on `nuts2`) → which regions co-patent, cross-regional openness;
- inventor **productivity** (patents per `person_id`), team size distribution, entry/exit of inventors over time;
- region-level regressors: average positions, share of connected inventors, share of extra-regional ties.

**Pitfalls** · 72% of the patents have all inventors in one region, so the region
projection sees no tie for them — keep within-region collaboration as a separate
variable. One application has 68 inventors and alone generates 2,278 ties.
`person_id` is algorithmic: a split identity breaks paths, a merged one invents
hubs.

### `pat_all_inventors_ITgreen.csv.gz` — 21,916 rows · 7,914 patents · 13,549 inventors

Same columns plus `green` (1/0: does this application carry a Y02/Y04S code).
All EPO applications 2010–2019 of the inventors above, green or not.

**What you can build** · the same network under a **wider boundary** (all
co-patenting, not only green) → how much of an inventor's position is green;
green vs non-green subnetworks; inventors who bridge the two (recombination).

### `pat_green_tech_IT.csv.gz` — 17,644 rows · 488 CPC4 classes

| variable | type | description | notes |
|---|---|---|---|
| `appln_id` | integer | patent application | joins the two files above |
| `prio_year` | integer | priority year | |
| `cpc4` | character | CPC 4-digit class, e.g. `Y02E`, `H01M` | one row per patent × class |
| `tech_gt` | integer | 1 if this class is itself a green class (Y02/Y04S) | a green patent also carries non-green classes |

**What you can build** · patent-level **technology co-classification** →
knowledge space; technological variety of a patent or a portfolio; green vs
non-green recombination.

### `green_region_tech_EU.csv.gz` — 120,120 rows · 1,329 regions · 631 classes

| variable | type | description | notes |
|---|---|---|---|
| `reg_code` | character | NUTS-3 / TL3 region | aggregate to NUTS-2 with `substr(., 1, 4)` |
| `ctry_code` | character | country, ISO-2 | EU27 + CH, NO, UK |
| `cpc4` | character | CPC 4-digit class | |
| `period` | character | `2010-2014` or `2015-2019` | two windows → entry/exit analysis |
| `n_pat` | numeric | **fractional** patent count | a patent is split across its regions and classes, so the world total is the patent count |

**What you can build** · region × technology matrix → RTA/Balassa → diversity and
ubiquity → relatedness, relatedness density, variety (related/unrelated),
coherence, complexity (KCI/TCI); **diversification**: entry into new
technologies between the two periods.

**Pitfalls** · fractional counts are the reason totals are not integers; never
mix them with full counts in the same regression. At NUTS-3 many cells are tiny —
aggregate before computing RTA.

### `green_tech_cooc_EU.csv.gz` (13,212 rows) and `green_tech_npat_EU.csv` (612 rows)

`cpc4_i`, `cpc4_j`, `n_cooc`: how often two classes appear on the same EU green
patent, 2015–2019 (upper triangle only, `n_cooc >= 2`). `cpc4`, `n_pat`: patents
per class, the denominator for association strength.

**What you can build** · the **knowledge space** (co-classification proximity),
technological relatedness matrices, betweenness of technologies (candidate
general-purpose technologies), communities of the green knowledge base.

### `cpc4_def.csv` (663) and `green_coinvention_country_benchmark.csv` (17)

Class labels; and network statistics of the green co-invention network for 17
countries, 2015–2019 (`patents`, `inventors`, `ties`, `avg_team`, `avg_degree`,
`density`, `giant_share`, `clustering`) — a ready-made comparison for putting one
country in perspective.

---

## 2. EU-funded projects — CORDIS, Horizon Europe (release 2026-08-06)

Unit: **one organisation in one project**. Subset: 4,604 projects flagged 100%
climate-relevant or classified under a sustainability euroSciVoc term.
Here the tie is **contractual, dated and priced** — and the population includes
the actors patents never see (universities, public bodies, NGOs, SMEs).

### `cordis_he_participants.csv.gz` — 44,644 rows · 16,116 organisations

| variable | type | description | notes |
|---|---|---|---|
| `project_id` | integer | project (the **event**) | joins `cordis_he_projects.csv` |
| `org_id` | integer | organisation id (the **actor**) | EC participant identification code |
| `org_name` | character | legal name | raw CORDIS: spacing errors, aggregated legal entities (CNRS, Fraunhofer) |
| `short_name` | character | acronym | often empty |
| `country` | character | ISO-2, 152 countries | includes non-EU partners |
| `nuts` | character | NUTS code | 5 characters = NUTS-3 (94% of rows), 2 = country only (non-EU), a few malformed |
| `activity_type` | character | `HES` higher education, `REC` research org, `PRC` private, `PUB` public body, `OTH` other | the single most useful attribute |
| `sme` | logical | SME flag | only meaningful for `PRC` |
| `role` | character | `coordinator`, `participant`, `thirdParty`, `associatedPartner` | coordination = a directed relation you can exploit |
| `ec_contrib_org` | numeric | EC contribution to this organisation, EUR | money on the tie, rare in network data |
| `start_year` | integer | project start, 2021–2027 | tie formation over time |
| `acronym` | character | project acronym | |
| `climate_100` | integer | 1 if the project carries a 100% climate policy marker | EU Rio-marker style tagging: 0 / 40 / 100 |

**What you can build**
- organisation collaboration network → hubs, brokers, communities, assortativity by country and by type;
- **region** (NUTS-2) and **country** networks from the same three lines;
- coordinator → participant **directed** network (who leads whom);
- weighted-by-money networks (`ec_contrib_org`), funding concentration;
- tie dynamics: new vs repeated partners across `start_year` windows;
- public–private and academia–industry mixing matrices.

**Pitfalls** · consortium size ranges from 1 to 96: mega-consortia turn into
cliques (use `max_size`). Aggregating to NUTS drops every non-European partner
(6.5% of participations). `org_id` mixes headquarters and subsidiaries: CNRS is
one node.

### `cordis_he_projects.csv` — 4,604 rows

| variable | type | description | notes |
|---|---|---|---|
| `project_id`, `acronym`, `title` | | project identity | |
| `status` | character | `SIGNED`, `CLOSED`, `TERMINATED` | terminated projects still have participants |
| `start_date`, `end_date` | date | contractual dates | |
| `ec_contrib`, `total_cost` | numeric | EU contribution and total budget, EUR | decimal comma in the raw file |
| `funding_scheme` | character | 27 instruments (`HORIZON-RIA`, `HORIZON-ERC`, `HORIZON-IA`, …) | instrument shapes consortium size |
| `topic` | character | call topic code | fine-grained thematic grouping |
| `legal_basis` | character | programme part (`HORIZON.2.5`, …) | **34 of 4,604 rows are mangled** by embedded newlines in the source CSV; unused by the scripts — use `programme` |
| `programme` | character | 15 readable programme names | the clean version of `legal_basis` |
| `start_year` | integer | | |
| `climate_100`, `sust_scivoc` | integer | selection flags: climate marker / sustainability euroSciVoc term | |

### `cordis_he_scivoc.csv.gz` — 11,785 rows · 686 fields

`project_id`, `sci_voc` (field name), `sci_voc_path` (full taxonomy path, e.g.
`natural sciences/chemical sciences/electrochemistry`).

**What you can build** · the **thematic space** of EU climate research (field
co-occurrence within projects), interdisciplinarity of a project or of an
organisation's portfolio, thematic communities, and — with the taxonomy path —
distance between fields at different levels of the hierarchy.

---

## 3. Publications — OpenAlex (open API, CC0)

Unit: **one author × institution × work**. The extract caches 2,000 articles
matching "circular economy" in title or abstract, 2020–2024, with ORCID.
The live API is the point: `03_publications.R` shows aggregate queries
(`group_by`) and cursor paging.

### `openalex_ce_authorships.csv.gz` — 10,517 rows · 1,999 works · 6,639 authors · 2,223 institutions

| variable | type | description | notes |
|---|---|---|---|
| `work_id` | character | OpenAlex work id (the **event**) | `W…` |
| `year` | integer | publication year | |
| `cited_by` | integer | citation count at extraction | a stock, not a flow: never comparable across years without normalisation |
| `author_id` | character | OpenAlex author id (the **actor**) | `A…`, algorithmically disambiguated |
| `author_name` | character | display name | |
| `position` | character | `first`, `middle`, `last` | proxies for contribution and seniority |
| `inst_id` | character | institution id | `I…`; a work × author can have several |
| `inst_name` | character | institution name | |
| `inst_country` | character | ISO-2 of the institution | |
| `inst_type` | character | `education`, `company`, `healthcare`, `government`, … | 10 values |

**What you can build** · co-authorship networks of **authors**,
**institutions** and **countries**; academia–industry collaboration (via
`inst_type`); first/last-author networks; and, with the API, citation networks,
topic distributions and time series without downloading records.

**Pitfalls** · a keyword query *is* a sample definition — report it. Author
disambiguation errors concentrate on common names. Institution ids mix
universities with their hospitals.

---

## 4. Trade in value added — OECD TiVA (FDVA table, SDMX API)

Unit: **value added of a source economy embodied in the final demand of a
destination economy**, million USD. 80 economies. This is the second family of
network data: the tie is **observed, directed and valued**, the graph is
complete, and TiVA is the output of an inter-country input–output model — the
network inherits its assumptions.

### `tiva_va_bilateral.csv.gz` — 25,600 rows (4 years × 80 × 80)

| variable | type | description | notes |
|---|---|---|---|
| `year` | integer | 1995, 2005, 2015, 2022 | benchmark years for structural change |
| `source` | character | ISO-3 of the economy **creating** the value added | |
| `destination` | character | ISO-3 of the economy whose **final demand** absorbs it | |
| `va_musd` | numeric | value added, million USD, current prices | current prices: deflate before comparing years in levels |

**What you can build** · directed weighted GVC network → in/out strength (VA
absorbed / supplied), net positions, dyadic imbalance, export-market
concentration (HHI); **filtering** (threshold, top-*k*, share, disparity filter)
and backbones; trade blocs by community detection, with and without
size normalisation; Gould–Fernandez brokerage between world regions; structural
change 1995 → 2022.

**Pitfalls** · the diagonal (`source == destination`) is domestic value added,
not a tie — drop it. Density is 1, so degree-based measures are meaningless.
Betweenness treats weights as distances: pass `1/w` or `weights = NA`.

### `tiva_va_by_industry.csv.gz` — 51,200 rows · 8 industries, 2022

Same columns plus `industry` (ISIC rev.4 aggregates, non-overlapping:
`A`, `C10T12`, `C24_25`, `C26_27`, `C29_30`, `D_E`, `J`, `M_N`; labels in
`tiva_industry_labels.csv`).

**What you can build** · one network per industry → how regionalised each value
chain is (modularity, NMI with geography), industry-specific hubs, comparison of
manufacturing against services.

### `country_regions.csv` — 80 rows

`iso3`, `region` (10 macro-regions). The group partition used for brokerage roles
and for testing whether communities are geographic.

---

## 5. Trade in goods — CEPII BACI (HS22, version V202501)

Unit: **exporter × importer × product × year**, thousand USD in the source,
converted to million USD here. Gross flows: what crosses the border, before any
value-added correction.

### `baci_bilateral_2023.csv.gz` — 16,360 pairs · 221 exporters

| variable | type | description | notes |
|---|---|---|---|
| `exporter`, `importer` | character | ISO-3 | special BACI codes (`S19`, …) removed |
| `exports_musd` | numeric | gross exports, million USD, 2023, all products | pairs below 1 MUSD dropped |

**What you can build** · the world trade web (directed, valued); the comparison
with TiVA that reveals **re-export hubs and assembly platforms** (value added per
dollar exported: 0.35 for Slovakia, 1.63 for Ireland); gravity-style normalised
intensities.

### `baci_country_product_2023.csv.gz` — 62,822 rows · 126 countries · 933 products

| variable | type | description | notes |
|---|---|---|---|
| `exporter` | character | ISO-3 | only exporters above 5 bn USD total: RCA on a tiny exporter is noise |
| `hs4` | character | HS 4-digit product, e.g. `0101`, `8703` | **read as character**: `fread` turns `"0101"` into `101` |
| `exports_musd` | numeric | exports of that product, million USD | products below 1 bn USD world exports dropped |

**What you can build** · country × product matrix → RCA → **product space**
(proximity), relatedness density, **ECI and PCI** (economic and product
complexity); export diversification and structural change; the same machinery as
the patent knowledge space, on a different category system.

### `hs4_labels.csv` — 933 rows

`hs4`, `label` (short heading), `hs2` (HS chapter). `hs2` also needs
`colClasses = "character"`: chapters 01–09 lose their leading zero otherwise.

---

## Cross-source summary

| source | tie | direction | value on the tie | levels available | best for |
|---|---|---|---|---|---|
| REGPAT | co-invention (inferred) | undirected | n. of shared patents | inventor, region, technology | knowledge production, recombination |
| CORDIS | co-participation (contract) | undirected (directed via `role`) | shared projects, EUR | organisation, region, country, topic | policy-made collaboration, funding |
| OpenAlex | co-authorship (inferred) | undirected | shared works | author, institution, country | science, science–industry links |
| TiVA | value-added flow (observed) | **directed** | million USD | country, country × industry | GVC position, structural change |
| BACI | gross export (observed) | **directed** | million USD | country, country × product | specialisation, complexity |

**One rule across all five**: the level of aggregation is not a technical
detail. An inventor network averaged by region and a region network built
directly from the same patents correlate at −0.02 in this data. They are
different objects answering different questions — pick the one your theory is
about, and say so.
