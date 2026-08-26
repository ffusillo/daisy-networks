#------------------------------------------------------------------------------#
#  DAISY 2026 - Hands-on: Network analysis with innovation data
#  BLOCK 2 of the session (~16 min) - CORDIS: EU-funded collaborative projects
#
#  Data: CORDIS Horizon Europe, open data (CC-BY), release 2026-08-06,
#        https://cordis.europa.eu/dataset  (project.csv, organization.csv,
#        policyPriorities.csv, euroSciVoc.csv ...)
#  Extract used here: the 4,604 HE projects tagged as 100% climate-relevant or
#  classified under a sustainability-related euroSciVoc term, and their 44,644
#  participations (16,116 distinct organisations).
#
#  Why this source: unlike patents and papers, here the collaboration is a
#  *contractual, funded and dated* relationship, with money attached to each
#  partner - and it covers the actors that patents miss (universities, public
#  bodies, NGOs, SMEs).
#------------------------------------------------------------------------------#
source("00_setup.R")

## ===========================================================================
## 1. The data
## ===========================================================================
part <- fread(daisy_data("cordis_he_participants.csv.gz"))
proj <- fread(daisy_data("cordis_he_projects.csv"))

part          # one row = one organisation in one project
proj[1:3, .(project_id, acronym, start_year, programme, ec_contrib)]

part[, .(projects = uniqueN(project_id), organisations = uniqueN(org_id),
         participations = .N)]

## Consortium size: the "event size" that will drive the projection
size <- part[, .(partners = .N), by = project_id]
size[, .(mean = mean(partners), median = as.double(median(partners)), max = max(partners))]
size[order(-partners)][1:5]
proj[project_id %in% size[order(-partners)][1:3]$project_id, .(acronym, title)]

## Who participates? (HES = higher education, REC = research org, PRC = private,
## PUB = public body, OTH = other)
part[, .(participations = .N, orgs = uniqueN(org_id),
         ec_meur = round(sum(ec_contrib_org, na.rm = TRUE) / 1e6)),
     by = activity_type][order(-participations)]

## Top countries by participation and by money
part[, .(participations = .N,
         ec_meur = round(sum(ec_contrib_org, na.rm = TRUE) / 1e6)),
     by = country][order(-ec_meur)][1:15]

## Projects per year
part[, .(projects = uniqueN(project_id)), by = start_year][order(start_year)]

## ===========================================================================
## 2. The organisation collaboration network
## ===========================================================================
## Same helper as for patents: the event is the project, the actor the partner.
## Mega-consortia (60+ partners) would dominate the topology, so we look at both.
pr_all <- proj_two_mode(part, "project_id", "org_id")
pr_lim <- proj_two_mode(part, "project_id", "org_id", max_size = 40)
c(ties_all = nrow(pr_all$edges), ties_below_40_partners = nrow(pr_lim$edges))

org_attr <- part[, .(org_name = org_name[1], country = country[1],
                     nuts = nuts[1], type = activity_type[1], sme = sme[1],
                     eur = sum(ec_contrib_org, na.rm = TRUE)), by = org_id]

g <- make_net(pr_lim, node_attr = org_attr, by = "org_id")
g

## Anatomy: EU funding networks look nothing like co-invention networks
cmp <- components(g)
c(nodes = vcount(g), edges = ecount(g), density = edge_density(g),
  components = cmp$no, giant_share = max(cmp$csize) / vcount(g),
  clustering = transitivity(g, type = "global"),
  avg_degree = mean(degree(g)))

## a connected, dense, high-clustering core: the "policy-made" network

## ===========================================================================
## 3. Who is central - and does centrality mean the same as money?
## ===========================================================================
V(g)$degree   <- degree(g)                       # distinct partners
V(g)$strength <- strength(g)                     # partner-projects
V(g)$betw     <- betweenness(g, weights = NA, normalized = TRUE)
V(g)$eigen    <- eigen_centrality(g, weights = NA)$vector
V(g)$constr   <- constraint(g)

nodes <- as.data.table(as_data_frame(g, what = "vertices"))
nodes[order(-degree)][1:15, .(org_name, country, type, n_events, degree, betw,
                              eur_meur = round(eur / 1e6, 1))]

## Money and network position are related but not the same thing
nodes[, round(cor(cbind(n_events, degree, strength, betw, eigen, eur),
                  use = "pairwise"), 2)]

## Brokers vs hubs: rank difference tells you who bridges rather than accumulates
nodes[, `:=`(r_deg = frankv(-degree), r_betw = frankv(-betw))]
nodes[degree > 20][order(r_betw - r_deg)][1:10,
      .(org_name, country, type, degree, betw = round(betw, 4))]

## ===========================================================================
## 4. Is EU research integrated, or nationally clustered?
## ===========================================================================
## (a) assortativity: do organisations partner with similar organisations?
assortativity_nominal(g, factor(V(g)$country))     # by country
assortativity_nominal(g, factor(V(g)$type))        # by type of organisation
assortativity_degree(g)                            # hubs with hubs?

## (b) Share of ties that cross national borders (compare with the
##     0.26-0.46 range we found for regions in the co-invention network)
el0 <- as.data.table(as_data_frame(g, what = "edges"))
ctry_of <- setNames(V(g)$country, V(g)$name)
el0[, cross := as.integer(ctry_of[from] != ctry_of[to])]
el0[, .(cross_border_share = weighted.mean(cross, weight))]

## (c) communities in the giant component and their national composition
##     (which algorithm? which resolution? how stable? -> 06_brokerage_communities.R,
##      which also computes Gould-Fernandez brokerage roles on this same network)
gg <- giant(g)
comm <- cluster_louvain(gg, weights = E(gg)$weight)
length(comm); sort(sizes(comm), decreasing = TRUE)[1:8]

memb <- data.table(org_id = V(gg)$name, country = V(gg)$country,
                   type = V(gg)$type, comm = membership(comm))
top_comm <- memb[, .N, by = comm][order(-N)][1:6]$comm

## how concentrated is each community by country? (HHI = 1 -> single country)
memb[comm %in% top_comm, .(
  orgs = .N,
  top_country = names(which.max(table(country))),
  top_share = round(max(prop.table(table(country))), 2),
  hhi = round(sum(prop.table(table(country))^2), 3)), by = comm][order(-orgs)]

## benchmark: concentration of the whole network
memb[, .(hhi_all = round(sum(prop.table(table(country))^2), 3))]
## => communities are thematic-institutional, not national: the opposite of what
##    we found for co-invention. Worth a slide in any paper on EU integration.

## ===========================================================================
## 5. Aggregate the same data at NUTS-2 level (and plot it)
## ===========================================================================
## The projection helper works at any level of aggregation: just change "actor".
## Regions are the level at which most of the innovation-policy literature works,
## and CORDIS geocodes every participant to a NUTS code.

## First look at what the geography column actually contains
part[, .N, by = .(nuts_length = nchar(nuts))][order(-N)]
## 5 characters = NUTS-3, 2 = country only (non-EU partners), a handful of odd
## ones. NUTS exists only for Europe: aggregating to regions silently DROPS
## every third-country partner. That is a change of population, not a detail.
part[, nuts2 := ifelse(nchar(nuts) >= 4, substr(nuts, 1, 4), NA_character_)]
part[, .(participations = .N, with_region = sum(!is.na(nuts2)),
         share_kept = round(mean(!is.na(nuts2)), 3))]

pr_reg  <- proj_two_mode(part[!is.na(nuts2)], "project_id", "nuts2")
reg_att <- unique(part[!is.na(nuts2), .(nuts2, country)], by = "nuts2")
g_reg   <- make_net(pr_reg, node_attr = reg_att, by = "nuts2")
g_reg

c(regions = vcount(g_reg), ties = ecount(g_reg),
  density = round(edge_density(g_reg), 3),
  giant_share = round(max(components(g_reg)$csize) / vcount(g_reg), 3))

## Which regions sit at the centre of EU climate research?
V(g_reg)$degree   <- degree(g_reg)
V(g_reg)$strength <- strength(g_reg)
V(g_reg)$betw     <- betweenness(g_reg, weights = NA, normalized = TRUE)
reg_nodes <- as.data.table(as_data_frame(g_reg, what = "vertices"))
setnames(reg_nodes, c("name", "n_events"), c("nuts2", "projects"))
reg_nodes[order(-strength)][1:12, .(nuts2, country, projects, degree, strength,
                                    betw = round(betw, 3))]

## Is regional collaboration national or European? (compare with the 26-46%
## extra-regional share of the co-invention network in block 1)
el_r <- as.data.table(as_data_frame(g_reg, what = "edges"))
ctry_of_reg <- setNames(V(g_reg)$country, V(g_reg)$name)
el_r[, cross := as.integer(ctry_of_reg[from] != ctry_of_reg[to])]
el_r[, .(cross_country_share = round(weighted.mean(cross, weight), 3))]

## Region-level indicators: raw weights favour big regions, so normalise.
## A simple revealed-collaboration index: observed ties / expected under
## independence given each region's number of participations.
tot <- data.table(nuts2 = V(g_reg)$name, part = V(g_reg)$n_events)
el_r <- merge(merge(el_r, tot, by.x = "from", by.y = "nuts2"),
              tot, by.x = "to", by.y = "nuts2", suffixes = c("_f", "_t"))
el_r[, rci := weight / (part_f * part_t / sum(tot$part))]
el_r[weight >= 10][order(-rci)][1:10, .(from, to, weight, rci = round(rci, 1),
                                        cross)]
## The strongest *relative* ties are pairs of regions in the same country, or in
## neighbouring ones: geography and institutional proximity survive even inside
## a supranational programme designed to overcome them.

## Draw the backbone: 323 regions are too many, keep the strongest ties
g_plot <- delete_edges(g_reg, E(g_reg)[weight < 40])
g_plot <- induced_subgraph(g_plot, V(g_plot)[degree(g_plot) > 0])
V(g_plot)$projects <- V(g_plot)$n_events

ggraph(g_plot, layout = "stress") +
  geom_edge_link0(aes(edge_width = weight), edge_colour = "grey80", edge_alpha = .8) +
  scale_edge_width(range = c(0.1, 2.5)) +
  geom_node_point(aes(size = projects, fill = country), shape = 21, colour = "white") +
  geom_node_text(aes(label = name), size = 2.6, repel = TRUE, max.overlaps = 20) +
  scale_size(range = c(2, 11)) +
  theme_graph(base_family = "sans") + theme(legend.position = "none") +
  labs(title = "NUTS-2 co-participation, Horizon Europe climate projects",
       subtitle = "ties with at least 40 shared projects; node size = participations")

## Exercise for later: the same three lines with actor = "country" give the
## country network. Compare the two rankings - Ile-de-France against France.

## ===========================================================================
## 6. IF WE HAVE TIME - tie formation: new or repeated partners?
## ===========================================================================
## A dynamic indicator you can build from any project database: how much of the
## collaboration in t is with partners already met before t?
early <- proj_two_mode(part[start_year <= 2022], "project_id", "org_id")$edges
late  <- proj_two_mode(part[start_year >= 2024], "project_id", "org_id")$edges
key   <- function(d) paste(pmin(d$from, d$to), pmax(d$from, d$to))
mean(key(late) %in% key(early))       # share of repeated ties
## Repetition rate = trust/lock-in vs renewal of the consortium ecosystem.

## ===========================================================================
## 7. And the same data as a *topic* network
## ===========================================================================
## euroSciVoc classifies every project into scientific fields, so the SAME
## projection turns projects into a map of what EU climate research is about.
## We build that map in 04_indicators.R, next to the patent knowledge space and
## the product space: it is the same construction applied to three different
## category systems (technologies, topics, products).
