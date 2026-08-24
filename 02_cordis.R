#------------------------------------------------------------------------------#
#  DAISY 2026 - Hands-on: Network analysis with innovation data
#  BLOCK 2 (~20 min) - CORDIS: networks of EU-funded collaborative projects
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

## (b) E-I index: share of ties that cross national borders (compare with the
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
## 5. Aggregate the same data at country level (and plot it)
## ===========================================================================
## The projection helper works at any level of aggregation: just change "actor".
pr_ctry <- proj_two_mode(part, "project_id", "country")
gc_ctry <- make_net(pr_ctry)
gc_ctry <- delete_edges(gc_ctry, E(gc_ctry)[weight < 30])   # keep readable
gc_ctry <- induced_subgraph(gc_ctry, V(gc_ctry)[degree(gc_ctry) > 0])

V(gc_ctry)$projects <- V(gc_ctry)$n_events
ggraph(gc_ctry, layout = "stress") +
  geom_edge_link0(aes(edge_width = weight), edge_colour = "grey80", edge_alpha = .8) +
  scale_edge_width(range = c(0.1, 3)) +
  geom_node_point(aes(size = projects), fill = "#2c7fb8", shape = 21, colour = "white") +
  geom_node_text(aes(label = name), size = 3, repel = TRUE) +
  scale_size(range = c(2, 12)) +
  theme_graph(base_family = "sans") + theme(legend.position = "none") +
  labs(title = "Country co-participation, Horizon Europe climate projects")

## Country-level indicators: raw weights favour big countries, so normalise.
## A simple revealed-collaboration index: observed ties / expected under
## independence given each country's number of participations.
el <- as.data.table(as_data_frame(gc_ctry, what = "edges"))
tot <- data.table(country = V(gc_ctry)$name, part = V(gc_ctry)$n_events)
el <- merge(merge(el, tot, by.x = "from", by.y = "country"),
            tot, by.x = "to", by.y = "country", suffixes = c("_f", "_t"))
el[, rci := weight / (part_f * part_t / sum(tot$part))]
el[order(-rci)][1:10, .(from, to, weight, rci = round(rci, 2))]
## Small neighbouring countries collaborate far above expectation: geography and
## institutional proximity survive even inside a supranational programme.

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
## 7. IF WE HAVE TIME - the same data as a *topic* network
## ===========================================================================
## euroSciVoc classifies each project into scientific fields: projecting the
## project x field matrix gives the thematic space of EU climate research -
## exactly the logic we use for patent technology classes in block 4.
sv <- fread(daisy_data("cordis_he_scivoc.csv.gz"))
pr_topic <- proj_two_mode(sv, "project_id", "sci_voc")
gt_net <- make_net(pr_topic)
gt_net <- delete_edges(gt_net, E(gt_net)[weight < 10])
gt_net <- induced_subgraph(gt_net, V(gt_net)[degree(gt_net) > 0])
sort(degree(gt_net), decreasing = TRUE)[1:15]

ggraph(gt_net, layout = "stress") +
  geom_edge_link0(aes(edge_width = weight), edge_colour = "grey85") +
  scale_edge_width(range = c(0.1, 2)) +
  geom_node_point(aes(size = n_events), fill = "#41ab5d", shape = 21, colour = "white") +
  geom_node_text(aes(label = name), size = 2.6, repel = TRUE, max.overlaps = 20) +
  scale_size(range = c(1, 9)) +
  theme_graph(base_family = "sans") + theme(legend.position = "none") +
  labs(title = "Thematic space of Horizon Europe climate projects (euroSciVoc)")
