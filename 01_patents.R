#------------------------------------------------------------------------------#
#  DAISY 2026 - Hands-on: Network analysis with innovation data
#  BLOCK 1 (~22 min) - PATENT DATA: the co-inventor network
#
#  Data: OECD REGPAT (May 2025), EPO applications, inventor-region file, joined
#  with CPC codes. Extract used here: all EPO patent applications with priority
#  year 2010-2019 that (i) have at least one Italian inventor and (ii) carry at
#  least one "green" CPC code (Y02/Y04S climate-change mitigation tagging).
#
#  What we do: from a raw patent-inventor table to (a) a co-invention network,
#  (b) inventor-level positional measures, (c) region-level indicators you can
#  put in a regression.
#------------------------------------------------------------------------------#
source("00_setup.R")

## ===========================================================================
## 1. The raw data: one row = one inventor on one patent
## ===========================================================================
inv <- fread(daisy_data("pat_green_inventors_IT.csv.gz"))
inv
str(inv)

## appln_id  : patent application (the "event")
## person_id : REGPAT disambiguated inventor id (the "actor")
## reg_code  : NUTS-3 region of residence of the inventor, ctry_code: country
## prio_year : priority year = closest to the moment of invention

## How much data do we have?
inv[, .(patents = uniqueN(appln_id), inventors = uniqueN(person_id))]

## Patents and inventors per year
by_year <- inv[, .(patents = uniqueN(appln_id), inventors = uniqueN(person_id)),
               by = prio_year][order(prio_year)]
by_year

ggplot(by_year, aes(prio_year, patents)) +
  geom_col(fill = "grey40") +
  labs(x = NULL, y = "green EPO applications with IT inventors") +
  theme_minimal()

## Team size: this is what drives the density of the co-invention network
team <- inv[, .(size = uniqueN(person_id)), by = appln_id]
team[, .(mean = mean(size), median = as.double(median(size)), max = max(size))]
table(team$size)

## Careful: many patents have ONE inventor -> they will be isolated nodes
## Careful #2: inventors can be counted in several regions (reg_share) and
## patents in several technologies; use fractional counts when you aggregate.

## Where are the inventors? (NUTS-3 -> NUTS-2 by truncation)
inv[, nuts2 := substr(reg_code, 1, 4)]
inv[ctry_code == "IT", .(inventors = uniqueN(person_id)), by = nuts2][order(-inventors)][1:10]

## ===========================================================================
## 2. From the two-mode (patent x inventor) to the co-invention network
## ===========================================================================
## The tie is *not observed*: we infer it from co-participation in a patent.
pr <- proj_two_mode(inv, event = "appln_id", actor = "person_id")

head(pr$edges)          # weight = number of patents co-invented
head(pr$nodes)          # n_events = number of patents of the inventor

## WARNING - the projection turns every team into a clique. The biggest patent
## in this sample has 68 inventors: that single document produces 68*67/2 =
## 2,278 ties. Whether to keep such documents is a research design choice.
team[order(-size)][1:5]
pr20 <- proj_two_mode(inv, "appln_id", "person_id", max_size = 20)
c(ties_all = nrow(pr$edges), ties_teams_below_20 = nrow(pr20$edges))

## Attach inventor attributes (one row per inventor: first region observed)
attr_inv <- inv[, .(inv_name = inv_name[1], ctry = ctry_code[1],
                    nuts2 = nuts2[1], reg = reg_code[1],
                    first_year = min(prio_year)), by = person_id]

g <- make_net(pr, node_attr = attr_inv, by = "person_id")
g

## Basic anatomy of the network
vcount(g); ecount(g)
edge_density(g)
mean(degree(g))
table(degree(g) == 0)                 # isolates = single-inventor patents only

## Components: co-invention networks are always highly fragmented
cmp <- components(g)
cmp$no                                # number of components
sort(cmp$csize, decreasing = TRUE)[1:10]
max(cmp$csize) / vcount(g)            # share of inventors in the giant component

gc_net <- giant(g)
gc_net

## Small-world-ness of the giant component
mean_distance(gc_net)
diameter(gc_net, weights = NA)
transitivity(gc_net, type = "global")   # very high: teams are cliques by construction
## benchmark against a random graph of the same size and density
rnd <- sample_gnm(vcount(gc_net), ecount(gc_net))
c(observed = transitivity(gc_net, type = "global"),
  random   = transitivity(rnd, type = "global"))

## Degree distribution: fat tailed, as usual in collaboration networks
ggplot(data.table(k = degree(g)), aes(k)) +
  geom_histogram(binwidth = 1, fill = "grey40") +
  scale_y_log10() + labs(x = "degree (number of distinct co-inventors)") +
  theme_minimal()

## ===========================================================================
## 3. Positions: who matters, and in which sense?
## ===========================================================================
V(g)$degree      <- degree(g)
V(g)$strength    <- strength(g)                       # weighted by n. of patents
V(g)$betw        <- betweenness(g, weights = NA, normalized = TRUE)
V(g)$eigen       <- eigen_centrality(g, weights = NA)$vector
V(g)$constraint  <- constraint(g)                     # Burt: LOW = structural holes
V(g)$clust       <- transitivity(g, type = "local", isolates = "zero")

cent <- as.data.table(as_data_frame(g, what = "vertices"))
setnames(cent, "name", "person_id")

## Top inventors by different criteria - they are NOT the same people
cent[order(-degree)][1:10, .(inv_name, nuts2, n_events, degree, strength, betw)]
cent[order(-betw)][1:10,   .(inv_name, nuts2, n_events, degree, betw, constraint)]

## Are patent counts and network position the same information?
cent[n_events > 0, round(cor(.SD, use = "pairwise"), 2),
     .SDcols = c("n_events", "degree", "strength", "betw", "eigen", "constraint")]

## => centrality is correlated with productivity but far from collinear: this is
##    why network position enters innovation regressions on its own.

## ===========================================================================
## 4. Visualise (only ever plot a subgraph you can actually read)
## ===========================================================================
sub <- giant(g)
sub <- induced_subgraph(sub, V(sub)[degree(sub) > 1])

ggraph(sub, layout = "stress") +
  geom_edge_link0(aes(edge_width = weight), edge_colour = "grey75") +
  scale_edge_width(range = c(0.2, 1.5)) +
  geom_node_point(aes(size = degree, fill = ctry), shape = 21, colour = "white") +
  scale_size(range = c(1, 6)) +
  theme_graph(base_family = "sans") +
  labs(title = "Giant component, green co-invention network (IT, 2010-2019)",
       fill = "inventor country")

## ===========================================================================
## 5. From nodes to variables: region-level network indicators
## ===========================================================================
## This is what usually ends up in an econometric model: aggregate inventor
## positions by region (or firm, or year window) and use them as regressors.
reg_ind <- cent[ctry == "IT" & nuts2 != "", .(
  inventors        = .N,
  patents          = sum(n_events),
  avg_degree       = mean(degree),
  avg_strength     = mean(strength),
  avg_betweenness  = mean(betw),
  avg_constraint   = mean(constraint, na.rm = TRUE),
  share_connected  = mean(degree > 0),
  top_inventor     = inv_name[which.max(degree)]
), by = nuts2][order(-patents)]
reg_ind[1:15]

## Cross-border openness: share of an inventor's ties that go outside the region
el <- as.data.table(as_data_frame(g, what = "edges"))
nuts_of <- setNames(V(g)$nuts2, V(g)$name)
el[, `:=`(n_from = nuts_of[from], n_to = nuts_of[to])]
ext <- rbind(el[, .(nuts2 = n_from, ext = as.integer(n_from != n_to), weight)],
             el[, .(nuts2 = n_to,   ext = as.integer(n_from != n_to), weight)])
open_reg <- ext[, .(external_tie_share = weighted.mean(ext, weight)), by = nuts2]
reg_ind <- merge(reg_ind, open_reg, by = "nuts2", all.x = TRUE)
reg_ind[order(-patents)][1:10, .(nuts2, patents, avg_degree, avg_constraint,
                                 share_connected, external_tie_share)]

fwrite(reg_ind, "output_region_network_indicators.csv")

## ===========================================================================
## 6. The network boundary is a decision, not a fact
## ===========================================================================
## A small function that runs the whole pipeline and returns summary statistics
net_stats <- function(dt) {
  p <- proj_two_mode(dt, "appln_id", "person_id")
  gg <- make_net(p)
  cmp <- components(gg)
  data.table(inventors = vcount(gg), ties = ecount(gg),
             density = edge_density(gg),
             avg_degree = mean(degree(gg)),
             giant_share = max(cmp$csize) / vcount(gg),
             clustering = transitivity(gg, type = "global"))
}

## So far a tie exists only if two inventors share a GREEN patent. But the same
## inventors also collaborate on non-green patents. Same actors, wider boundary:
inv_all <- fread(daisy_data("pat_all_inventors_ITgreen.csv.gz"))
inv_all[, .(patents = uniqueN(appln_id), green = uniqueN(appln_id[green == 1]))]

boundary <- rbind(
  `green ties only` = net_stats(inv),
  `all co-patenting ties` = net_stats(inv_all), idcol = "boundary")
boundary
## Connectivity, average degree and the giant component all move: any statement
## about "the" position of an inventor is conditional on this choice.

## ... and so is the *population* boundary. Same pipeline, run for you on the
## full REGPAT for 17 countries (green patents, 2015-2019):
bench <- fread(daisy_data("green_coinvention_country_benchmark.csv"))
bench[order(-patents)]
## Note (i) how small the giant component is everywhere over a 5-year window,
## (ii) FI: avg_degree of 21 driven by a handful of very large teams - always
## look for the mega-document before interpreting a "dense" network.

## ===========================================================================
## 7. IF WE HAVE TIME - does the network change over time?
## ===========================================================================
## We reuse net_stats() defined above on moving windows.
windows <- list(`2010-2013` = 2010:2013, `2014-2016` = 2014:2016,
                `2017-2019` = 2017:2019)
evo <- rbindlist(lapply(windows, function(y) net_stats(inv[prio_year %in% y])),
                 idcol = "window")
evo

## Discussion: fragmentation, densification, and the sensitivity of ALL of this
## to the length of the time window - a modelling choice, not a data property.
