#------------------------------------------------------------------------------#
#  DAISY 2026 - Hands-on: Network analysis with innovation data
#  BLOCK 6 - BROKERAGE AND COMMUNITY DETECTION: a toolbox with hints
#  (reference material: the ideas are used live in blocks 2 and 4, the menu of
#   algorithms and benchmarks is here for when you write your own paper)
#
#  Two questions come back in every seminar: "who is the broker?" and "how do I
#  find groups?". Both have several answers, and the answers do not agree. This
#  script is the reference: what each measure actually counts, how to compute it
#  when igraph has no function for it, and how to report it so that a referee
#  cannot ask "what would happen with another algorithm?".
#
#  It runs on the networks built in blocks 2 and 5, but every function here works
#  on any igraph object.
#------------------------------------------------------------------------------#
source("00_setup.R")

## ===========================================================================
## PART A - BROKERAGE
## ===========================================================================
## Four different ideas travel under the same word:
##
##   1. BETWEENNESS         - how often you sit on shortest paths (global flow)
##   2. STRUCTURAL HOLES    - Burt: are your contacts disconnected from each
##                            other? (constraint, effective size, efficiency)
##   3. GOULD-FERNANDEZ     - what KIND of gap do you bridge, given a group
##      brokerage roles       partition (coordinator, gatekeeper, representative,
##                            consultant, liaison)
##   4. E-I INDEX           - how outward-looking are your ties (group level)
##
## They answer different questions. Pick the one that matches your theory, and
## say why. A paper that reports "betweenness" when the argument is about
## bridging two communities is answering the wrong question.

part <- fread(daisy_data("cordis_he_participants.csv.gz"))
g_org <- make_net(proj_two_mode(part, "project_id", "org_id", max_size = 40),
                  node_attr = part[, .(org_name = org_name[1], country = country[1],
                                       type = activity_type[1]), by = org_id],
                  by = "org_id")

## Work on a manageable subgraph. Brokerage needs to look at every 2-path, so
## cost grows with degree^2: on the full 15,000-organisation network the roles
## below take minutes, on the 2,400 most connected ones a few seconds.
g <- induced_subgraph(g_org, V(g_org)[degree(g_org) >= 50])
g <- giant(g)
c(nodes = vcount(g), edges = ecount(g), mean_degree = round(mean(degree(g))))

## --- A1. Burt's ego-network measures --------------------------------------- ##
## igraph gives you constraint(); effective_size() is one of the helpers defined
## in 00_setup.R - open it, it is five lines of sparse-matrix algebra.
V(g)$constraint <- constraint(g)                     # LOW  = many holes
V(g)$eff_size   <- effective_size(g)                 # HIGH = many non-redundant
V(g)$efficiency <- V(g)$eff_size / degree(g)         # per-contact yield
V(g)$betw       <- betweenness(g, weights = NA, normalized = TRUE)

brok <- as.data.table(as_data_frame(g, what = "vertices"))
brok[, degree := degree(g)]
round(cor(brok[, .(degree, betw, constraint, eff_size, efficiency)]), 2)
## note: effective size correlates ~1 with degree, efficiency does not. If you
## want "brokerage net of size", use efficiency or constraint, not effective size.

brok[order(constraint)][1:10, .(org_name, country, type, degree,
                                constraint = round(constraint, 3),
                                efficiency = round(efficiency, 2))]

## --- A2. Gould & Fernandez brokerage roles --------------------------------- ##
## A node v brokers a 2-path i -> v -> j when i and j are NOT directly tied.
## Given a partition into groups, the role depends on where i, v and j sit:
##
##   coordinator     i, v, j all in the same group        (within-group broker)
##   gatekeeper      i outside, v and j inside            (controls entry)
##   representative  i and v inside, j outside            (controls exit)
##   consultant      i and j in the SAME other group      (itinerant broker)
##   liaison         i, v, j all in different groups      (bridges two others)
##
## Not in igraph. sna::brokerage() does it (with intergraph::asNetwork), but
## brokerage_roles() in 00_setup.R is transparent, has no dependency and runs on
## 10^3-10^4 nodes. Read the function before using it: the two corrections it
## applies (dropping i == j, and dropping pairs that are already tied) are the
## whole definition.

## NOTE on undirected networks: with i -> v -> j read in both directions,
## "gatekeeper" and "representative" are the same count by construction (see the
## table below). Only directed data - citations, trade flows, supply chains -
## separates controlling entry from controlling exit.

## groups = country: who mediates between national research systems?
roles <- brokerage_roles(g, V(g)$country)
roles <- merge(roles, brok[, .(name, org_name, country, type, degree)], by = "name")
roles[order(-liaison)][1:10, .(org_name, country, type, degree,
                               coordinator, gatekeeper, representative,
                               consultant, liaison)]

## --- A3. Raw counts are useless without a benchmark ------------------------ ##
## Every count grows with degree. Compare each node with what it would score if
## the group labels were reshuffled at random (keeping the network fixed).
brokerage_z <- function(g, group, reps = 10) {
  obs <- brokerage_roles(g, group)
  cols <- setdiff(names(obs), "name")
  sims <- lapply(seq_len(reps), function(r)
    as.matrix(brokerage_roles(g, sample(group))[, ..cols]))
  mu <- Reduce(`+`, sims) / reps
  sdv <- sqrt(Reduce(`+`, lapply(sims, function(x) (x - mu)^2)) / max(reps - 1, 1))
  z <- (as.matrix(obs[, ..cols]) - mu) / pmax(sdv, 1e-9)
  cbind(obs[, .(name)], as.data.table(round(z, 2)))
}
## (10 permutations to keep the class moving - about 30 seconds; use 500+ in a
##  paper, and run it once overnight rather than in a loop you watch)
zz <- brokerage_z(g, V(g)$country, reps = 10)
zz <- merge(zz, brok[, .(name, org_name, country, degree)], by = "name")
zz[order(-liaison)][1:10, .(org_name, country, degree, gatekeeper,
                            representative, liaison)]
## Now "liaison" means "more than expected by chance", which is what the theory
## is about. The identity of the top brokers usually changes: show both tables.

## --- A4. Choosing ---------------------------------------------------------- ##
## - fragmented network (co-invention)? betweenness is dominated by component
##   structure; prefer ego-level measures (constraint, efficiency).
## - projected affiliation network? every event is a clique, so constraint is
##   mechanically high; compare against events of the same size.
## - theory about categories (countries, sectors, public/private)? that is
##   exactly Gould-Fernandez; report the z-scores, not the raw counts.
## - weighted networks: betweenness treats weights as DISTANCES (a strong tie is
##   a long detour!). Pass weights = 1/w, or weights = NA to ignore them.

## ===========================================================================
## PART B - COMMUNITY DETECTION
## ===========================================================================
## No algorithm "finds the true communities": each optimises a different
## objective. What matters is that your groups are (i) reproducible,
## (ii) interpretable, (iii) not an artefact of one arbitrary choice.

gc <- giant(g)
W  <- E(gc)$weight

## --- B1. The main families ------------------------------------------------- ##
algos <- list(
  louvain      = function(x) cluster_louvain(x, weights = E(x)$weight),
  leiden_cpm   = function(x) cluster_leiden(x, objective_function = "CPM",
                                            resolution = 0.05, weights = E(x)$weight),
  fast_greedy  = function(x) cluster_fast_greedy(x, weights = E(x)$weight),
  walktrap     = function(x) cluster_walktrap(x, weights = E(x)$weight),
  infomap      = function(x) cluster_infomap(x, e.weights = E(x)$weight),
  label_prop   = function(x) cluster_label_prop(x, weights = E(x)$weight)
)
comp <- rbindlist(lapply(names(algos), function(a) {
  t0 <- proc.time()[["elapsed"]]
  cl <- algos[[a]](gc)
  data.table(algorithm = a, communities = length(unique(membership(cl))),
             modularity = round(modularity(gc, membership(cl), weights = W), 3),
             largest_share = round(max(table(membership(cl))) / vcount(gc), 2),
             seconds = round(proc.time()[["elapsed"]] - t0, 1))
}))
comp[order(-modularity)]
## In this dense projected network infomap collapses almost everything into one
## module and label propagation into two: both are designed for sparser graphs.
## That is information, not failure - report it instead of hiding it.
## Read this table before choosing: infomap tends to many small communities,
## modularity-based methods to few large ones, label propagation is fast but
## unstable, edge-betweenness (not run here) is O(n*m) - forget it above ~1,000 nodes.

## --- B2. Weights and direction matter -------------------------------------- ##
c(weighted   = modularity(gc, membership(cluster_louvain(gc, weights = W)), weights = W),
  unweighted = modularity(gc, membership(cluster_louvain(gc, weights = NA))))
## Decide explicitly: is a tie of weight 8 eight times "closer" than a tie of 1?
## For counts of shared documents, usually yes; for correlations, usually no.

## --- B3. The resolution parameter is a research choice --------------------- ##
res_scan <- rbindlist(lapply(c(0.5, 1, 2, 4), function(r) {
  cl <- cluster_louvain(gc, weights = W, resolution = r)
  data.table(resolution = r, communities = length(unique(membership(cl))),
             modularity = round(modularity(gc, membership(cl), weights = W), 3),
             largest_share = round(max(table(membership(cl))) / vcount(gc), 2))
}))
res_scan
## Modularity has a RESOLUTION LIMIT: it cannot see communities smaller than
## ~sqrt(2m). Scanning the resolution and reporting the range is the honest move.

## --- B4. Stability: run it again ------------------------------------------- ##
## Louvain and Leiden are stochastic. Are your communities a property of the
## network or of the seed?
parts <- lapply(1:20, function(s) { set.seed(s); membership(cluster_louvain(gc, weights = W)) })
mods  <- sapply(parts, function(m) modularity(gc, m, weights = W))
c(min = round(min(mods), 4), max = round(max(mods), 4),
  n_comm_min = min(sapply(parts, function(m) length(unique(m)))),
  n_comm_max = max(sapply(parts, function(m) length(unique(m)))))

## agreement between runs: NMI = 1 means identical partitions
nmi <- outer(seq_along(parts), seq_along(parts), Vectorize(function(i, j)
  compare(parts[[i]], parts[[j]], method = "nmi")))
round(c(mean_nmi = mean(nmi[upper.tri(nmi)]), min_nmi = min(nmi[upper.tri(nmi)])), 3)
## compare() also gives "adjusted.rand" and "vi" (variation of information).

## --- B5. Consensus: keep what survives ------------------------------------- ##
## co-membership frequency across runs, then cluster the consensus matrix
co <- Reduce(`+`, lapply(parts, function(m) outer(m, m, "==") * 1)) / length(parts)
mean(co[upper.tri(co)] > 0 & co[upper.tri(co)] < 1)   # share of unstable pairs
g_cons <- graph_from_adjacency_matrix(co * (co >= 0.9), mode = "undirected",
                                      weighted = TRUE, diag = FALSE)
cons <- components(g_cons)
c(consensus_groups = cons$no, largest = max(cons$csize))
## Nodes that never travel together are the ones you can safely talk about.

## --- B6. Interpretation is the actual result ------------------------------- ##
## A community is only useful if you can say what it IS. Cross it with attributes.
cl <- cluster_louvain(gc, weights = W)
memb <- data.table(name = V(gc)$name, country = V(gc)$country,
                   type = V(gc)$type, comm = as.integer(membership(cl)))
top <- memb[, .N, by = comm][order(-N)][1:5]$comm
memb[comm %in% top, .(orgs = .N,
                      top_country = names(which.max(table(country))),
                      country_hhi = round(sum(prop.table(table(country))^2), 3),
                      pct_private = round(100 * mean(type == "PRC")),
                      pct_academic = round(100 * mean(type %in% c("HES", "REC")))),
     by = comm][order(-orgs)]

## --- B6b. Look at them: the same partition, drawn --------------------------- ##
## A table of community sizes hides what a picture shows immediately: whether the
## communities are compact blocks or a hairball cut arbitrarily in five.
## Detect on the FULL network, then draw a readable subgraph of it - never detect
## on the subgraph you happen to be able to plot.
cl_wt <- cluster_walktrap(gc, weights = W)          # a second opinion
c(louvain = length(unique(membership(cl))),
  walktrap = length(unique(membership(cl_wt))),
  nmi = round(compare(membership(cl), membership(cl_wt), method = "nmi"), 2))

## a drawable backbone: strong ties only, and no leftover isolates
g_plot <- delete_edges(gc, E(gc)[weight < 4])
g_plot <- induced_subgraph(g_plot, V(g_plot)[degree(g_plot) > 0])
g_plot <- giant(g_plot)
c(nodes = vcount(g_plot), edges = ecount(g_plot))

## attach everything we want to draw BEFORE building the layout, and keep only
## the largest communities as named colours (a legend with 31 entries is noise)
top6 <- function(x) {
  keep <- names(sort(table(x), decreasing = TRUE))[1:6]
  factor(ifelse(x %in% keep, paste0("c", x), "other"),
         levels = c(paste0("c", keep), "other"))
}
V(g_plot)$deg      <- degree(g_plot)
V(g_plot)$louvain  <- top6(as.integer(membership(cl))[match(V(g_plot)$name, V(gc)$name)])
V(g_plot)$walktrap <- top6(as.integer(membership(cl_wt))[match(V(g_plot)$name, V(gc)$name)])

## ONE layout object, reused by both plots: the nodes then sit in exactly the
## same place and the only thing that changes between the two figures is colour
set.seed(42)
lay_comm <- create_layout(g_plot, layout = "stress")

comm_plot <- function(lay, colour_var, title, subtitle) {
  ggraph(lay) +
    geom_edge_link0(aes(edge_width = weight), edge_colour = "grey85", edge_alpha = .6) +
    scale_edge_width(range = c(0.1, 1.2), guide = "none") +
    geom_node_point(aes(size = deg, fill = .data[[colour_var]]), shape = 21,
                    colour = "white", stroke = 0.3) +
    scale_size(range = c(1.5, 8), guide = "none") +
    theme_graph(base_family = "sans") +
    labs(title = title, subtitle = subtitle, fill = "community")
}

comm_plot(lay_comm, "louvain",
          "Louvain communities, CORDIS organisations (ties with 4+ shared projects)",
          "communities detected on the full network; the six largest are coloured")

## The same nodes, the same positions, a different algorithm. Flip between the
## two: where the colours disagree, your "communities" are a property of the
## algorithm rather than of the network.
comm_plot(lay_comm, "walktrap",
          "The same network and layout, walktrap instead of Louvain",
          paste0("agreement with Louvain: NMI = ",
                 round(compare(membership(cl), membership(cl_wt), method = "nmi"), 2)))

## And the community-level view: communities as nodes, ties between them
## aggregated. Useful when the node-level picture is unreadable - which, above a
## couple of thousand nodes, it usually is.
g_meta <- contract(gc, membership(cl), vertex.attr.comb = "ignore")
g_meta <- simplify(g_meta, edge.attr.comb = list(weight = "sum"))
V(g_meta)$members <- as.integer(table(membership(cl)))
V(g_meta)$name <- as.character(seq_len(vcount(g_meta)))
g_meta <- induced_subgraph(g_meta, V(g_meta)[members >= 20])

ggraph(g_meta, layout = "stress") +
  geom_edge_link0(aes(edge_width = weight), edge_colour = "grey80") +
  scale_edge_width(range = c(0.2, 3.5), guide = "none") +
  geom_node_point(aes(size = members), fill = "#2c7fb8", shape = 21, colour = "white") +
  geom_node_text(aes(label = paste0("c", name, "\n", members)), size = 2.6) +
  scale_size(range = c(4, 18), guide = "none") +
  theme_graph(base_family = "sans") +
  labs(title = "The same partition seen from above: communities as nodes",
       subtitle = "node size = organisations in the community, edge width = ties between them")

## --- B7. What to report ---------------------------------------------------- ##
## In the paper, one sentence must contain: algorithm + implementation and
## version + weights used + resolution + seed/number of runs + modularity +
## number and size distribution of communities + one robustness check
## (another algorithm, or the consensus above). Anything less is not replicable.
##
## And remember what modularity cannot do: it always finds a partition, even in
## a random graph. Benchmark against a degree-preserving rewiring:
rnd <- rewire(gc, with = keeping_degseq(niter = 10 * ecount(gc)))
c(observed = round(modularity(gc, membership(cluster_louvain(gc, weights = NA))), 3),
  rewired  = round(modularity(rnd, membership(cluster_louvain(rnd))), 3))
