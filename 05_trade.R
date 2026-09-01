#------------------------------------------------------------------------------#
#  DAISY 2026 - Network analysis: indicators, mapping and empirics
#  with innovation data
#  BLOCK 4 of the session - TRADE AND GLOBAL VALUE CHAINS:
#  observed, directed, valued networks
#
#  Everything so far was an *inferred* tie: two actors shared a document. Trade
#  data are the opposite case, and they break most of our habits:
#
#    - the tie is OBSERVED and has a VALUE (millions of USD), not a count;
#    - it is DIRECTED (i sells to j is not j sells to i);
#    - the network is essentially COMPLETE: almost every pair trades something.
#      Density carries no information; the whole signal is in the weights.
#      Filtering stops being cosmetic and becomes part of the method.
#
#  Data: OECD TiVA (FDVA table) - value added of country i embodied in the final
#  demand of country j, 80 economies, 1995/2005/2015/2022, total and by industry.
#  Caveat to state out loud: TiVA is not a measurement, it is the output of an
#  inter-country input-output model. The network inherits its assumptions.
#------------------------------------------------------------------------------#
source("00_setup.R")

## ===========================================================================
## 1. Value added flows between countries
## ===========================================================================
va  <- fread(daisy_data("tiva_va_bilateral.csv.gz"))
reg <- fread(daisy_data("country_regions.csv"))
va

## the diagonal is domestic value added in domestic final demand: not a tie
va <- va[source != destination]
va[year == 2022, .(flows = .N, total_bn = round(sum(va_musd) / 1000))]
va[year == 2022][order(-va_musd)][1:8]

## how much of world trade in value added is a handful of pairs?
va[year == 2022, .(top10_share = round(100 * sum(sort(va_musd, decreasing = TRUE)[1:10]) /
                                        sum(va_musd)))]

## ===========================================================================
## 2. A directed, weighted graph
## ===========================================================================
g22 <- graph_from_data_frame(va[year == 2022, .(source, destination, weight = va_musd)],
                             directed = TRUE,
                             vertices = reg[, .(name = iso3, region)])
g22
c(nodes = vcount(g22), edges = ecount(g22), density = round(edge_density(g22), 3))
## density ~ 1: the topology is a complete digraph, so degree is useless here...

## ... and strength is everything
V(g22)$out_str <- strength(g22, mode = "out")      # value added SOLD abroad
V(g22)$in_str  <- strength(g22, mode = "in")       # foreign VA ABSORBED
nodes <- as.data.table(as_data_frame(g22, what = "vertices"))
nodes[, `:=`(net = in_str - out_str,
             share_out = round(100 * out_str / sum(out_str), 1))]
nodes[order(-out_str)][1:12, .(name, region, out_bn = round(out_str/1000),
                               in_bn = round(in_str/1000), share_out)]

## Concentration: is a country's foreign demand diversified or hostage to one market?
el22 <- as.data.table(as_data_frame(g22, what = "edges"))
hhi <- el22[, .(hhi = sum((weight / sum(weight))^2),
                top_market = to[which.max(weight)],
                top_share = round(100 * max(weight) / sum(weight))), by = from]
hhi[order(-hhi)][1:10]
hhi[from %in% c("DEU","ITA","CHN","USA","MEX","IRL")]
## Mexico and Canada are structurally exposed to one market; Germany and Italy
## are not. Same network, a country-level risk indicator.

## Dyadic asymmetry: who is upstream of whom
dy <- merge(el22, el22, by.x = c("from","to"), by.y = c("to","from"))
dy <- dy[from < to, .(a = from, b = to, a_to_b = weight.x, b_to_a = weight.y)]
dy[, imbalance := round((a_to_b - b_to_a) / (a_to_b + b_to_a), 2)]
dy[a_to_b + b_to_a > 50000][order(-imbalance)][1:8]

## ===========================================================================
## 3. Filtering a valued network: four options, four different networks
## ===========================================================================
## (a) absolute threshold - simple, arbitrary, biased against small countries
f_abs <- el22[weight >= 5000]

## (b) top-k destinations of each country - guarantees every node survives
f_topk <- el22[order(from, -weight)][, head(.SD, 5), by = from]

## (c) share of the source's total exports of value added
f_share <- el22[, .SD[weight / sum(weight) >= 0.05], by = from]

## (d) DISPARITY FILTER (Serrano, Boguna & Vespignani 2009, PNAS) - keeps the
##     links that are significantly stronger than a random allocation of a
##     node's strength across its ties.
disparity_filter <- function(el, alpha = 0.05) {
  d <- copy(as.data.table(el))
  d[, `:=`(k = .N, s = sum(weight)), by = from]
  d[, p_ij := weight / s]
  d[, alpha_ij := (1 - p_ij)^(k - 1)]              # p-value under the null
  d[k > 1 & alpha_ij < alpha]
}
f_disp <- disparity_filter(el22, alpha = 0.01)

rbindlist(list(
  data.table(filter = "none",            edges = nrow(el22)),
  data.table(filter = "weight >= 5bn",   edges = nrow(f_abs)),
  data.table(filter = "top 5 per country", edges = nrow(f_topk)),
  data.table(filter = ">= 5% of exports", edges = nrow(f_share)),
  data.table(filter = "disparity a=0.01", edges = nrow(f_disp))))

## Does the choice change the answer? Compare betweenness rankings.
bt <- function(el) {
  gg <- graph_from_data_frame(el[, .(from, to, weight)], directed = TRUE,
                              vertices = reg[, .(name = iso3)])
  betweenness(gg, weights = NA)
}
btw <- data.table(country = reg$iso3, abs5 = bt(f_abs), topk = bt(f_topk),
                  share5 = bt(f_share), disparity = bt(f_disp))
round(cor(btw[, -1], method = "spearman"), 2)
btw[order(-disparity)][1:10]

## The absolute threshold makes small open economies disappear; the disparity
## filter keeps them. State the filter in the paper - it IS a modelling choice.

## ===========================================================================
## 4. Map it (you can only draw a valued network after filtering it)
## ===========================================================================
gp <- graph_from_data_frame(f_disp[, .(from, to, weight)], directed = TRUE,
                            vertices = reg[, .(name = iso3, region)])
gp <- induced_subgraph(gp, V(gp)[degree(gp) > 0])
V(gp)$out_str <- strength(gp, mode = "out")

ggraph(gp, layout = "stress") +
  geom_edge_fan(aes(edge_width = weight, edge_alpha = weight),
                edge_colour = "grey55",
                arrow = arrow(angle = 15, length = unit(0.09, "inches"),
                              type = "closed"),
                start_cap = circle(2, "mm"), end_cap = circle(3, "mm")) +
  scale_edge_width(range = c(0.1, 1.6)) + scale_edge_alpha(range = c(0.15, 0.7)) +
  geom_node_point(aes(size = out_str, fill = region), shape = 21, colour = "white") +
  geom_node_text(aes(label = name), size = 2.6, repel = TRUE) +
  scale_size(range = c(2, 12)) +
  theme_graph(base_family = "sans") +
  labs(title = "Value added embodied in foreign final demand, 2022",
       subtitle = "OECD TiVA, disparity-filter backbone (alpha = 0.01)",
       fill = "region") + guides(size = "none", edge_width = "none",
                                 edge_alpha = "none")

## ===========================================================================
## 5. Communities = trade blocs, and how regional they are
## ===========================================================================
## Louvain needs an undirected graph: symmetrise the flows (i<->j = i->j + j->i)
symmetrise <- function(el) {
  d <- copy(as.data.table(el))
  d[, `:=`(a = pmin(from, to), b = pmax(from, to))]
  d[, .(weight = sum(weight)), by = .(from = a, to = b)]
}

sym22 <- symmetrise(el22)
gu <- graph_from_data_frame(sym22, directed = FALSE,
                           vertices = reg[, .(name = iso3, region)])
cl_raw <- cluster_louvain(gu, weights = E(gu)$weight)
c(blocs = length(unique(membership(cl_raw))),
  modularity = round(modularity(cl_raw), 3),
  largest_share = round(max(table(membership(cl_raw))) / vcount(gu), 2))

## Two blocks, and they are essentially "around the USA" and "around Germany".
## On a COMPLETE VALUED network, modularity is driven by the size of the nodes:
## big economies trade a lot with everybody, so the partition mostly recovers
## who is big. Before looking for structure, take size out.

## Revealed trade intensity: observed flow / flow expected from the two
## countries' sizes (the same normalisation logic as the CORDIS index in block 2)
normalise <- function(sym) {
  d <- copy(sym)
  str <- rbind(d[, .(c = from, w = weight)], d[, .(c = to, w = weight)])[
    , .(s = sum(w)), by = c]
  W <- sum(d$weight)
  d <- merge(merge(d, str, by.x = "from", by.y = "c"),
             str, by.x = "to", by.y = "c", suffixes = c("_f", "_t"))
  d[, weight := weight / (s_f * s_t / (2 * W))]
  d[, .(from, to, weight)]
}

blocs <- function(yr, normalised = TRUE) {
  e <- va[year == yr, .(from = source, to = destination, weight = va_musd)]
  sym <- symmetrise(e)
  if (normalised) sym <- normalise(sym)
  gg <- graph_from_data_frame(sym, directed = FALSE,
                              vertices = reg[, .(name = iso3, region)])
  cl <- cluster_louvain(gg, weights = E(gg)$weight)
  list(g = gg, cl = cl,
       stats = data.table(year = yr, normalised = normalised,
                          blocs = length(unique(membership(cl))),
                          modularity = round(modularity(cl), 3),
                          nmi_with_geography = round(
                            compare(membership(cl), as.integer(factor(V(gg)$region)),
                                    method = "nmi"), 3)))
}

rbind(blocs(1995, normalised = FALSE)$stats, blocs(2022, normalised = FALSE)$stats,
      blocs(1995)$stats, blocs(2022)$stats)

## With raw weights: two blocks, no change in 27 years - the size effect swamps
## everything. With normalised weights: more, smaller blocks that align much more
## closely with geography, and now the 1995 -> 2022 comparison is informative.
## This is exactly the question in Fusillo, Montresor & Vittucci Marzetti (2024):
## have national and regional boundaries really faded away?

b22 <- blocs(2022)
memb <- data.table(iso3 = V(b22$g)$name, region = V(b22$g)$region,
                   bloc = as.integer(membership(b22$cl)))
va_out <- va[year == 2022, .(va = sum(va_musd)), by = .(iso3 = source)]
memb <- merge(memb, va_out, by = "iso3")
memb[, .(members = .N, va_bn = round(sum(va) / 1000),
         regions = uniqueN(region),
         core = paste(head(iso3[order(-va)], 4), collapse = " ")),
     by = bloc][order(-va_bn)]

## Communities on the DIRECTED graph instead (infomap follows the flow):
im <- cluster_infomap(g22, e.weights = E(g22)$weight)
c(louvain_blocs = length(unique(membership(b22$cl))),
  infomap_blocs = length(unique(membership(im))),
  agreement_nmi = round(compare(membership(b22$cl), membership(im), method = "nmi"), 3))

## Infomap puts everything in one module: a random walk on a complete weighted
## digraph never gets trapped anywhere. Not a bug - a property of the data.
## See 06_brokerage_communities.R for how to choose an algorithm and report it.

## ===========================================================================
## 6. Brokerage in value chains: who intermediates between regions?
## ===========================================================================
## Here the data are DIRECTED, so "gatekeeper" (controls what enters my region)
## and "representative" (controls what leaves it) are finally different things.
## Run it on the disparity backbone: raw counts on a complete graph are meaningless.
## a slightly looser backbone (alpha = 0.05) leaves enough 2-paths to classify
gb <- graph_from_data_frame(disparity_filter(el22, alpha = 0.05)[, .(from, to, weight)],
                            directed = TRUE, vertices = reg[, .(name = iso3, region)])
gb <- induced_subgraph(gb, V(gb)[degree(gb) > 0])
c(nodes = vcount(gb), edges = ecount(gb))

roles <- brokerage_roles(gb, V(gb)$region)
roles <- merge(roles, data.table(name = V(gb)$name, region = V(gb)$region,
                                 out_str = strength(gb, mode = "out")), by = "name")
head(roles[order(-liaison)], 10)[, .(name, region, coordinator, gatekeeper,
                                     representative, consultant, liaison)]

## normalise by degree: who brokers MORE than their size implies?
roles[, total := coordinator + gatekeeper + representative + consultant + liaison]
head(roles[total >= 10][order(-liaison / total)], 10)[,
     .(name, region, share_liaison = round(liaison / total, 2),
       share_gatekeeper = round(gatekeeper / total, 2),
       share_coordinator = round(coordinator / total, 2), total)]
## Small open economies and re-export hubs (NLD, BEL, SGP, HKG, MEX) live off
## intermediation; the large ones broker within their own bloc.

## ===========================================================================
## 7. IF WE HAVE TIME - one method, eight industries
## ===========================================================================
## The same pipeline, run per industry: the geography of value chains is not the
## same for food, cars, electronics and business services.
ind <- fread(daisy_data("tiva_va_by_industry.csv.gz"))
lab <- fread(daisy_data("tiva_industry_labels.csv"))
ind <- ind[source != destination]

by_ind <- rbindlist(lapply(unique(ind$industry), function(i) {
  e  <- ind[industry == i, .(from = source, to = destination, weight = va_musd)]
  gu <- graph_from_data_frame(normalise(symmetrise(e)), directed = FALSE,
                             vertices = reg[, .(name = iso3, region)])
  cl <- cluster_louvain(gu, weights = E(gu)$weight)
  s  <- e[, .(w = sum(weight)), by = from][match(V(gu)$name, from), w]
  data.table(industry = i,
             va_bn = round(sum(e$weight) / 1000),
             top3 = paste(V(gu)$name[order(-s)][1:3], collapse = " "),
             hhi = round(sum((s / sum(s))^2), 3),
             blocs = length(unique(membership(cl))),
             modularity = round(modularity(cl), 3),
             nmi_geography = round(compare(membership(cl),
                                   as.integer(factor(V(gu)$region)), "nmi"), 3))
}))
merge(by_ind, lab, by = "industry")[order(-va_bn)]
## Compare "modularity" (how bloc-structured the industry is) with
## "nmi_geography" (whether those blocs are geographic). Electronics and
## transport equipment are the regionalised ones; services are not.

## ===========================================================================
## 8. IF WE HAVE TIME - gross trade or value added? (BACI vs TiVA)
## ===========================================================================
## BACI (CEPII) records gross bilateral flows of goods by HS6 product: what
## crosses the border. TiVA records where the value was actually created. For
## some countries the two tell very different stories - and the difference IS
## the network position.
bac <- fread(daisy_data("baci_bilateral_2023.csv.gz"))
bac[order(-exports_musd)][1:6]

gross <- bac[exporter != importer, .(gross_bn = sum(exports_musd) / 1000), by = .(iso3 = exporter)]
vadd  <- va[year == 2022, .(va_bn = sum(va_musd) / 1000), by = .(iso3 = source)]
cmp   <- merge(gross, vadd, by = "iso3")
cmp[, ratio := round(va_bn / gross_bn, 2)]          # VA generated per $ exported
cor(cmp$gross_bn, cmp$va_bn, method = "spearman")

cmp[gross_bn > 100][order(ratio)][1:10]             # pure transit / assembly
cmp[gross_bn > 100][order(-ratio)][1:10]            # value created at home
## Low ratio = you ship a lot but much of the value is foreign (Vietnam, Mexico,
## Belgium, Netherlands: assembly platforms and re-export hubs). Interpreting a
## gross-trade network as a network of "who produces what" is a measurement error
## with a name: double counting.

## The same contrast at the level of a single tie
gr_pairs <- bac[, .(a = pmin(exporter, importer), b = pmax(exporter, importer),
                    gross = exports_musd)][, .(gross = sum(gross)), by = .(a, b)]
va_pairs <- dy[, .(a, b, va = a_to_b + b_to_a)]
pairs <- merge(gr_pairs, va_pairs, by = c("a", "b"))
pairs[va > 20000][, ratio := round(va / gross, 2)][order(ratio)][1:8]

## ===========================================================================
## 9. And the product space
## ===========================================================================
## The exporter x product matrix of BACI feeds exactly the machinery of
## 04_indicators.R: revealed comparative advantage, proximity between products,
## the product space, and the complexity indices. We build it there, so that the
## three category systems of this session - technologies (patents), topics
## (CORDIS) and products (trade) - sit side by side in one script.
