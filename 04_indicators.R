#------------------------------------------------------------------------------#
#  DAISY 2026 - Hands-on: Network analysis with innovation data
#  BLOCK 5 of the session (~20 min) - NETWORKS AS MEASUREMENT DEVICES
#  ("indirect" uses: when you do not study the network, you use it to build a
#   variable)
#
#  Idea: a knowledge base has a *co-relational* structure. Represent technologies
#  as nodes and their joint use as links, and the network becomes a measurement
#  instrument for concepts that have no direct observable counterpart:
#  relatedness, variety, coherence, complexity, diversification potential.
#
#  Data: OECD REGPAT green patents (Y02/Y04S), EU NUTS regions, CPC 4-digit
#  classes, two periods (2010-2014, 2015-2019).
#------------------------------------------------------------------------------#
source("00_setup.R")

## ===========================================================================
## 1. The regional technology portfolios
## ===========================================================================
rt  <- fread(daisy_data("green_region_tech_EU.csv.gz"))
def <- fread(daisy_data("cpc4_def.csv"))              # CPC4 labels
rt

## fractional counts: a patent is split across its regions and its CPC classes,
## so that every patent contributes exactly 1 to the world total
rt[, sum(n_pat), by = period]

## NUTS-3 -> NUTS-2 (the level at which regional innovation is usually studied)
rt[, nuts2 := substr(reg_code, 1, 4)]
reg_tech <- rt[, .(n_pat = sum(n_pat)), by = .(nuts2, ctry_code, cpc4, period)]

## keep the recent period and regions/technologies with enough mass
d <- reg_tech[period == "2015-2019"]
big_reg  <- d[, .(tot = sum(n_pat)), by = nuts2][tot >= 20, nuts2]
big_tech <- d[, .(tot = sum(n_pat)), by = cpc4][tot >= 20, cpc4]
d <- d[nuts2 %in% big_reg & cpc4 %in% big_tech]
d[, .(regions = uniqueN(nuts2), technologies = uniqueN(cpc4))]

## region x technology matrix
X <- as.matrix(dcast(d, nuts2 ~ cpc4, value.var = "n_pat", fill = 0),
               rownames = "nuts2")
dim(X)

## ===========================================================================
## 2. Revealed Technological Advantage: from counts to specialisation
## ===========================================================================
## RTA_rt = (X_rt / X_r.) / (X_.t / X_..)   ("Balassa index")
RTA <- (X / rowSums(X)) / rep(colSums(X) / sum(X), each = nrow(X))
M   <- (RTA >= 1) * 1                        # binary specialisation matrix
mean(M)                                      # density of the two-mode network

diversity <- rowSums(M)                      # n. of technologies of a region
ubiquity  <- colSums(M)                      # n. of regions with that technology
sort(diversity, decreasing = TRUE)[1:10]
merge(data.table(cpc4 = names(ubiquity), ubiquity), def,
      by = "cpc4")[order(ubiquity)][1:8]                # most exclusive classes

## ===========================================================================
## 3. Two ways of measuring RELATEDNESS between technologies
## ===========================================================================
## (A) co-classification inside patents (the "knowledge space" proper):
##     two classes are related if inventors combine them in the same document
cooc <- fread(daisy_data("green_tech_cooc_EU.csv.gz"))
npat <- fread(daisy_data("green_tech_npat_EU.csv"))
techs <- colnames(X)
cooc  <- cooc[cpc4_i %in% techs & cpc4_j %in% techs]
C <- sparseMatrix(i = match(cooc$cpc4_i, techs), j = match(cooc$cpc4_j, techs),
                  x = cooc$n_cooc, dims = c(length(techs), length(techs)),
                  dimnames = list(techs, techs), symmetric = FALSE)
C <- as.matrix(C + t(C))
n_t <- setNames(npat$n_pat, npat$cpc4)[techs]
## association strength (Van Eck & Waltman): observed / expected co-occurrence
Phi_pat <- C / outer(n_t, n_t) * sum(n_t) / 2
Phi_pat[!is.finite(Phi_pat)] <- 0

## (B) co-specialisation across regions (Hidalgo et al. "proximity"):
##     two technologies are related if the same regions are good at both
Co <- t(M) %*% M
Phi_reg <- Co / outer(ubiquity, ubiquity, pmax)      # min conditional probability
Phi_reg[!is.finite(Phi_reg)] <- 0
diag(Phi_reg) <- 0

## Do the two measures agree? (they answer different questions!)
iu <- upper.tri(Phi_reg)
cor(Phi_pat[iu], Phi_reg[iu], method = "spearman")

## ===========================================================================
## 4. The green knowledge space, drawn
## ===========================================================================
## keep the strongest links only, otherwise the map is a hairball
thr <- quantile(Phi_pat[iu], 0.98)
A <- Phi_pat * (Phi_pat >= thr)
g_ks <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
g_ks <- induced_subgraph(g_ks, V(g_ks)[degree(g_ks) > 0])
V(g_ks)$patents <- n_t[V(g_ks)$name]
V(g_ks)$label   <- def$label[match(V(g_ks)$name, def$cpc4)]
V(g_ks)$comm    <- membership(cluster_louvain(g_ks, weights = E(g_ks)$weight))
g_ks

## which technologies bridge the green knowledge space? (candidate GPTs)
sort(betweenness(g_ks, weights = NA), decreasing = TRUE)[1:10]

ggraph(g_ks, layout = "stress") +
  geom_edge_link0(aes(edge_width = weight), edge_colour = "grey85") +
  scale_edge_width(range = c(0.1, 1.5)) +
  geom_node_point(aes(size = patents, fill = factor(comm)), shape = 21, colour = "white") +
  geom_node_text(aes(label = name), size = 2.4, repel = TRUE, max.overlaps = 30) +
  scale_size(range = c(1, 10)) +
  theme_graph(base_family = "sans") + theme(legend.position = "none") +
  labs(title = "Green knowledge space (CPC4 co-classification, EU 2015-2019)")

## ===========================================================================
## 5. Region-level indicators built ON the network
## ===========================================================================
sh <- X / rowSums(X)                  # patent shares of each region

## (a) VARIETY = entropy of the portfolio, decomposed into related (within
##     3-digit CPC groups) and unrelated (between groups) variety
grp <- substr(colnames(X), 1, 3)
H <- function(p) { p <- p[p > 0]; -sum(p * log2(p)) }
variety <- apply(sh, 1, H)
grp_sh  <- t(rowsum(t(sh), grp))                      # shares by 3-digit group
unrelated <- apply(grp_sh, 1, H)                      # between-group entropy
related   <- variety - unrelated                      # within-group entropy

## (b) COHERENCE = average relatedness of the technologies a region holds,
##     weighted by their importance in the portfolio (Nesta & Saviotti)
coherence <- as.numeric(rowSums((sh %*% Phi_pat) * sh))

## (c) RELATEDNESS DENSITY of the technologies the region is NOT (yet) in:
##     how "close" is a new technology to what the region already does?
dens <- (M %*% Phi_reg) / rep(colSums(Phi_reg), each = nrow(M)) * 100
avg_density_out <- rowSums(dens * (1 - M)) / rowSums(1 - M)
## CAVEAT: averaged over a region, relatedness density is almost collinear with
## diversity (see the correlation matrix below). The informative variation is at
## the region x technology level - which is exactly how it is used in section 6.

## (d) COMPLEXITY of the regional portfolio
##     Hidalgo & Hausmann (2009) - applied to technologies by Balland & Rigby
##     (2017). Two equivalent readings of the same bipartite network:
##
##  d.1 METHOD OF REFLECTIONS - iterate "average of my neighbours' average"
reflections <- function(M, iter = 2) {
  kr <- rowSums(M); kt <- colSums(M)
  for (i in seq_len(iter)) {
    kr_new <- as.numeric((M %*% kt) / rowSums(M))
    kt_new <- as.numeric((t(M) %*% kr) / colSums(M))
    kr <- kr_new; kt <- kt_new
  }
  list(kr = setNames(kr, rownames(M)), kt = setNames(kt, colnames(M)))
}
## careful with the parity of the iteration: odd orders measure average
## UBIQUITY (high = simple), even orders average DIVERSITY (high = complex)
r1 <- reflections(M, 1); r2 <- reflections(M, 2)
c(order1_vs_diversity = cor(r1$kr, diversity),
  order2_vs_diversity = cor(r2$kr, diversity))

##  d.2 EIGENVECTOR FORM (what the Atlas of Economic Complexity computes)
##      Mtilde = D^-1 M U^-1 M' ; complexity = 2nd eigenvector, standardised.
##      complexity() in 00_setup.R does it for both sides of the matrix at once -
##      open it and read the SIGN conventions, they are where mistakes happen.
cx  <- complexity(M)
kci <- cx$actor        # regional knowledge complexity
tci <- cx$category     # technological complexity
c(kci_vs_reflections = round(cor(kci[names(r2$kr)], r2$kr), 2),
  tci_vs_ubiquity    = round(cor(tci, ubiquity[names(tci)]), 2))
## the second correlation must be NEGATIVE: complex technologies are held by few
## regions. If it is positive, your sign convention is upside down.

## most and least complex green technologies - and their ubiquity, to show that
## complexity is NOT just the inverse of ubiquity: it is second-order. A class
## held by few regions that are themselves poorly diversified is not complex.
tech_cx <- merge(data.table(cpc4 = names(tci), tci,
                            ubiquity = ubiquity[names(tci)]), def, by = "cpc4")
tech_cx[order(-tci)][1:6, .(cpc4, tci = round(tci, 2), ubiquity, label)]
tech_cx[order(tci)][1:6,  .(cpc4, tci = round(tci, 2), ubiquity, label)]

## put everything together: one row per region, ready for a regression
ind <- data.table(nuts2 = rownames(X),
                  patents = rowSums(X),
                  diversity, variety, related_variety = related,
                  unrelated_variety = unrelated,
                  coherence, relatedness_density = avg_density_out,
                  complexity = kci[rownames(X)])
ind[, country := substr(nuts2, 1, 2)]
ind[order(-complexity)][1:12]
round(cor(ind[, .(patents, diversity, variety, related_variety,
                  unrelated_variety, coherence, relatedness_density, complexity)]), 2)

fwrite(ind, "output_region_knowledge_indicators.csv")

ggplot(ind, aes(log(patents), complexity, label = nuts2)) +
  geom_point(aes(size = variety), alpha = .6, colour = "#2c7fb8") +
  geom_text(size = 2.4, vjust = -1, check_overlap = TRUE) +
  labs(x = "log green patents", y = "complexity of the portfolio (KCI)",
       size = "variety") + theme_minimal()

## ===========================================================================
## 6. IF WE HAVE TIME - does relatedness predict diversification?
## ===========================================================================
## The canonical evolutionary-economic-geography test: regions enter new
## technologies that are related to what they already do. We have two periods.
d0 <- reg_tech[period == "2010-2014" & nuts2 %in% rownames(X) & cpc4 %in% colnames(X)]
X0 <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
X0[cbind(d0$nuts2, d0$cpc4)] <- d0$n_pat
RTA0 <- (X0 / pmax(rowSums(X0), 1)) / rep(colSums(X0) / sum(X0), each = nrow(X0))
M0 <- (RTA0 >= 1) * 1
M0[!is.finite(M0)] <- 0

Co0 <- t(M0) %*% M0
Phi0 <- Co0 / outer(pmax(colSums(M0), 1), pmax(colSums(M0), 1), pmax)
Phi0[!is.finite(Phi0)] <- 0; diag(Phi0) <- 0
dens0 <- (M0 %*% Phi0) / rep(pmax(colSums(Phi0), 1e-9), each = nrow(M0)) * 100

entry <- data.table(
  nuts2 = rep(rownames(M), times = ncol(M)),
  cpc4  = rep(colnames(M), each = nrow(M)),
  had   = as.vector(M0), has = as.vector(M),
  density0 = as.vector(dens0))
entry <- entry[had == 0]                       # only technologies NOT held before
entry[, entered := as.integer(has == 1)]
entry[, mean(entered), by = .(density_bin = cut(density0, breaks = c(-1, 5, 10, 20, 100)))][order(density_bin)]

summary(glm(entered ~ density0, data = entry, family = binomial))$coefficients
## => the probability of entering a new green technology increases with the
##    relatedness density of that technology to the regional portfolio: the
##    network is the measurement device behind the "principle of relatedness".
