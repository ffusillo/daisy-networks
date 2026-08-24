#------------------------------------------------------------------------------#
#  DAISY Summer School 2026 - Hands-on session
#  NETWORK ANALYSIS WITH INNOVATION AND SUSTAINABILITY DATA
#  Fabrizio Fusillo (University of Turin)
#
#  00_setup.R : packages, data access, and the helper functions we reuse in every
#               block. Run this first (in VS Code, RStudio or Colab).
#
#  The scripts are numbered by DATA SOURCE; the session runs them in this order:
#     01_patents -> 02_cordis -> 03_publications -> 05_trade -> 04_indicators
#  (trade data come before the indicator block because they introduce valued
#   networks; 06_brokerage_communities is a reference toolbox for self-study)
#------------------------------------------------------------------------------#

## ---------------------------------------------------------------------------
## 1. Packages
## ---------------------------------------------------------------------------
## On Linux (Google Colab) we point R to the Posit binary repository, otherwise
## every install compiles from source and takes ages.
if (Sys.info()[["sysname"]] == "Linux") {
  codename <- tryCatch({
    os <- readLines("/etc/os-release", warn = FALSE)
    sub('.*=', '', grep("^VERSION_CODENAME=", os, value = TRUE))
  }, error = function(e) "jammy")
  if (length(codename) == 0 || codename == "") codename <- "jammy"
  options(repos = c(CRAN = sprintf(
    "https://packagemanager.posit.co/cran/__linux__/%s/latest", codename)))
} else {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
options(timeout = 1800)   # the 60s default is not enough to download bulk data

pkgs <- c("data.table",   # fast data handling (the workhorse for raw big files)
          "igraph",       # network analysis
          "Matrix",       # sparse matrices: two-mode -> one-mode projections
          "ggplot2",      # plots
          "ggraph",       # network visualisation, ggplot2 grammar
          "jsonlite",     # REST APIs (OpenAlex)
          "R.utils")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

## optional packages, only used in clearly marked "if you have time" chunks
## install.packages(c("sna", "intergraph", "graphlayouts"))

setDTthreads(0)           # use all available cores
set.seed(20260907)        # layouts and community detection are stochastic

## ---------------------------------------------------------------------------
## 2. Where is the data?
## ---------------------------------------------------------------------------
## The session works with small pre-processed extracts (~5 MB in total) of
## four sources. Locally they sit in ./data ; in Colab they are downloaded once
## from the course repository. Everything is read through daisy_data().
DAISY_URL <- Sys.getenv("DAISY_DATA_URL",
  "https://raw.githubusercontent.com/ffusillo/daisy-networks/main/data/")

daisy_data <- function(file) {
  local <- file.path("data", file)
  if (file.exists(local)) return(local)
  local <- file.path("lesson", "data", file)
  if (file.exists(local)) return(local)
  cache <- file.path(tempdir(), "daisy_data")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(cache, file)
  if (!file.exists(dest)) {
    message("downloading ", file, " ...")
    download.file(paste0(DAISY_URL, file), dest, mode = "wb", quiet = TRUE)
  }
  dest
}

## ---------------------------------------------------------------------------
## 3. HELPER 1 - from affiliation (two-mode) data to a one-mode network
## ---------------------------------------------------------------------------
## Almost all innovation network data are *indirectly observed*: we do not see
## the tie, we see two actors sharing an event (a patent, a project, a paper).
## The event x actor incidence matrix B gives the one-mode projection
##      A = t(B) %*% B
## where A[i,j] = number of events shared by actors i and j, and A[i,i] = number
## of events of actor i. Sparse matrices make this cheap even for 10^5 actors.
##
##   dt     : data.table in long format, one row = one actor in one event
##   event  : name of the event column  (patent, project, publication ...)
##   actor  : name of the actor column  (inventor, organisation, institution...)
##   max_size: drop events with more actors than this (huge events create huge
##             cliques: 1 project with 200 partners = 19,900 edges)
proj_two_mode <- function(dt, event, actor, max_size = Inf) {
  d <- unique(as.data.table(dt)[, .(ev = get(event), ac = get(actor))])
  d <- d[!is.na(ev) & !is.na(ac) & ev != "" & ac != ""]
  if (is.finite(max_size)) {
    big <- d[, .N, by = ev][N > max_size, ev]
    if (length(big)) message("dropping ", length(big), " events with > ",
                             max_size, " actors")
    d <- d[!ev %in% big]
  }
  d[, `:=`(ev = as.factor(ev), ac = as.factor(ac))]
  B <- sparseMatrix(i = as.integer(d$ev), j = as.integer(d$ac), x = 1,
                    dims = c(nlevels(d$ev), nlevels(d$ac)),
                    dimnames = list(levels(d$ev), levels(d$ac)))
  A <- Matrix::crossprod(B, B)                 # actor x actor
  n_ev <- diag(A)                              # events per actor
  diag(A) <- 0
  A <- Matrix::drop0(A)
  tri <- Matrix::summary(Matrix::triu(A))      # upper triangle -> edge list
  edges <- data.table(from = colnames(A)[tri$i],
                      to   = colnames(A)[tri$j],
                      weight = tri$x)
  list(edges = edges,
       nodes = data.table(name = colnames(A), n_events = as.numeric(n_ev)),
       incidence = B)
}

## HELPER 2 - assemble an igraph object with node attributes attached
make_net <- function(proj, node_attr = NULL, by = "name") {
  nodes <- proj$nodes
  if (!is.null(node_attr)) {
    node_attr <- copy(as.data.table(node_attr))
    node_attr[, (by) := as.character(get(by))]   # node names are always character
    node_attr <- unique(node_attr, by = by)
    nodes <- merge(nodes, node_attr, by.x = "name", by.y = by, all.x = TRUE)
  }
  graph_from_data_frame(proj$edges, directed = FALSE, vertices = nodes)
}


## ---------------------------------------------------------------------------
## HELPER 4 - Burt's effective size (igraph has constraint(), not this one)
## ---------------------------------------------------------------------------
effective_size <- function(g) {
  A <- as_adjacency_matrix(g, sparse = TRUE); A <- (A > 0) * 1
  deg <- Matrix::rowSums(A)
  redundancy <- Matrix::rowSums((A %*% A) * A)      # 2 x ties among my contacts
  as.numeric(deg - redundancy / pmax(deg, 1))
}

## ---------------------------------------------------------------------------
## HELPER 5 - Gould & Fernandez (1989) brokerage roles
## ---------------------------------------------------------------------------
## v brokers the 2-path i -> v -> j when i and j are NOT directly tied. Given a
## group partition, the role depends on where i, v and j sit:
##   coordinator  i, v, j same group      gatekeeper      i outside, v, j inside
##   representative  i, v inside, j out   consultant      i, j in one other group
##   liaison      i, v, j all different
## In an UNDIRECTED network gatekeeper == representative by construction.
## Block 6 discusses how to read the output; the cost grows with degree^2, so
## filter the network first.
brokerage_roles <- function(g, group) {
  A <- as_adjacency_matrix(g, sparse = TRUE); A <- (A > 0) * 1
  if (!is_directed(g)) A <- ((A + Matrix::t(A)) > 0) * 1
  grp <- factor(group); k <- nlevels(grp); gi <- as.integer(grp); n <- vcount(g)
  nm <- if (is.null(V(g)$name)) as.character(seq_len(n)) else V(g)$name
  out <- matrix(0, n, 5, dimnames = list(nm,
    c("coordinator", "gatekeeper", "representative", "consultant", "liaison")))
  for (v in seq_len(n)) {
    I <- which(A[, v] > 0); O <- which(A[v, ] > 0)
    if (!length(I) || !length(O)) next
    M <- outer(tabulate(gi[I], k), tabulate(gi[O], k))     # all (g_i, g_j) pairs
    both <- intersect(I, O)                                 # drop i == j
    if (length(both)) diag(M) <- diag(M) - tabulate(gi[both], k)
    sub <- A[I, O, drop = FALSE]                            # drop direct i -> j
    if (sum(sub)) {
      GI <- sparseMatrix(seq_along(I), gi[I], dims = c(length(I), k))
      GO <- sparseMatrix(seq_along(O), gi[O], dims = c(length(O), k))
      M <- M - as.matrix(Matrix::t(GI) %*% sub %*% GO)
    }
    gv <- gi[v]; own <- rep(FALSE, k); own[gv] <- TRUE
    other <- M[!own, !own, drop = FALSE]
    out[v, ] <- c(M[gv, gv], sum(M[!own, gv]), sum(M[gv, !own]),
                  sum(diag(other)), sum(other) - sum(diag(other)))
  }
  as.data.table(out, keep.rownames = "name")
}

## ---------------------------------------------------------------------------
## HELPER 6 - economic / knowledge complexity (Hidalgo & Hausmann 2009)
## ---------------------------------------------------------------------------
## Input: a binary ACTOR x CATEGORY matrix M (countries x products, regions x
## technologies, ...). Returns the two complexity indices, i.e. the second
## eigenvectors of the two "method of reflections" operators:
##      Mcc = D^-1 M U^-1 M'      (actor side, ECI/KCI)
##      Mpp = U^-1 M' D^-1 M      (category side, PCI/TCI)
## Signs are the whole difficulty. Conventions used here:
##   - actor index increases with DIVERSITY (diversified actors are complex);
##   - category index is aligned with the actor index (a complex category is one
##     that only complex actors have) - equivalently it DECREASES with ubiquity.
## Getting this backwards silently returns the ranking upside down, which is the
## most common mistake with these measures: always sanity-check the extremes.
complexity <- function(M) {
  d <- rowSums(M); u <- colSums(M)
  M <- M[d > 0, u > 0, drop = FALSE]; d <- rowSums(M); u <- colSums(M)
  ev2 <- function(A) as.numeric(scale(Re(eigen(A)$vectors[, 2])))
  aci <- ev2((M / d) %*% t(sweep(M, 2, u, "/")))          # actor side
  if (cor(aci, d) < 0) aci <- -aci                        # diversified = complex
  cci <- ev2(t(sweep(M, 2, u, "/")) %*% (M / d))          # category side
  avg_actor <- as.numeric(t(M) %*% aci / u)               # mean ACI of holders
  if (cor(cci, avg_actor) < 0) cci <- -cci
  list(actor = setNames(aci, rownames(M)),
       category = setNames(cci, colnames(M)),
       diversity = d, ubiquity = u)
}

## HELPER 7 - the giant (largest) component, we often work on it
giant <- function(g) {
  cmp <- components(g)
  induced_subgraph(g, V(g)[cmp$membership == which.max(cmp$csize)])
}

cat("Setup complete -", R.version.string, "| igraph", as.character(packageVersion("igraph")), "\n")
