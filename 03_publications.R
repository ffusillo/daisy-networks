#------------------------------------------------------------------------------#
#  DAISY 2026 - Hands-on: Network analysis with innovation data
#  BLOCK 3 of the session (~8 min) - PUBLICATION DATA: the OpenAlex API
#
#  Publications get a session of their own in this school, so here we only look
#  at (i) how to pull relational data out of the OpenAlex API in three lines,
#  (ii) how the very same projection logic gives co-authorship networks, and
#  (iii) what changes when the actors are institutions rather than people.
#
#  OpenAlex: fully open (CC0), no key, ~250M works. Be polite: add your e-mail
#  ("polite pool"), max 10 requests/second, 100k/day.
#  R packages worth knowing: openalexR (wrapper), rcrossref, europepmc,
#  bibliometrix. Here we use the raw API so you see what happens.
#------------------------------------------------------------------------------#
source("00_setup.R")

MAIL <- "your.name@your.university.it"      # <- put YOUR address here
oa <- function(path, ...) {
  url <- paste0("https://api.openalex.org/", path, "&mailto=", MAIL)
  fromJSON(URLencode(url), simplifyVector = FALSE)
}

## ===========================================================================
## 1. Aggregate queries: indicators without downloading any record
## ===========================================================================
## group_by returns counts, not records: perfect for descriptive statistics.
res <- oa(paste0("works?filter=title_and_abstract.search:circular economy,",
                 "publication_year:2015-2024&group_by=publication_year"))
by_year <- rbindlist(lapply(res$group_by, function(x)
  data.table(year = as.integer(x$key), works = x$count)))[order(year)]
by_year

## Which countries publish on it?
res <- oa(paste0("works?filter=title_and_abstract.search:circular economy,",
                 "publication_year:2020-2024&group_by=authorships.countries"))
rbindlist(lapply(res$group_by, function(x)
  data.table(country = x$key_display_name, works = x$count)))[1:15]

## ===========================================================================
## 2. Record-level download: works with their authorships
## ===========================================================================
## 200 records per page, cursor paging. select= keeps the payload small.
fetch_works <- function(n_pages = 2) {
  q <- paste0("works?filter=title_and_abstract.search:circular economy,",
              "publication_year:2020-2024,type:article,has_orcid:true",
              "&select=id,display_name,publication_year,cited_by_count,authorships",
              "&per-page=200")
  out <- list(); cursor <- "*"
  for (i in seq_len(n_pages)) {
    r <- oa(paste0(q, "&cursor=", cursor))
    out[[i]] <- r$results; cursor <- r$meta$next_cursor
    if (is.null(cursor)) break
  }
  unlist(out, recursive = FALSE)
}

## flatten works x authors x institutions into a long table
flatten_authorships <- function(works) rbindlist(lapply(works, function(w)
  rbindlist(lapply(w$authorships, function(a) {
    ins <- if (length(a$institutions)) a$institutions else list(list())
    rbindlist(lapply(ins, function(s) data.table(
      work_id      = sub(".*/", "", w$id),
      year         = w$publication_year,
      cited_by     = w$cited_by_count,
      author_id    = sub(".*/", "", a$author$id %||% NA_character_),
      author_name  = a$author$display_name %||% NA_character_,
      inst_id      = sub(".*/", "", s$id %||% NA_character_),
      inst_name    = s$display_name %||% NA_character_,
      inst_country = s$country_code %||% NA_character_,
      inst_type    = s$type %||% NA_character_)))
  }), fill = TRUE)), fill = TRUE)

## live if the wifi cooperates, otherwise the cached extract (2,000 works)
aut <- tryCatch(flatten_authorships(fetch_works(2)),
                error = function(e) {
                  message("API unreachable, using the cached extract")
                  fread(daisy_data("openalex_ce_authorships.csv.gz"))
                })
if (nrow(aut) < 1000) aut <- fread(daisy_data("openalex_ce_authorships.csv.gz"))

aut[, .(works = uniqueN(work_id), authors = uniqueN(author_id),
        institutions = uniqueN(inst_id), countries = uniqueN(inst_country))]

## ===========================================================================
## 3. Three networks out of one table (same helper as blocks 1 and 2)
## ===========================================================================
## (a) co-authorship between researchers
aut_attr <- unique(aut[, .(author_id, author_name, country = inst_country)],
                   by = "author_id")
g_aut <- make_net(proj_two_mode(aut, "work_id", "author_id", max_size = 30),
                  node_attr = aut_attr, by = "author_id")
cmp <- components(g_aut)
c(authors = vcount(g_aut), ties = ecount(g_aut),
  giant_share = max(cmp$csize) / vcount(g_aut))

## (b) collaboration between institutions - the level most used in economics
g_ins <- make_net(proj_two_mode(aut, "work_id", "inst_id", max_size = 30),
                  node_attr = unique(aut[!is.na(inst_id),
                                         .(inst_id, inst_name, inst_country, inst_type)]),
                  by = "inst_id")
g_ins
V(g_ins)$degree <- degree(g_ins)
V(g_ins)$betw   <- betweenness(g_ins, weights = NA, normalized = TRUE)
ins <- as.data.table(as_data_frame(g_ins, what = "vertices"))
ins[order(-degree)][1:12, .(inst_name, inst_country, inst_type,
                            papers = n_events, degree, betw = round(betw, 3))]

## (c) country co-publication network
g_ctry <- make_net(proj_two_mode(unique(aut[, .(work_id, inst_country)]),
                                 "work_id", "inst_country"))
sort(strength(g_ctry), decreasing = TRUE)[1:12]

el <- as.data.table(as_data_frame(g_ctry, what = "edges"))
el[order(-weight)][1:10]

ggraph(delete_edges(g_ctry, E(g_ctry)[weight < 5]), layout = "stress") +
  geom_edge_link0(aes(edge_width = weight), edge_colour = "grey80") +
  scale_edge_width(range = c(0.2, 3)) +
  geom_node_point(aes(size = n_events), fill = "#d95f0e", shape = 21, colour = "white") +
  geom_node_text(aes(label = name), size = 3, repel = TRUE) +
  theme_graph(base_family = "sans") + theme(legend.position = "none") +
  labs(title = "Country co-publication network, circular economy research")

## ===========================================================================
## 4. What to keep in mind (and what connects this to the other blocks)
## ===========================================================================
## - Author disambiguation: OpenAlex ids are algorithmic. Same problem as
##   REGPAT person_id and CORDIS organisationID: measurement error in the NODES
##   propagates to every network statistic. Check your top-degree actors by hand.
## - Coverage/selection: our query is a keyword search; a different query is a
##   different network. Prefer topic/concept ids or a validated keyword list, and
##   always report the query in the paper.
## - Linking science and technology: patent front-page and non-patent-literature
##   citations (PATSTAT TLS214, Lens.org, Reliance-on-Science) let you build
##   *directed* science -> technology networks. That is where publication and
##   patent data meet.
