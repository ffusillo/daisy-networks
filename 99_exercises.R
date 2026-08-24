#------------------------------------------------------------------------------#
#  DAISY 2026 - Hands-on: Network analysis with innovation data
#  EXERCISES - to try during the last minutes of the session, or afterwards on
#  your own data. Hints are given; the point is that every task is a small
#  variation of the same three-step pipeline:
#     long table  ->  proj_two_mode()  ->  igraph object  ->  measure/aggregate
#------------------------------------------------------------------------------#
source("00_setup.R")

## ---------------------------------------------------------------------------
## 1. HOW FRAGILE ARE THE RANKINGS? (patents)
## ---------------------------------------------------------------------------
## Rebuild the green co-invention network keeping only patents with at most 10
## inventors, and compare the top-10 inventors by betweenness with the full
## network. How many names survive?
## Hint: proj_two_mode(inv, "appln_id", "person_id", max_size = 10)

## ---------------------------------------------------------------------------
## 2. WHO BRIDGES GREEN AND NON-GREEN TECHNOLOGY? (patents)
## ---------------------------------------------------------------------------
## Using pat_all_inventors_ITgreen.csv.gz, classify each inventor as green-only,
## non-green-only or mixed, and test whether "mixed" inventors have higher
## betweenness in the overall network. This is the empirical core of several
## papers on green technology recombination.
## Hint: inv_all[, .(green_share = mean(green)), by = person_id] then merge on
## the vertex table and compare distributions.

## ---------------------------------------------------------------------------
## 3. A NATIONAL SUBNETWORK (CORDIS)
## ---------------------------------------------------------------------------
## Take the Horizon Europe organisation network, extract the subgraph of Italian
## organisations, and find (a) the most central ones, (b) the share of their ties
## that stay inside Italy. Repeat for another country and compare.
## Hint: induced_subgraph(g, V(g)[country == "IT"]) for (a); for (b) work on the
## full edge list and use the country attribute of both endpoints.

## ---------------------------------------------------------------------------
## 4. MONEY AND POSITION (CORDIS)
## ---------------------------------------------------------------------------
## Aggregate the organisation network at country level and check whether
## betweenness in the country network is correlated with the EC contribution
## received per participation. Who punches above its weight?

## ---------------------------------------------------------------------------
## 5. YOUR OWN LITERATURE (OpenAlex)
## ---------------------------------------------------------------------------
## Change the query in 03_publications.R to the topic of your PhD, rebuild the
## institution network and identify the 10 most central institutions. Then look
## at them by hand: do you recognise duplicates or aggregation problems?

## ---------------------------------------------------------------------------
## 6. PERSISTENCE OF REGIONAL KNOWLEDGE STRUCTURES (indicators)
## ---------------------------------------------------------------------------
## Recompute variety, coherence and complexity for the period 2010-2014 and
## correlate them with the 2015-2019 values. Which indicator is most persistent?
## Then regress the growth of green patents 2015-2019 on the 2010-2014
## indicators. (Careful: this is a descriptive exercise, not a causal claim.)

## ---------------------------------------------------------------------------
## 7. STAY IN TWO MODES (advanced)
## ---------------------------------------------------------------------------
## Everything we did projected the two-mode network into one mode. Try instead
## to work directly on the bipartite graph: build it with the incidence matrix
## returned by proj_two_mode() (element $incidence), set the "type" attribute,
## and compute bipartite degree and clustering.
## Hint:
##   B  <- pr$incidence
##   gb <- graph_from_biadjacency_matrix(B)
##   table(V(gb)$type)
## Which measures still make sense? Which ones do not?

## ---------------------------------------------------------------------------
## CHECKLIST - the five questions to ask before you trust a network result
## ---------------------------------------------------------------------------
## 1. NODES     - are the actor identifiers disambiguated? (inventor names,
##                organisation ids, author ids). Errors in nodes are not noise:
##                they systematically split hubs and destroy paths.
## 2. TIES      - what does the tie mean? co-participation is not interaction.
##                Are mega-events (200-partner projects, 68-inventor patents)
##                creating cliques that drive your topology?
## 3. BOUNDARY  - which actors/events are in the population, and why? Country,
##                technology, sector and time-window choices all move the result.
## 4. TIME      - is the network a snapshot, a cumulative window, or a moving
##                window? Centrality is not comparable across window lengths.
## 5. INFERENCE - are you describing or estimating? Network measures are
##                generated regressors: they are endogenous to the same process
##                you are explaining, and they are correlated across units by
##                construction (no independence).
