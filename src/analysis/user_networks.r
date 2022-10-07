library(data.table)
library(stringr)
library(igraph)
library(Matrix)
library(corpustools)
library(glue)

#make sure to export the graphs before making isolate deletion False
#remember to make a note under the figures that isolates were deleted for the images
#-----------------------------------
prepare_data <- function(fname) {
  d = fread(fname, header = T)
  data.table(channel= gsub('.*/', '', d$url),
             user = d$ids,
             date =as.Date(d$date))
}

create_matrix <- function(rows, columns, values) {
  urows = unique(rows)
  ucolumns = unique(columns)
  m = spMatrix(length(urows), length(ucolumns), 
               match(rows, urows), match(columns, ucolumns), values)
  dimnames(m) = list(urows, ucolumns)
  m
}

cross_matrix <- function(rows, columns) {
  d = data.table(rows=rows, columns=columns)
  d = d[,list(n=.N), by=c('rows','columns')]
  create_matrix(d$rows, d$columns, d$n)
}

cosine_sim <- function(m) {
  mag = sqrt(colSums(m^2))
  m2 = t(t(m) / mag)
  crossprod(m2)
}

conprob <- function(m) {
  m[m>0] = 1 ## dichotomizes m, so it's like jaccard
  cp = crossprod(m)
  cp / colSums(m)
}


d = prepare_data('df_users.csv')
d 

#user overlap over time - only works with cosine similarity
get_user_network <- function(d, min_date=NULL, max_date=NULL) {
  ds = d
  if (!is.null(min_date)) ds = ds[ds$date >= min_date,]
  if (!is.null(max_date)) ds = ds[ds$date < max_date,]
  m = cross_matrix(ds$user, ds$channel)
  g = graph_from_adjacency_matrix(cosine_sim(m), mode = 'undirected', weighted = T, diag = F)
  V(g)$size = log(1+colSums(m))
  g
}


plot_network <- function(gs, labelsize=1) {
  ## use plot_semnet with custom clustering
  
  ## cluster data.  
  cl = igraph::cluster_louvain(gs)   
  ## Set cluster membership as vertex attribute.
  V(gs)$cluster = cl$membership
  ## Use attribute in plot for colors
  plot_semnet(gs, labelsize_coef = labelsize, vertexcolor_attr = 'cluster', reduce_labeloverlap = T)
}


plot_network2 <- function(gs) {
  ## use plot_semnet with custom clustering
  ## cluster data.  
  cl = igraph::cluster_louvain(gs)   
  ## Set cluster membership as vertex attribute.
  V(gs)$cluster = cl$membership
  ## Use attribute in plot for colors
  plot_semnet(gs, labelsize_coef = NaN, vertexcolor_attr = 'cluster', vertex.label=NA)
}


gu = get_user_network(d)
gsu = backbone_filter(gu, alpha = 0.01, delete_isolates = T) # delete_isolates =F#isolates to False to have measures consider those as well
plot_network(gsu, labelsize = 0.9)


## per period
p1 = list(start='2017-03-18', end='2020-03-11')
p2 = list(start='2020-03-12', end='2021-01-19')
p3 = list(start='2021-01-20', end='2021-06-18')

g1u = get_user_network(d, min_date=p1$start, max_date=p1$end)
g2u = get_user_network(d, min_date=p2$start, max_date=p2$end)
g3u = get_user_network(d, min_date=p3$start, max_date=p3$end)

gs1u = backbone_filter(g1u, alpha = 0.01, delete_isolates = T) 
gs2u = backbone_filter(g2u, alpha = 0.01, delete_isolates = T)
gs3u = backbone_filter(g3u, alpha = 0.01, delete_isolates = T) 

## network per period
plot_network2(gs1u)
plot_network2(gs2u)
plot_network2(gs3u)

#network properties
print(gu) #174 nodes - 3499 edges
print(gsu) #174 nodes - 287 edges/ 89 - 287
print(gs1u) #22 nodes - 3 edges/ 6 - 3
print(gs2u) #99 nodes - 75 edges/48 - 75
print(gs3u) #164 nodes - 205 edges/88 - 205

#measure modularity gu 96 clusters/12 (84 isolates)
wt <- cluster_louvain(gsu)
members <- membership(wt)
modularity(wt) #0.63
wt


#modularity gs1u 19 clusters/3 (16 isolates)
wt <- cluster_louvain(gs1u)
members <- membership(wt)
modularity(wt) #0.12
wt


#modularity gs2u 58 clusters/7 (51 isolates)
wt <- cluster_louvain(gs2u)
members <- membership(wt)
modularity(wt) #0.75 
wt


#modularity gs3u 87 clusters /12 (75 isolates)
wt <- cluster_louvain(gs3u)
members <- membership(wt)
modularity(wt) #0.64
wt

#----------These measure need to be corrected in paper

#edge density
#edge_density(gu, loops = FALSE)#0.1152
#edge_density(g1u, loops = FALSE) #0.0315
#edge_density(g2u, loops = FALSE) #0.0829
#edge_density(g3u, loops = FALSE) #0.1196


#edge density
edge_density(gsu, loops = FALSE)#0.019
edge_density(gs1u, loops = FALSE) #0.013
edge_density(gs2u, loops = FALSE) #0.016
edge_density(gs3u, loops = FALSE) #0.015


#transitivity - 
transitivity(gsu,type="global") #0.2757
transitivity(gs1u,type="global") #NAN
transitivity(gs2u,type="global") #0.2832168
transitivity(gs3u,type="global") #0.3043

#assortativity
# positive value means high degree connects with high and low with low
#negative value means that high degree nodes connect to low degree ones
assortativity_degree(gsu, directed = FALSE) # - 0.0607
assortativity_degree(gs1u, directed = FALSE) #nan
assortativity_degree(gs2u, directed = FALSE) # -0.2113
assortativity_degree(gs3u, directed = FALSE) #0.09758


#centralization
centralize <- function(x) {
  print(names(x))
  centralization <- centr_degree(x,
                                 mode = "all",
                                 loops = FALSE, 
                                 normalized = TRUE)$centralization
  print('normalized centralization =')
  print(centralization)
  
}

#network measures (centralization) - sensitive to isolates removal
centralize(gsu) #0.2505225
centralize(gs1u) #0.00
centralize(gs2u) #0.1970398
centralize(gs3u) #0.1804