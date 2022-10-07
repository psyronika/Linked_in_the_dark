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
  data.table(channel= d$source,
             user = as.numeric(d$ids, '(?<=\\(id=)[0-9]+'),
             date =as.Date(d$date),
             link = d$expanded_url)
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



get_link_network <- function(d, min_date=NULL, max_date=NULL) {
  ds = d[!is.na(d$link) & !d$link == '',]
  if (!is.null(min_date)) ds = ds[ds$date >= min_date,]
  if (!is.null(max_date)) ds = ds[ds$date < max_date,]
  m = cross_matrix(ds$link, ds$channel)
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
  plot_semnet(gs, labelsize_coef = labelsize, vertexcolor_attr = 'cluster', reduce_labeloverlap = 10)
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

d = prepare_data('df_urls.csv')
d 


gl = get_link_network(d)
gsl = backbone_filter(gl, alpha = 0.01, delete_isolates = T)
plot_network(gsl, labelsize = 0.9) 


## per period
p1 = list(start='2017-03-18', end='2020-03-11')
p2 = list(start='2020-03-12', end='2021-01-19')
p3 = list(start='2021-01-20', end='2021-06-18')

g1l = get_link_network(d, min_date=p1$start, max_date=p1$end)
g2l = get_link_network(d, min_date=p2$start, max_date=p2$end)
g3l = get_link_network(d, min_date=p3$start, max_date=p3$end)

gs1l = backbone_filter(g1l, alpha = 0.01, delete_isolates = T)
gs2l = backbone_filter(g2l, alpha = 0.01, delete_isolates = T)
gs3l = backbone_filter(g3l, alpha = 0.01, delete_isolates = T)

#network properties
print(gsl) #160 -174/ 118 - 174
print(gs1l)  #20 - 1/2-1
print(gs2l) #87 - 44/52 - 44
print(gs3l) #148 -139/104 - 139


## network per period
plot_network2(gs1l)
plot_network2(gs2l)
plot_network2(gs3l)


#measure modularity gl 61 clusters/19
wt <- cluster_louvain(gsl)
members <- membership(wt)
modularity(wt) #0.71
wt
#largest community
x <- which.max(sizes(wt))
subg <- induced.subgraph(gsl, which(membership(wt) == x))
print(subg)

V(gsl)[wt$membership == 1]#12 nodes
V(gsl)[wt$membership == 2]
V(gsl)[wt$membership == 3]
V(gsl)[wt$membership == 4]
V(gsl)[wt$membership == 5]#8 nodes
E(gsl)[wt$membership == 6] #17 nodes - green 26 edges 
V(gsl)[wt$membership == 7] #11 nodes - light green
V(gsl)[wt$membership == 8]#2 nodes
V(gsl)[wt$membership == 9]#10 nodes




community <- wt$membership
v <- as_ids(V(gsl))
vcom <- data.frame(v,community)
vcom

#modularity gsl1 #19 clusters/1 (18 isolates)
wt <- cluster_louvain(gs1l)
members <- membership(wt)
modularity(wt) #.0
wt

V(gs1l)[wt$membership == 1]
V(gs1l)[wt$membership == 2]
V(gs1l)[wt$membership == 3]
V(gs1l)[wt$membership == 4]


community <- wt$membership
v <- as_ids(V(gs1l))
vcom <- data.frame(v,community)
vcom


#modularity gsl2- 51/16 (35 isolates)
wt <- cluster_louvain(gs2l)
members <- membership(wt)
modularity(wt) #0.85
wt


V(gs2l)[wt$membership == 1]
V(gs2l)[wt$membership == 2]
V(gs2l)[wt$membership == 3]
V(gs2l)[wt$membership == 4]


community <- wt$membership
v <- as_ids(V(gs2l))
vcom <- data.frame(v,community)
vcom

x <- which.max(sizes(wt))
subg <- induced.subgraph(gs2l, which(membership(wt) == x))
subg

#modularity gsl3 #64 /20 (44 isolates)
wt <- cluster_louvain(gs3l)
members <- membership(wt)
modularity(wt) #0.76
print(wt)


#edge density
#edge_density(gl, loops = FALSE)#0.1152
#edge_density(g1l, loops = FALSE) #0.0315
#edge_density(g2l, loops = FALSE) #0.0829
#edge_density(g3l, loops = FALSE) #0.1196


#edge density - to report - noticed some mistakes in paper
edge_density(gsl, loops = FALSE)#0.01367(1.4%)
edge_density(gs1l, loops = FALSE) #0.00526(0.5%)
edge_density(gs2l, loops = FALSE) #0.01176(1.2%)
edge_density(gs3l, loops = FALSE) #0.01277(1.3%)


#transitivity - 
transitivity(gsl,type="global") #0.358
transitivity(gs1l,type="global") #NAN
transitivity(gs2l,type="global") #0.4218
transitivity(gs3l,type="global") #0.4646

#assortativity
# positive value means high degree connects with high 
#and low with low
#negative value means that high degree nodes connect to low degree ones
assortativity_degree(gsl, directed = FALSE) # 0.1442661
assortativity_degree(gs1l, directed = FALSE) #nan
assortativity_degree(gs2l, directed = FALSE) #0.2693358
assortativity_degree(gs3l, directed = FALSE) #0.2242294


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

#network measures (centralization)
centralize(gsl) #0.107/0.139552
centralize(gs1l) #0.0526/NAN
centralize(gs2l) #0.088/0.08784314
centralize(gs3l) #0.132/0.1319246
