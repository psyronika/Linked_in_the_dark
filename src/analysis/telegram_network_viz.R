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
             user = as.numeric(str_extract(d$user_data, '(?<=\\(id=)[0-9]+')),
             date =as.Date(d$date),
             link = d$links)
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


d = prepare_data('clean_df.csv')
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

plot_network <- function(gs) {
  V(gs)$size = V(gs)$size*0.7
  V(gs)$cluster = igraph::cluster_walktrap(gs)$membership
  colors = rainbow(max(V(gs)$cluster))
  V(gs)$color = colors[V(gs)$cluster]
  par(mar=c(0,0,0,0))
  plot(gs, labelsize_coef = 0.7, reduce_labeloverlap = 2, edge.arrow.size=0.1)
}

plot_network2 <- function(gs, labelsize=1) {
  ## use plot_semnet with custom clustering
  
  ## cluster data.  
  cl = igraph::cluster_louvain(gs)   
  ## Set cluster membership as vertex attribute.
  V(gs)$cluster = cl$membership
  ## Use attribute in plot for colors
  plot_semnet(gs, labelsize_coef = labelsize, vertexcolor_attr = 'cluster', reduce_labeloverlap = T)
}

gub = get_user_network(d)
gsub = backbone_filter(gub, alpha = 0.01) # delete_isolates =F#isolates to False to have measures consider those as well
plot_network2(gsub, labelsize = 0.6)

print(gsub)

## per period
p1 = list(start='2017-03-18', end='2020-03-11')
p2 = list(start='2020-03-12', end='2021-01-19')
p3 = list(start='2021-01-20', end='2021-06-18')

g1ub = get_user_network(d, min_date=p1$start, max_date=p1$end)
g2ub = get_user_network(d, min_date=p2$start, max_date=p2$end)
g3ub = get_user_network(d, min_date=p3$start, max_date=p3$end)

gs1ub = backbone_filter(g1ub, alpha = 0.01) #only delete isolates for visualization!!!
gs2ub = backbone_filter(g2ub, alpha = 0.01)
gs3ub = backbone_filter(g3ub, alpha = 0.01)

## network per period
plot_network2(gs1ub, labelsize = 0.6)
plot_network2(gs2ub, labelsize = 0.6)
plot_network2(gs3ub, labelsize = 0.6)


#network properties
print(gs1ub) 
print(gs2ub)
print(gs3ub)


#---------------------------------------------------------------------
#link overlap over time

get_link_network <- function(d, min_date=NULL, max_date=NULL) {
  ds = d[!is.na(d$link) & !d$link == '',]
  if (!is.null(min_date)) ds = ds[ds$date >= min_date,]
  if (!is.null(max_date)) ds = ds[ds$date < max_date,]
  m = cross_matrix(ds$link, ds$channel)
  g = graph_from_adjacency_matrix(cosine_sim(m), mode = 'undirected', weighted = T, diag = F)
  V(g)$size = log(1+colSums(m))
  g
}

glb = get_link_network(d)
gslb = backbone_filter(glb, alpha = 0.01)
plot_network2(gslb, 0.6) 


## per period
p1 = list(start='2017-03-18', end='2020-03-11')
p2 = list(start='2020-03-12', end='2021-01-19')
p3 = list(start='2021-01-20', end='2021-06-18')

g1lb = get_link_network(d, min_date=p1$start, max_date=p1$end)
g2lb = get_link_network(d, min_date=p2$start, max_date=p2$end)
g3lb = get_link_network(d, min_date=p3$start, max_date=p3$end)

gs1lb = backbone_filter(g1lb, alpha = 0.01)
gs2lb = backbone_filter(g2lb, alpha = 0.01)
gs3lb = backbone_filter(g3lb, alpha = 0.01)

## network per period
plot_network2(gs1lb, 0.6)
plot_network2(gs2lb, 0.6)
plot_network2(gs3lb, 0.6)


#difference between graphs over time
#plots the difference between two graphs
#make them with isolate deletion though

dif1 <- difference(gs2lb,gs1lb, byname = "auto")
plot_network2(dif1, 0.6)
dif1
dif2 <- difference(gs3lb,gs2lb, byname = "auto")
plot_network2(dif2, 0.6)
dif2
dif3 <- difference(gs2ub,gs1ub, byname = "auto")
plot_network2(dif3, 0.6)
dif3
dif4 <- difference(gs3ub,gs2ub, byname = "auto")
plot_network2(dif4, 0.6)
dif4

#between users and links graphs
dif5 <- difference(gsub,gslb, byname = "auto")
plot_network2(dif5, 0.6)
dif5

print(dif1)
print(dif2)
print(dif3)
print(dif4)
print(dif5)