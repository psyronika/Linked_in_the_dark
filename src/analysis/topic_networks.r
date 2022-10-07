library(data.table)
library(stringr)
library(igraph)
library(Matrix)
library(corpustools)
library(glue)


#make sure to export the graphs before making isolate deletion False
#remember to make a note under the figures that isolates were deleted for the images
#-----------------------------------------------------------------------------------
prepare_data <- function(fname) {
  d = fread(fname, header = T)
  data.table(channel= d$Class,
             topic = d$Name,
             topic_n = d$Topic,
             date =as.Date(d$date),
             words = d$Words)
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


d = prepare_data('df_topics.csv')
d 

#topic overlap over time -using the cosine similarity of top words of each topic
get_topic_network <- function(d, min_date=NULL, max_date=NULL) {
  ds = d
  if (!is.null(min_date)) ds = ds[ds$date >= min_date,]
  if (!is.null(max_date)) ds = ds[ds$date < max_date,]
  m = cross_matrix(ds$words, ds$channel)
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


gt = get_topic_network(d)
gst = backbone_filter(gt, alpha = 0.1, delete_isolates = T) # delete_isolates =F#isolates to False to have measures consider those as well
plot_network(gst, labelsize = 0.9)


## per period
p1 = list(start='2017-03-18', end='2020-03-11')
p2 = list(start='2020-03-12', end='2021-01-19')
p3 = list(start='2021-01-20', end='2021-06-18')

g1t = get_topic_network(d, min_date=p1$start, max_date=p1$end)
g2t = get_topic_network(d, min_date=p2$start, max_date=p2$end)
g3t = get_topic_network(d, min_date=p3$start, max_date=p3$end)

gs1t = backbone_filter(g1t, alpha = 0.1, delete_isolates = T) #only delete isolates for visualization!!!
gs2t = backbone_filter(g2t, alpha = 0.1, delete_isolates = T)
gs3t = backbone_filter(g3t, alpha = 0.1, delete_isolates = T)


#network properties
print(gst) #168 - 306 / 112 (56) - 306 
print(gs1t) #168- 461/ 122 (46) -  461
print(gs2t) #168 - 281/ 112(56) - 281
print(gs3t) #168 -182/ 99 (69) - 182

## network per period
plot_network2(gs1t)
plot_network2(gs2t)
plot_network2(gs3t)

#measure modularity gst 69 clusters /13 clusters (56 isolates)
wt <- cluster_louvain(gst)
members <- membership(wt)
modularity(wt) #0.58
wt




V(gst)[wt$membership == 1]#5 nodes
V(gst)[wt$membership == 2]#4
V(gst)[wt$membership == 3]#21 nodes -52 edges yellow cluster
V(gst)[wt$membership == 4]#9
V(gst)[wt$membership == 5]#14 - green nodes
V(gst)[wt$membership == 6] #19 nodes - 55 edges -  green cluster
V(gsl)[wt$membership == 7] #13 nodes -
V(gsl)[wt$membership == 8]#17 nodes
V(gsl)[wt$membership == 9]#2 nodes


#modularity gt1 54 without /8 clusters (46 isolates)
wt <- cluster_louvain(gs1t)
members <- membership(wt)
modularity(wt) #0.16
wt


#modularity gt2 66 /10 clusters (56 isolates)
wt <- cluster_louvain(gs2t)
members <- membership(wt)
modularity(wt) #0.53
wt


#modularity gt3 82/ 13 clusters (69 isolates)
wt <- cluster_louvain(gs3t)
members <- membership(wt)
modularity(wt) #0.73
wt

#these measures are not incluced by the presence of isolates 

#edge density
edge_density(gt, loops = FALSE)#0.2251
edge_density(g1t, loops = FALSE) #0.2233
edge_density(g2t, loops = FALSE) #0.2251
edge_density(g3t, loops = FALSE) #0.2251


#edge density - to report 
edge_density(gst, loops = FALSE)#0.0492278
edge_density(gs1t, loops = FALSE) #0.06245766
edge_density(gs2t, loops = FALSE) #0.04520592
edge_density(gs3t, loops = FALSE) #0.03751804


#transitivity - 
transitivity(gst,type="global") #0.5468
transitivity(gs1t,type="global") #0.7117
transitivity(gs2t,type="global") #0.49857
transitivity(gs3t,type="global") #0.50107

#assortativity
# positive value means high degree connects with high 
#and low with low
#negative value means that high degree nodes connect to low degree ones
assortativity_degree(gst, directed = FALSE) #0.4533
assortativity_degree(gs1t, directed = FALSE) #0.5286
assortativity_degree(gs2t, directed = FALSE) #0.5683
assortativity_degree(gs3t, directed = FALSE) #0.4568


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
centralize(gst) #0.099/0.1333333
centralize(gs1t) #0.100/0.1213499
centralize(gs2t) #0.083/ 0.1099099
centralize(gs3t) #0.060/ 0.0866821




