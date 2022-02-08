library(data.table)
library(stringr)
library(igraph)
library(Matrix)
library(corpustools)
library(glue)
library(papaja)

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

gu = get_user_network(d)
gsu = backbone_filter(gu, alpha = 0.01, delete_isolates = F) # delete_isolates =F#isolates to False to have measures consider those as well
plot_network2(gsu, labelsize = 0.6)


## per period
p1 = list(start='2017-03-18', end='2020-03-11')
p2 = list(start='2020-03-12', end='2021-01-19')
p3 = list(start='2021-01-20', end='2021-06-18')

g1u = get_user_network(d, min_date=p1$start, max_date=p1$end)
g2u = get_user_network(d, min_date=p2$start, max_date=p2$end)
g3u = get_user_network(d, min_date=p3$start, max_date=p3$end)

gs1u = backbone_filter(g1u, alpha = 0.01, delete_isolates = F) #only delete isolates for visualization!!!
gs2u = backbone_filter(g2u, alpha = 0.01, delete_isolates = F)
gs3u = backbone_filter(g3u, alpha = 0.01, delete_isolates = F)


#network properties
print(gsu) #174 - 287
print(gs1u) #22 - 3
print(gs2u) #99 - 75
print(gs3u) #164 - 205

## network per period
plot_network2(gs1u, labelsize = 0.6)
plot_network2(gs2u, labelsize = 0.6)
plot_network2(gs3u, labelsize = 0.6)



#density
edge_density(gsu, loops = FALSE)#0.0190685
edge_density(gs1u, loops = FALSE) #0.01298701
edge_density(gs2u, loops = FALSE) #0.01546073
edge_density(gs3u, loops = FALSE) #0.01533742


#transitivity
transitivity(gsu,type="global") #0.2757186
transitivity(gs1u,type="global") #NAN
transitivity(gs2u,type="global") #0.2832168
transitivity(gs3u,type="global") # 0.3043478

#assortativity
# positive value means high degree connects with high 
#and low with low
#negative value means that high degree nodes connect to low degree ones
assortativity_degree(gsu, directed = FALSE) #-0.06070338
assortativity_degree(gs1u, directed = FALSE) #nan
assortativity_degree(gs2u, directed = FALSE) #-0.2113975
assortativity_degree(gs3u, directed = FALSE) #0.09758617

#measure modularity gu
wt <- cluster_louvain(gsu)
members <- membership(wt)
modularity(wt) #0.6269263
wt


x <- which.max(sizes(wt))
subg <- induced.subgraph(gsu, which(membership(wt) == x))
print_all(subg)

#modularity gu1
wt <- cluster_louvain(gs1u)
members <- membership(wt)
modularity(wt) #0.1208581
wt

#modularity gu2
wt <- cluster_louvain(gs2u)
members <- membership(wt)
modularity(wt) #0.7518795
#modularity gu3
wt

wt <- cluster_louvain(gs3u)
members <- membership(wt)
modularity(wt) #0.6367708


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

gl = get_link_network(d)
gsl = backbone_filter(gl, alpha = 0.01, delete_isolates = F)
plot_network2(gsl, 0.6) 



## per period
p1 = list(start='2017-03-18', end='2020-03-11')
p2 = list(start='2020-03-12', end='2021-01-19')
p3 = list(start='2021-01-20', end='2021-06-18')

g1l = get_link_network(d, min_date=p1$start, max_date=p1$end)
g2l = get_link_network(d, min_date=p2$start, max_date=p2$end)
g3l = get_link_network(d, min_date=p3$start, max_date=p3$end)

gs1l = backbone_filter(g1l, alpha = 0.01, delete_isolates = F)
gs2l = backbone_filter(g2l, alpha = 0.01, delete_isolates = F)
gs3l = backbone_filter(g3l, alpha = 0.01, delete_isolates = F)

#network properties
print(gsl) #160 -170
print(gs1l)  #20- 1
print(gs2l) #87 - 41
print(gs3l) #148 - 134


## network per period
plot_network2(gs1l, 0.6)
plot_network2(gs2l, 0.6)
plot_network2(gs3l, 0.6)



#edge density
edge_density(gl, loops = FALSE)#0.01336478
edge_density(g1l, loops = FALSE) #0.005263158
edge_density(g2l, loops = FALSE) #0.01095964
edge_density(g3l, loops = FALSE) # 0.01231844

#transitivity
transitivity(gsl,type="global") #0.3642565
transitivity(gs1l,type="global") #NAN
transitivity(gs2l,type="global") #0.4210526
transitivity(gs3l,type="global") #0.4704684

#assortativity
# positive value means high degree connects with high 
#and low with low
#negative value means that high degree nodes connect to low degree ones
assortativity_degree(gsl, directed = FALSE) #0.1530868
assortativity_degree(gs1l, directed = FALSE) #nan
assortativity_degree(gs2l, directed = FALSE) #0.2556282
assortativity_degree(gs3l, directed = FALSE) #0.2133746


#measure modularity gl
wt <- cluster_louvain(gsl)
members <- membership(wt)
modularity(wt) #0.7034536

wt
#largest community
x <- which.max(sizes(wt))
subg <- induced.subgraph(gsl, which(membership(wt) == x))
print_all(subg)

V(gsl)[wt$membership == 1]
V(gsl)[wt$membership == 2]
V(gsl)[wt$membership == 3]
V(gsl)[wt$membership == 4]


community <- wt$membership
v <- as_ids(V(gsl))
vcom <- data.frame(v,community)
vcom

#modularity gl1
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


#modularity gl2
wt <- cluster_louvain(gs2l)
members <- membership(wt)
modularity(wt) #0.8456044
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

#modularity gl3
wt <- cluster_louvain(gs3l)
members <- membership(wt)
modularity(wt) #0.7709169
print(wt)

#degree
ddl<- degree_distribution(gsl)
ddu<- degree_distribution(gsu)
plot(ddu)
plot(ddl)



#properties of graphs before edge filtering
print(gl) #nodes = 160, edges = 4339
print(gu)#nodes = 174, edges = 3499


#properties after edge filtering alpha = 0.01, no isolate deletion
print(gsl) #nodes = 160, edges = 170
print(gsu) #nodes = 174, edges  = 287

#properties after edge filtering alpha = 0.01, after isolate deletion
print(gsl) #nodes = 118, edges = 170
print(gsu) #nodes = 89, edges  = 287



#clusters per user overlap and link overlap
#link
#components(gl, mode = c("weak", "strong"))#146 -  largest: 15 
components(gsl, mode = c("weak", "strong")) #bone 56 - 84
#components(g1l, mode = c("weak", "strong")) #11  -  8
components(gs1l, mode = c("weak", "strong")) #bones 19 - 2
#components(g2l, mode = c("weak", "strong")) #11  -  77
components(gs2l, mode = c("weak", "strong")) #bones 52 - 10
#components(g3l, mode = c("weak", "strong")) #11  -  138
components(gs3l, mode = c("weak", "strong")) #bones 61- 57 
#user
#components(gu, mode = c("weak", "strong")) #78 - largest:97
components(gsu, mode = c("weak", "strong"))# bone 87 - 87 
#components(g1u, mode = c("weak", "strong")) #15  -  8
components(gs1u, mode = c("weak", "strong")) #bones 19 - 2
#components(g2u, mode = c("weak", "strong")) #49  -  51
components(gs2u, mode = c("weak", "strong")) #bones 55 - 40
#components(g3u, mode = c("weak", "strong")) #71  -  94
components(gs3u, mode = c("weak", "strong")) #bones 79 -  84 


#distance
mean_distance(gsu, directed = FALSE, unconnected = TRUE) #2.840994
mean_distance(gsl, directed = FALSE, unconnected = TRUE) #4.918159

#look more into this
diversity(gsu, weights = NULL, vids = V(gsu))




glue("Average path length in  link network: ", mean_distance(gsl, directed = F))
glue("Distance between DeBataafseRepubliek and FVDNL in network:", distances(gsl, v="DeBataafseRepubliek", to="FVDNL", weights = NULL))

glue("Average path length in  user network: ", mean_distance(gsu, direct = F))
glue("Distance between DeBataafseRepubliek and FVDNL in network:", distances(gsu, v="DeBataafseRepubliek", to="FVDNL", weights = NULL))



glue("Neighbors of FVDNL in Link Network:")
neighbors(gsl, V(gsl)["FVDNL"], mode="all")
neighbors(gsl, V(gsl)["DeBataafseRepubliek"], mode="all")
neighbors(gsl, V(gsl)["leefbewust"], mode="all") #CLOSEST NEIGHBOR IS a wappie channel
neighbors(gsl, V(gsl)["avondkloknl"], mode="all")
neighbors(gsl, V(gsl)["blckbxtv"], mode="all")
neighbors(gsl, V(gsl)["DeDagelijkseStandaard"], mode="all") #alternative news clisest to wilderspvv 
neighbors(gsl, V(gsl)["verkiezingen2021"], mode="all") #channel on elections closest neighbor is a red pill channel


#difference between graphs over time
#plots the difference between two graphs
#make them with isolate deletion though

dif1 <- difference(gs2l,gs1l, byname = "auto")
plot_network2(dif1, 0.6)
dif1
dif2 <- difference(gs3l,gs2l, byname = "auto")
plot_network2(dif2, 0.6)
dif2
dif3 <- difference(gs2u,gs1u, byname = "auto")
plot_network2(dif3, 0.6)
dif3
dif4 <- difference(gs3u,gs2u, byname = "auto")
plot_network2(dif4, 0.6)
dif4

#between users and links graphs
dif5 <- difference(gsu,gsl, byname = "auto")
plot_network2(dif5, 0.6)
dif5

# Density measure
percent_density <- function(x) {
  (ecount(x)/((vcount(x)*(vcount(x)-1)/2)))*100
}
network_statistics <- function(x) {
  print(names(x))
  print('Node count')
  print(vcount(x))
  print('Edge count')
  print(ecount(x))
  density <- percent_density(x)
  centralization <- centr_degree(x,
                                 mode = "all",
                                 loops = FALSE, 
                                 normalized = TRUE)$centralization
  print('normalized centralization =')
  print(centralization)
  print('density =')
  print(density)
  print('assortativity degree =')
  print(assortativity_degree(x, directed = FALSE))
  
}


#network measures
network_statistics(gsu)
network_statistics(gs1u)
network_statistics(gs2u)
network_statistics(gs3u)

network_statistics(gsl)
network_statistics(gs1l)
network_statistics(gs2l)
network_statistics(gs3l)




dif1 <- gs2l %m% gs1l
print_all(dif1)
dif2 <- gs2u %m% gs1u
print_all(dif2)

dif3 <- gs3l %m% gs2l
print_all(dif3)
dif4 <- gs3u %m% gs2u
print_all(dif4)



