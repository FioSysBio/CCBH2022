############################################################################
#
# Study of the basic topology of the Pseudomonas aeruginosa Gene Regulation Network
#
############################################################################


library(dplyr)
library(tibble)
library(readr)
library(igraph)
library(scales)


dados <- 
  read_csv2("GRN.csv")


# Data pre-processing
c1 <- 
  dados$`Regulator (TF or sigma)` %>%
  strsplit(" ") %>%
  unlist()

c1.TF <- c1[gtools::odd(1:length(c1))]
c1.ortologo.TF <- c1[gtools::even(1:length(c1))]

dados$`Regulator (TF)` <- c1.TF
dados$`orthologs of TF` <- c1.ortologo.TF


nrow(dados)  ==
  sum(paste(dados$`Regulator (TF)`,dados$`orthologs of TF`) == dados$`Regulator (TF or sigma)`)

dados <- 
  dados %>%
  select(`Regulator (TF)`,
         `Target gene`,
         `mode of regulation`, 
         `orthologs of TF`,
         `Ortholog of the target gene`
  )

dados$`orthologs of TF` <- gsub("\\(|\\)", "",
                                dados$`orthologs of TF`)

rm(c1, c1.ortologo.TF,c1.TF)

auxTF <- dados[,c(1,4)] %>% setNames(c("gene_CCBH4851", "orthologs"))
auxTarget <- dados[,c(2,5)] %>% setNames(c("gene_CCBH4851", "orthologs"))

vert <- 
  dplyr::union(auxTF,auxTarget) %>%
  filter(!is.na(gene_CCBH4851) )

vert$rotulo <- paste0(vert$gene_CCBH4851,
                      " (",
                      vert$orthologs, 
                      ")")      

arestas <- 
  dados[,c(1,2,3)] %>%
  filter(!is.na(`Target gene`))


##### Finalmente: REDE COM IGRAPH

Rede <- graph_from_data_frame(d = arestas,
                              directed = TRUE
                              # vertices = vert
)


# Exploring the Network

V(Rede)$color <- ifelse(V(Rede)$name %in% auxTF$gene_CCBH4851,
                        "yellow", "lightblue")

codificacao <- '"+" = "green" ; "-" = "red" ; "?" = "black" ; "d" = "darkgrey"'
E(Rede)$color <- car::Recode(E(Rede)$`mode of regulation`,
                             codificacao)


V(Rede)$size <- 2+ log(1+degree(Rede, mode = "out"))

V(Rede)$name <- V(Rede)$rotulo 


# Topology

V(Rede)$name <- vert$gene_CCBH4851 


nrow(vert)

nrow(arestas)

sum(V(Rede)$color == "yellow")

sum(V(Rede)$color == "lightblue")

scientific(graph.density(Rede, loops=TRUE))

###
p.kin <- degree_distribution(Rede, mode="in")
p.kin.na <- ifelse(p.kin == 0,NA,p.kin)
min.kin <- min(degree(Rede, mode="in"))
max.kin <-  max(degree(Rede, mode="in"))
plot(min.kin:max.kin, p.kin.na,
     xlab= "k-in (grau input)", ylab= "P(k-in)", type="h")
title(sub="Figure 3: Distribution of Degree k-in", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")

### Power Law ajusts kin

x.in <- log10(min.kin:max.kin)
y.in <- log10(p.kin.na)
data.in <- 
  data.frame(X=x.in,Y=y.in) %>%
  filter(!is.na(X) & !is.na(Y) & X != -Inf)
ajuste.in <- lm(Y~X,data=data.in)
log.A.in <- ajuste.in$coefficients[1]
A.in <- 10^(log.A.in) 
gama.in <- -ajuste.in$coefficients[2]  
#
plot(x.in, y.in,
     xlab= "log(k-in)", ylab= "log P(k-in)")
plotrix::ablineclip(log.A.in, -gama.in, x1= 0,x2=log10(max.kin))

###
# output
p.kout <- degree_distribution(Rede, mode="out")
p.kout.na <- ifelse(p.kout == 0,NA,p.kout)
min.kout <- min(degree(Rede, mode="out"))
max.kout <- max(degree(Rede, mode="out"))
#
plot(min.kout:max.kout, p.kout.na, 
     xlab= "k-out (grau output)", ylab= "P(k-out)", type="h")
title(sub="Figure 3: Distribution of degree k-out", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")
#
# Axis y in log10
plot(min.kout:max.kout, p.kout.na, 
     xlab= "k-out (Degree output)", ylab= "P(k-out)", type="h",log="y")
title(sub="Figure 4: Distribution of Degree k-out", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")


### Power Law ajusts kout
x.out <- log10(min.kout:max.kout)
y.out <- log10(p.kout.na)
data.out <- 
  data.frame(X=x.out,Y=y.out) %>%
  filter(!is.na(X) & !is.na(Y) & X != -Inf)
ajuste.out <- lm(Y~X,data=data.out)
log.A.out <- ajuste.out$coefficients[1]
A.out <- 10^(log.A.out) 
gama.out <- -ajuste.out$coefficients[2]  
#
plot(x.out, y.out,
     xlab= "log(k-out)", ylab= "log P(k-out)")
title(sub="Figura 5: Distribuição de Grau k-out log-log", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")
plotrix::ablineclip(log.A.out, -gama.out, x1= 0,x2=log10(max.kout))



# clustering coefficient


# global
CoeffCluster.global <- scientific(transitivity(Rede,type="globalundirected"))

# mean
CoeffCluster.medio <- scientific(transitivity(Rede,type="average"))

# local
CoeffCluster.i <- transitivity(Rede,type="localundirected",
                               isolates = "NaN",
                               vids = NULL,
                               weights = NULL)

# total distribution
hist(CoeffCluster.i,
     xlab= "coeficiente de clusterização local", ylab= "frequência", main=NULL)
title(sub="Figura 6: Distribuição Total de Coef. Clusterização", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")

# proportion of vertices with zero coef
propCzero <-table(CoeffCluster.i)[1]/nrow(vert)


# proportion of vertices with coef 1:
propChum <-table(CoeffCluster.i)[nrow(table(CoeffCluster.i))]/nrow(vert)


# distribution without nulls and without 1
hist(ifelse(CoeffCluster.i ==0 | CoeffCluster.i ==1, NA, CoeffCluster.i),
     xlab= "coeficiente de clusterização local", ylab= "frequência", main=NULL)
title(sub="Figura 2: Distribuição Parcial de Coef. Clusterização", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")

################
# Scatter plot between clustering coefficient by degree
k.i <- degree(Rede,mode="all")

C.k.i <- CoeffCluster.i 
#
plot(k.i, C.k.i,
     xlab="k (grau total)", ylab= "C(k) (coef. cluster. por grau k)")
title(sub="Figura 7: Coef. Clusterização em função do grau k dos vértices ", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")

# k =0 ou k=1 absent 
C.k.i.filtrado <- transitivity(Rede,type="localundirected",
                               isolates = "NaN")
#
plot(k.i, C.k.i.filtrado,
     xlab="k (grau total)", ylab= "C(k) (coef. cluster. por grau k)")
title(sub="Figura 7: Coef. Clusterização em função do grau k dos vértices ", 
      cex.sub = 0.75, font.sub = 3, col.sub = "black")

##################################

# connectivity

# clusters
count_components(Rede)


components(Rede)$csize

# histogram
plot(components(Rede)$csize)

reguladores.por.grupo <- c()
for(i in 1:count_components(Rede)){
  grupo <- names(components(Rede)$membership)[components(Rede)$membership == i]
  reguladores.por.grupo<- 
    c(reguladores.por.grupo, 
      sum(grupo %in% auxTF$gene_CCBH4851))
}

###

propor.reguladores.por.grupo <- reguladores.por.grupo/components(Rede)$csize

###
probab.n.reguladores.por.grupo <-
  table(reguladores.por.grupo)/count_components(Rede)

###
plot(probab.n.reguladores.por.grupo,
     log= "xy",
     ylim=c(0.01,1))


#####################################################

# Estudo de atividade dominante: ativacao ou desativação?

tab.modoRegula <- 
  table(arestas$`mode of regulation`)
tab.modoRegula <- tab.modoRegula[order(tab.modoRegula)]

# Dominant activity study: activation or non-activation
#
auto_regul <- 
  filter(arestas, `Regulator (TF)` == `Target gene`)
#
tab.modo.Auto_regul <- 
  table(auto_regul$`mode of regulation`)
tab.modo.Auto_regul <- tab.modo.Auto_regul[order(tab.modo.Auto_regul)]

#####################
# Path Lengths

### The diameter of a graph is the length of the longest geodesic.  
diameter(Rede,directed = TRUE,unconnected = TRUE)
diameter(Rede,directed = FALSE,unconnected = TRUE)


### the average path size:
mean_distance(Rede, directed = TRUE, unconnected = TRUE)
mean_distance(Rede, directed = FALSE, unconnected = TRUE)


####################
# Motifs de 3 vertices


triad_census(Rede)
triad_census(Rede)[9]
triad_census(Rede)[13]


Rede2 <- graph_from_data_frame(d = filter(arestas, 
                                          `mode of regulation` == "+" |
                                            `mode of regulation` == "-"))
los  "+" e "-" somente
triad_census(Rede2)
triad_census(Rede2)[9]
triad_census(Rede2)[13]

Rede3 <- graph_from_data_frame(d = filter(arestas, 
                                          `mode of regulation` == "+" )
                               # directed = TRUE, vertices = vert
)

triad_census(Rede3)
triad_census(Rede3)[9]
triad_census(Rede3)[13]


Rede4 <- graph_from_data_frame(d = filter(arestas, 
                                          `mode of regulation` == "-" ))

triad_census(Rede4)
triad_census(Rede4)[9]
triad_census(Rede4)[13]

####################
# Motifs of 2 vertices 

dyad_census(Rede)
dyad_census(Rede2)
dyad_census(Rede3)

#####
V(Rede)$name <- V(Rede)$rotulo 

# Kleinberg's hub centrality scores
# ref: https://www.sci.unich.it/~francesc/teaching/network/kleinberg.html

hubs.em.ordem.dec <-(hub_score(Rede)$vector)[order((hub_score(Rede)$vector), 
                                                   decreasing = TRUE)]

authority.em.ordem.dec <-
  authority_score(Rede)$vector[order((authority_score(Rede)$vector), 
                                     decreasing = TRUE)]


# top 10 (Kleinberg's hub centrality scores) 
hubs.em.ordem.dec[1:30]

# top 10 (Kleinberg's authority centrality scores) 
authority.em.ordem.dec[1:10]

## usando  k-out
k.hubs.em.ordem.dec <- degree(Rede, mode="out")[order(degree(Rede, mode="out"),
                                                      decreasing = TRUE)]

# top 10 k-out
k.hubs.em.ordem.dec[1:30]











