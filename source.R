################################################################# 
#####            Project1: Bayesian Networks                #####
#################################################################

# 0) Setup environment
# 
### Install the pakages "deal" and "yeastCC"
# install.packages("deal",repos="http://lib.stat.cmu.edu/R/CRAN",dependencies=TRUE)
# BiocManager::install("yeastCC")

### Load the packages
library("deal")
library("yeastCC")


# 1) Data Preprocessing
#
### Load the data
expr <- as.data.frame(t(exprs(yeastCC)[orf800,]))
cat("Observations:", nrow(expr), "\n")
cat("Genes:", ncol(expr), "\n")

### Replace missing data with gene expression median
# for-loop
for (j in 1:ncol(expr)) {
  index_na <- which(is.na(expr[,j]))  # indices of NAs
  expr[index_na, j] <- median(expr[,j], na.rm=T)  # replace by median
}

# apply() alternative 
medianfill <- function(exprCol) {
  exprCol[which(is.na(exprCol))] <- median(exprCol, na.rm=T)
  return(exprCol)
}
expr <- as.data.frame(apply(expr, 2, medianfill))

### Filter genes based on inter quartile range iqr
# for-loop
iqr <- vector(length=ncol(expr))
for (j in 1:ncol(expr)) {
  iqr[j] <- quantile(expr[,j], c(0.75)) - quantile(expr[,j], c(0.25))
}

# apply() alternative
gene.iqr <- function(exprCol) {
  quantile(exprCol, c(0.75)) - quantile(exprCol, c(0.25))
}
iqr <- apply(expr, 2, gene.iqr)

# Keep only genes with iqr > 1.6
expr = expr[, iqr > 1.6]   # keep only genes with large variation
names(expr)  # print selected genes
#YBR054W	YRO2
#YBR088C	POL30
#YER124C	DSE1
#YGL028C	SCW11
#YLR286C	CTS1
#YHR143W	DSE2
#YNL327W	EGT2
#YGR108W	CLB1
#YNR067C	DSE4
#YOL007C	CSI2

# deal package seems not to work if all nodes are continuous.
# Thus, we need to add a discrete "dummy" node to work around this issue. 
# Adding this node will not have any influence on the final results.
expr$dummy <- factor(rep("1", nrow(expr)))
head(expr)

#
# 2-6) Build Bayesian Network
#
### 2) Create prior structure 
# EITHER: Using default 
#G0 <- network(expr) # prior structure with no edges

# OR: Specify your own prior network manually:
# To insert an arrow from node 'A' to node 'B',
# first click node 'A' and then click node 'B'.
# Add an arrow from YOL007C to YBR088C, 
# and from YNL327W to YER124C, YHR143W and YNR067C.
# When the graph is finished, click 'stop',
# Then, inspect the local probability distribution of node i,
# by clicking on node i and then on the background.,
# Finish by clicking the center button 'Stop'.
G0  <- network(expr, specifygraph=TRUE, inspectprob=TRUE)

# We don't want any arrows starting from the "dummy" node, thus we construct a list of banned dependencies:
banlist(G0) <- matrix(c(11,11,11,11,11,11,11,11,11,11,1,2,3,4,5,6,7,8,9,10),ncol=2)
plot(G0)


### 3) Show local probability distribution
localprob(G0)
localprob(G0)$YBR088C

### 4) Compute joint prior distribution
prior <- jointprior(G0, 5)  # equivalent to imaginary sample size = 5

### 5) Learn the initial network

G0 <- getnetwork(learn(G0, expr, prior))
print(G0$score)

### 6) Search for optimal network (takes some time)
nwSearch <- autosearch(G0, expr, prior, removecycles=FALSE, trace=FALSE)
G <- getnetwork(nwSearch)
plot(G)

### Function for building an optimal network from expression data (from above)
build.optimal.network <- function(exprData, N=5) { 
  N0 <- network(exprData)
  print(head(exprData))
  banlist(N0) <- matrix(c(11,11,11,11,11,11,11,11,11,11,1,2,3,4,5,6,7,8,9,10),ncol=2)
  plot(N0)
  prior <- jointprior(N0, N)
  N0 <- getnetwork(learn(N0, exprData, prior))
  nwasarch <- autosearch(N0, exprData, prior, removecycles=FALSE, trace=FALSE)
  getnetwork(nwasarch)
}

### Custom function for plotting a BN 
plot.bn <- function(BN, file=NULL) {
  par(mar=c(0,0,0,0))
  plot(BN, cexscale=13, unitscale=27, arrowlength=0.1, xr=c(0, 350), yr=c(20,370))
  if (!is.null(file)) {
    plt <- recordPlot()
    pdf(file)
    replayPlot(plt)
    dev.off()
  }
}
BN <- build.optimal.network(expr)
plot(BN)
plot.bn(BN, file="~/Desktop/BNstar.pdf")


for (i in colnames(expr)){
  print(var(expr[i]))
}


# Task 7
var_table = lapply(expr, var)


#Task 8 
perturbate_data = function(data){
  var_table = lapply(data, var)
  for (i in colnames(data)){
    if (!is.numeric(data[,i])){
      next()
    }
    data[,i] = data[,i] + rnorm(nrow(data),0,as.numeric(var_table[i])/10)
  }
  return(data)
}

expr_perturbate_list = list()
for (i in seq(1,30)){
  expr_perturbate_list[[i]] = perturbate_data(expr)
}


#Task 9
select_experiment = function(data,position){
  return(data$YHR143W[position])
}

expr_data_pert_list = list()
for (i in seq(1,77)){
  expr_data_pert_list[[i]] = sapply(expr_perturbate_list, select_experiment , i)
}

require(ggplot2)
require(reshape2)

exp_data_frame = as.data.frame(expr_data_pert_list)
colnames(exp_data_frame) = rownames(expr)

ggplot(data = melt(exp_data_frame), aes(x=variable, y=value)) + theme_bw(base_size = 16) + theme(axis.text.x=element_text(angle=90,hjust=1)) + geom_boxplot()

pbn_list = list()
for (i in seq(1,30)){
  pbn_list[[i]] = build.optimal.network(expr_perturbate_list[[i]])
}

plot(build.optimal.network(expr_perturbate_list[[5]]))
