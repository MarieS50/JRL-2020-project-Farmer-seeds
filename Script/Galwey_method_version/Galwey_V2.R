# test du seuil de signification par la methode de Galwey

# 1. Relative directory and libraries
# The directory will be set at the source of the script
# 1. Verify that it is ok on your computer
rm(list=ls())

path<-dirname(rstudioapi::getSourceEditorContext()$path)
last_slash<-unlist(lapply(strsplit(path, ''), function(x) which(x == "/")))
last_slash<-tail(last_slash, n=1)

# setwd(substr(path, 1 , last_slash))
setwd(path)
getwd()

library(anyLib)
anyLib(c("data.table", "apercu", "mlmm", "corpcor","plyr"))


# fonction meff
meff<- function(R, eigen = FALSE, method, ...) {
  
  # match method argument
  method <- match.arg(method, c("nyholt", "liji", "gao", "galwey"))
  
  if (eigen) {
    # can pass eigenvalues directly to function if eigen is TRUE
    
    if (!class(R) %in% c("numeric", "integer"))
      stop("eigenvalues are not numeric or integer.")
    
    evs <- R
    abs_evs <- abs(evs)
  } else {
    # check that R is symmetric
    if (!isSymmetric(R))
      stop("R is not symmetric.")
    
    # ensure that the correlation matrix is positive semi-definite
    #R <- as.matrix(nearPD(R)$mat)
    
    # get eigenvalues and absolute eigenvalues of R matrix
    evs <- eigen(R)$values
    abs_evs <- abs(evs)
  }
  
  k <- length(evs)
  
  if (method == "nyholt") {
    # effective number of tests (based on Nyholt, 2004)
    m <- 1 + (k - 1) * (1 - var(evs) / k)
  }
  
  if (method == "liji") {
    # effective number of tests (based on Li & Ji, 2005)
    # adding a small value to the eigenvalues to overcome numerical imprecisions
    abs_evs <- abs_evs + sqrt(.Machine$double.eps)
    m <- sum(ifelse(abs_evs >= 1, 1, 0) + (abs_evs - floor(abs_evs)))
  }
  
  if (method == "gao") {
    # effective number of tests (based on Gao, 2008)
    m <- which(cumsum(sort(abs_evs, decreasing = TRUE)) / sum(abs_evs) > 0.995)[1]
  }
  
  if (method == "galwey") {
    # effective number of tests (based on Galwey, 2009)
    evs[evs < 0] <- 0
    m <- sum(sqrt(evs))^2 / sum(evs)
  }
  
  # always round down estimated value
  m <- floor(m)
  
  return(m)
  
}


# installer snpStats
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("snpStats")

library(snpStats)

# le fichier de données
load("X_NA.Rdata")

G<-X_NA
rm(X_NA)

# SNP Physical position on  Zavitan
load("BREEDWHEAT_on_durum_physic_WEW2.Rdata")
head(names(BLAST))


# only snp with a blast
liste<-which(colnames(G) %in% BLAST[,1])
length(liste)

G<-G[,liste]
dim(G)

G<-as.matrix(G)

map<-BLAST[which(BLAST[,1] %in% colnames(G)),c(1,2,6)]
names(map)<-c("SNP","Chr","Pos")
head(map)

class(map)
dim(map)

map <- map[map$SNP %in% colnames(G), ]
dim(map)

map[,2] <-mapvalues(map[,2], from=c( "chrom", 
                                     "chr1A", "chr1B", "chr2A", "chr2B",
                                     "chr3A", "chr3B", "chr4A", "chr4B",
                                     "chr5A", "chr5B", "chr6A", "chr6B",
                                     "chr7A", "chr7B")
                    , to=c("chrom", 1:14))

map[,2]<-as.numeric(map[,2])

dim(map)
head(map)
tail(map)

# Chromosome & position ordering
map <- map[order(map$Pos), ]
map <- map[order(map$Chr), ]
head(map)
tail(map)


# example on chromosome 1 (you can develop a function from it instead of this horrific loop

Galwey<-0

for ( i in 1:14)
{ 
liste<-map[which(map$Chr==i),1]
geno<-G[,liste]
dim(geno)

geno<-apply(geno,2, round)

geno<-as.matrix(geno)
geno<-apply(geno, 2, as.numeric)
rownames(geno)<-rownames(G)
dim(geno)

# geno étant la matrice avec les individus en ligne et les SNPs en colonne.
# On la passe juste au format "SnpMatrix"
test <- as(geno,"SnpMatrix") 
test[1:2,1:2]

#Linkage disequiloibrium, be careful to modify test[,1:100] by test for the whole set of data
LDmat <- as.matrix(snpStats::ld(test, depth = ncol(test) - 1, stats = "R"))

# LDMat is a square matrix  with only the UPPER part filled 
# need to symetrize it

LDmat[lower.tri(LDmat)] = t(LDmat)[lower.tri(LDmat)] # on complète la partie inférieure
diag(LDmat) <- 1 # on complète la diagonale

# Application of the meff function
Galwey[i]<-meff(LDmat, method = "galwey") # on calcule le nombre effectif de tests indépendants avec la méthode de
print (paste("chr :",i, ":", Galwey[i]))

}

# on the whole set of chromosomes
Me<-sum(Galwey)
capture.output(Me, file="Me.txt")
