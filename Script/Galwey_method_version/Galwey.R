# test du seuil de signification par la methode de Galwey

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
library(snpStats)

# le fichier de données
load("X_NA.Rdata")

# 
# geno<-as.matrix(X_NA)
# rownames(geno)<-rownames(X_NA)
# 
# dim(geno)
# 
# geno<-apply(geno, 2, as.numeric)
# dim(geno)

geno<-as.matrix((X_NA))
geno<-apply(geno, 2, as.numeric)
geno<-round(geno,digits=0)
rownames(geno)<-rownames(X_NA)
#on cree la SNP matrix
test <- as(geno,"SnpMatrix")

# geno étant la matrice avec les individus en ligne et les SNPs en colonne.
# On la passe juste au format "SnpMatrix"

LDmat <- as.matrix(snpStats::ld(test, depth = ncol(test) - 1, stats = "R",symmetric	
=TRUE))

LDmat[1:12,1:12]

# on calcule la corrélation entre les SNPs deux à deux. 
# L'output est une matrice carrée mais seule la partie supérieure à la diagonale est calculée (gain de temps)

LDmat[lower.tri(LDmat)] = t(LDmat)[lower.tri(LDmat)] # on complète la partie inférieure

diag(LDmat) <- 1 # on complète la diagonale

Me<-meff(LDmat, method = "galwey") # on calcule le nombre effectif de tests indépendants avec la méthode de

capture.output(Me, file = "Me.txt")
sink(file = "genotype.txt")