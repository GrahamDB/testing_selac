GetLikelihoodSAC_CodonForSingleCharGivenOptimum <- function(charnum=1, codon.data, phy, Q_codon, root.p=NULL, scale.factor, anc.indices, return.all=FALSE) {
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  
  nl <- nrow(Q_codon[[1]])
  #Now we need to build the matrix of likelihoods to pass to dev.raydisc:
  liks <- matrix(0, nb.tip + nb.node, nl)
  if(all(codon.data[,charnum+1] < 65)){
    #no need to subset
    liks[cbind(1:nb.tip,codon.data[,charnum+1])] <- 1
  } else {
    key<-codon.data[,charnum+1] < 65
    liks[cbind(which(key),codon.data[which(key),charnum+1])] <- 1
    liks[which(!key),] <- 1
    if(nl > 4){
      liks[which(!key),c(49, 51, 57)] <- 0
    }
  }
  
  #The result here is just the likelihood:
  result <- -FinishLikelihoodCalculation(phy=phy, liks=liks, Q=Q_codon, root.p=root.p, anc=anc.indices)
  if(return.all) stop("return all not currently implemented");
  return(result)
}



GetLikelihoodSAC_CodonForManyCharVaryingBySite <- function(codon.data, phy, Q_codon_array, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, aa.optim_array, codon_mutation_matrix, Ne, rates, numcode, diploid, n.cores.by.gene.by.site=1) {
  
  nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  final.likelihood.vector <- rep(NA, nsites.unique)
  #unique.aa <- GetMatrixAANames(numcode)
  
  #We rescale the codon matrix only:
  diag(codon_mutation_matrix) = 0
  diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
  scale.factor <- -sum(diag(codon_mutation_matrix) * codon.freq.by.gene, na.rm=TRUE)
  codon_mutation_matrix_scaled = codon_mutation_matrix * (1/scale.factor)
  #Finish the Q_array codon mutation matrix multiplication here:
  for(k in 1:21){
    if(diploid == TRUE){
      Q_codon_array[,,.unique.aa[k]] = (2 * Ne) * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
    }else{
      Q_codon_array[,,.unique.aa[k]] = Ne * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
    }
    diag(Q_codon_array[,,.unique.aa[k]]) = 0
    diag(Q_codon_array[,,.unique.aa[k]]) = -rowSums(Q_codon_array[,,.unique.aa[k]])
  }
  
  #Put the na.rm=TRUE bit here just in case -- when the amino acid is a stop codon, there is a bunch of NaNs. Should be fixed now.
  #scale.factor <- -sum(Q_codon_array[DiagArray(dim(Q_codon_array))] * equilibrium.codon.freq, na.rm=TRUE)
  phy <- reorder(phy, "pruningwise")
  
  ## This is obviously not very elegant, but not sure how else to code it to store this stuff in this way -- WORK IN PROGRESS:
  tempGetAAExpQt <- local({
    p0=phy
    Qca=Q_codon_array
    r0=rates
    function(aa) {
      GetExpQt(phy=p0, Q=Qca[,,aa], scale.factor=NULL, rates=r0)
    } })
  expQt <- NULL
  expQt <- mclapply(c("K", "N", "T", "R", "S", 
                      "I", "M", "Q", "H", "P", 
                      "L", "E",  "D", "A", "G", 
                      "V", "Y", "C", "W", "F"),
                    FUN=tempGetAAExpQt, 
                    mc.cores=n.cores.by.gene.by.site)
  names(expQt) <- c("K", "N", "T", "R", "S", 
                    "I", "M", "Q", "H", "P", 
                    "L", "E",  "D", "A", "G", 
                    "V", "Y", "C", "W", "F")
  #Generate matrix of root frequencies for each optimal AA:
  root.p_array <- matrix(codon.freq.by.aa, nrow=dim(Q_codon_array)[2], ncol=21)
  root.p_array <- t(root.p_array)
  root.p_array <- root.p_array / rowSums(root.p_array)
  rownames(root.p_array) <- .unique.aa
  
  phy.sort <- reorder(phy, "pruningwise")
  anc.indices <- unique(phy.sort$edge[,1])
  MultiCoreLikelihoodBySite <- function(i){
    tmp <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=i, 
                                                           codon.data=codon.data$unique.site.patterns, 
                                                           phy=phy.sort, 
                                                           Q_codon=expQt[[aa.optim_array[i]]], 
                                                           root.p=root.p_array[aa.optim_array[i],], 
                                                           scale.factor=scale.factor, 
                                                           anc.indices=anc.indices, 
                                                           return.all=FALSE)
    return(tmp)
  }
  final.likelihood.vector.mc <- unlist(mclapply(1:nsites.unique, MultiCoreLikelihoodBySite, mc.cores=n.cores.by.gene.by.site))

  return(final.likelihood.vector.mc)
}


GetLikelihoodNucleotideForManyCharVaryingBySite <- function(nuc.data, phy, nuc.mutation.rates, include.gamma=FALSE, rates.k=NULL, ncats=NULL, root.p_array=NULL, n.cores.by.gene.by.site=1) {
  nsites.unique <- dim(nuc.data$unique.site.patterns)[2]-1
  final.likelihood.vector <- rep(NA, nsites.unique)
  if(is.null(root.p_array)) {
    #Generate matrix of equal frequencies for each site:
    root.p_array <- rep(0.25, 4)
  }
  #Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
  diag(nuc.mutation.rates) = 0
  nuc.mutation.rates = t(nuc.mutation.rates * root.p_array)
  diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
  scale.factor <- -sum(diag(nuc.mutation.rates) * root.p_array)
  
  expQt <- GetExpQt(phy=phy, Q=nuc.mutation.rates, scale.factor=scale.factor, rates=rates.k)
  phy.sort <- reorder(phy, "pruningwise")
  anc.indices <- unique(phy.sort$edge[,1])
  
  MultiCoreLikelihoodBySite <- function(nsite.index){
    tmp <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=nsite.index, codon.data=nuc.data$unique.site.patterns, phy=phy, Q_codon=expQt, root.p=root.p_array, scale.factor=scale.factor, anc.indices=anc.indices, return.all=FALSE)
    return(tmp)
  }
  final.likelihood.vector <- unlist(mclapply(1:nsites.unique, MultiCoreLikelihoodBySite, mc.cores=n.cores.by.gene.by.site))
  return(final.likelihood.vector)
}


GetExpQt <- function(phy, Q, scale.factor, rates=NULL){
  
  if(!is.null(scale.factor)){
    Q.scaled = Q * (1/scale.factor)
  }else{
    Q.scaled = Q
  }
  if(!is.null(rates)){
    Q.scaled = Q.scaled * rates
  }
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  expQt <- as.list(numeric(nb.tip + nb.node))
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  #phy <- reorder(phy, "pruningwise")
  #Obtain an object of all the unique ancestors
  anc <- unique(phy$edge[,1])
  for (i  in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows <- which(phy$edge[,1]==focal)
    desNodes <- phy$edge[desRows,2]
    for (desIndex in sequence(length(desRows))){
      expQt[[desNodes[desIndex]]] <- internal_expm(Q.scaled * phy$edge.length[desRows[desIndex]])
    }
  }
  return(expQt)
}

internal_expm <- function (x, order = 8, 
                           trySym = TRUE, tol = .Machine$double.eps)
{
  stopifnot(is.numeric(x) || (isM <- inherits(x, "dMatrix")) || 
              inherits(x, "mpfrMatrix"))
  if (length(d <- dim(x)) != 2) 
    stop("argument is not a matrix")
  if (d[1] != d[2]) 
    stop("matrix not square")
  method <- "Higham08.b"
  preconditioning = "2bal"
  checkSparse <- !nzchar(Sys.getenv("R_EXPM_NO_DENSE_COERCION"))
  isM <- !is.numeric(x) && isM
  if (isM && checkSparse) {
    if (is(x, "sparseMatrix")) {
      x <- as(x, "denseMatrix")
    }
  }
  stopifnot(is.matrix(x))
  res <- expm.Higham08(x, balancing = TRUE)
  return(res)
}



#Step 2: Finish likelihood by taking our already exponentiated Q down the tree and simply re-traverse the tree and multiply by the observed likelihood.
FinishLikelihoodCalculation <- function(phy, liks, Q, root.p, anc){
  
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  if(any(root.p < 0) | any(is.na(root.p))){
    return(1000000)
  }
  #Obtain an object of all the unique ancestors
  for (i  in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows <- which(phy$edge[,1]==focal)
    desNodes <- phy$edge[desRows,2]
    v <- 1
    for (desIndex in desNodes){
      if(desIndex <= nb.tip){
        if(sum(liks[desIndex,]) < 2){
          v <- v * (Q[[desIndex]] %*% liks[desIndex,])
        }
      }else{
        v <- v * (Q[[desIndex]] %*% liks[desIndex,])
      }
    }
    comp[focal] <- sum(v)
    liks[focal,] <- v/comp[focal]
  }
  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
    return(1000000)
  }
  else{
    loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
    if(is.infinite(loglik)){return(1000000)}
  }
  loglik
}

