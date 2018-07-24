#selac_wrap_methods

OptimizeEdgeLengths <- function(x, par.mat, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, n.cores.by.gene, n.cores.by.gene.by.site=1, estimate.importance=FALSE, neglnl=FALSE, HMM=FALSE) {
  
  if(logspace) {
    x <- exp(x)
  }
  phy$edge.length = x
  if(is.null(aa.optim_array)){ stop("skipping this implementation") 
  } else {
    
    if(nuc.model == "JC"){
      max.par = 6
    }
    if(nuc.model == "GTR"){
      max.par = 6 + 5
    }
    if(nuc.model == "UNREST"){
      max.par = 3 + 11
    }
    if(include.gamma == TRUE){
      max.par = max.par + 1
    }
    if(k.levels > 0){
      max.par = max.par + 2
    }
    MultiCoreLikelihood <- function(partition.index){
      codon.data = NULL
      codon.data$unique.site.patterns = codon.site.data[[partition.index]]
      codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
      likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
      return(likelihood.tmp)
    }
    #This orders the nsites per partition in decreasing order (to increase efficiency):
    partition.order <- 1:n.partitions
    likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
  }
  return(likelihood)
}

