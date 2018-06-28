#Profile basic test


source("setup.R")
local({
  if(!require(profvis,lib.loc = user_path)) {
    install.packages("profvis",lib = user_path)
    if(!require(profvis,lib.loc = user_path)) {
      stop("Failed to install profiler")
    }
  }
  library(htmlwidgets,lib.loc = user_path)
  library(jsonlite,lib.loc = user_path)
  library(yaml,lib.loc = user_path)
  invisible(T)
})
print(selac_release)
setup_selac_for_profiling()


# basic loader to build further tests
load_inputs <- function(){
  tree <- read.tree("rokasYeast.tre")
  phy <- drop.tip(tree, "Calb")
  
  yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
  yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
  
  chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
  codon.data <- chars[phy$tip.label,]
  
  aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
  aa.optim <- apply(aa.data[, -1], 2, selac:::GetMaxName) #starting values for all, final values for majrule
  aa.optim.full.list <- aa.optim
  codon.freq.by.aa <- selac:::GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
  codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
  aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
  colnames(aa.optim.frame.to.add) <- colnames(codon.data)
  codon.data <- rbind(codon.data, aa.optim.frame.to.add)
  codon.data <- selac:::SitePattern(codon.data, includes.optimal.aa=TRUE)
  aa.optim = codon.data$optimal.aa
  
  codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
  
}

lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams
test_selac.gamma.quadrature <- function(){
  lSAC.c4mc.full(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5), 5)), 
                 codon.data=codon.data, phy=phy, aa.optim_array=aa.optim, 
                 codon.freq.by.aa=codon.freq.by.aa,  codon.freq.by.gene=codon.freq.by.gene, 
                 numcode=1, diploid=TRUE, aa.properties=NULL, 
                 volume.fixed.value=0.0003990333, 
                 nuc.model="GTR", 
                 codon.index.matrix=codon.index.matrix, 
                 include.gamma=TRUE, 
                 gamma.type="quadrature", 
                 ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, 
                 n.cores.by.gene.by.site=1)
  
}


test_selac.gamma.median <- function(){
  lSAC.c4mc.full(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5), 5)), 
                 codon_data=codon.data, phy=phy, aa.optim_array=aa.optim, 
                 codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, 
                 numcode=1, diploid=TRUE, aa.properties=NULL, 
                 volume.fixed.value=0.0003990333, 
                 nuc.model="GTR", 
                 codon.index.matrix=codon.index.matrix, 
                 include.gamma=TRUE, 
                 gamma.type="median", 
                 ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE)
  
}
test_selac.unrest <- function(){
  lSAC.c4mc.full(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, rep(1,11))), 
                 codon.data, phy, aa.optim_array=aa.optim, 
                 codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, 
                 numcode=1, diploid=TRUE, aa.properties=NULL, 
                 volume.fixed.value=0.0003990333, 
                 nuc.model="UNREST", 
                 codon.index.matrix, 
                 include.gamma=FALSE, 
                 ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, 
                 n.cores.by.gene.by.site=1)
  
}
test_selac.gtr <- function(){
  lSAC.c4mc.full(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5))), 
                 codon.data, phy, aa.optim_array=aa.optim, 
                 codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, 
                 numcode=1, diploid=TRUE, aa.properties=NULL, 
                 volume.fixed.value=0.0003990333, 
                 nuc.model="GTR", 
                 codon.index.matrix, 
                 include.gamma=FALSE, 
                 ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, 
                 n.cores.by.gene.by.site=1)
  
}

std.params = c(C.Phi.q.Ne = 4*4e-7*.5*5e6,
               alpha=1.829272, 
               beta=0.101799)
std.gamma=0.0003990333
std.base.freq = c(A=0.25,C=0.25,G=0.25)
std.poly.params = c(NA,NA)
std.gamma.shape = 5
## Notes on nuc.mutation.params:
# used as rates value in CreateNucleotideMutationMatrix(rates, model, base.freqs)->res
# either: length(base.freqs) == 4 && sum(base.freqs) == 1
# or: is.null(base.freqs) == TRUE
# dim(res) == c(4,4)
# rowSums(res) == rep(1,4)
## CreateNucleotideMutationMatrix with JC model
# rates = rates[1] (ie just uses first value)
## CreateNucleotideMutationMatrix with GTR model
# rates = rates[1:5] (ie just uses first 5 values)
## CreateNucleotideMutationMatrix with HKY model
# rates = rates[1:2] (ie just uses first 2 values)
## CreateNucleotideMutationMatrix with UNREST model
# rates = rates[1:11] (ie just uses first 11 values)
# 
# std.nuc.mutation.paramsA = c(1,1,1,1,1)
# std.nuc.mutation.paramsB = rep(1,11)
# std.nuc.mutation.paramsC = c(1,1,1,1,1)
std.nuc.params = mapply(rep_len,length.out=c(JC=1,GTR=5,HKY=2,UNREST=11),x=rep(1,4),
                        USE.NAMES =  T,SIMPLIFY = F)
test_selac_std <- function(phy, codon.data, 
                           nuc.model=c("JC", "GTR", "HKY", "UNREST"),
                           gamma.type=c("none", "median","quadrature","lognormal" ),
                           nCores=1){
  nuc.model=match.arg(nuc.model)
  gamma.type=mathc.arg(gamma.type)
  if(nuc.model == "HKY") stop("HKY model not implemented for GetLikelihoodSAC_CodonForManyCharGivenAllParams.")
  include.gamma = (gamma.type != "none")
  
  aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
  aa.optim <- apply(aa.data[, -1], 2, selac:::GetMaxName) #starting values for all, final values for majrule
  aa.optim.full.list <- aa.optim
  codon.freq.by.aa <- selac:::GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
  codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
  aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
  colnames(aa.optim.frame.to.add) <- colnames(codon.data)
  codon.data <- rbind(codon.data, aa.optim.frame.to.add)
  codon.data <- selac:::SitePattern(codon.data, includes.optimal.aa=TRUE)
  aa.optim = codon.data$optimal.aa
  
  codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
  model.params = std.params
  if(nuc.model != "UNREST")
    model.params=c(model.params,std.base.freq)
  model.params=c(model.params,std.nuc.params[[nuc.model]])
  if(include.gamma){
    model.params=c(model.params,std.gamma.shape)
    lSAC.c4mc.full(log(model.params), 
                   codon.data=codon.data, phy=phy, aa.optim_array=aa.optim, 
                   codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, 
                   numcode=1, diploid=TRUE, aa.properties=NULL, 
                   volume.fixed.value=std.gamma, 
                   nuc.model=nuc.model, 
                   codon.index.matrix=codon.index.matrix, 
                   include.gamma=TRUE, gamma.type=gamma.type,
                   ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, 
                   n.cores.by.gene.by.site=nCores)->res
    
  }else{
    lSAC.c4mc.full(log(model.params), 
                   codon.data=codon.data, phy=phy, aa.optim_array=aa.optim, 
                   codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, 
                   numcode=1, diploid=TRUE, aa.properties=NULL, 
                   volume.fixed.value=std.gamma, 
                   nuc.model=nuc.model, 
                   codon.index.matrix=codon.index.matrix, 
                   include.gamma=FALSE, 
                   ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, 
                   n.cores.by.gene.by.site=nCores)->res
  }
  return(res)
}
#round(selac.gtr, 3)

#get_test_key <- function(phy.source, nuc.model, gamma.type, nCores, seed)

