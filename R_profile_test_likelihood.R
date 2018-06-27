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
profile_prefix=paste0("likelihood_test_",selac_release)
print(selac_release)

try({
prof_obj <- profvis({
  ## Test 1
  set.seed(4)
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
  selac.gtr <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5))), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1)
  comparison <- identical(round(selac.gtr, 3), -7066.477)
  print(comparison)
  
  
}, prof_output = paste0(profile_prefix,"_GTR_rokasYeast.Rprof"))
htmlwidgets::saveWidget(prof_obj, paste0(profile_prefix,"_GTR_rokasYeast.Rprofvis"))
})
try({
prof_obj <- profvis({
  ## Test 2
  set.seed(4)
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
  selac.unrest <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, rep(1,11))), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="UNREST", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1)
  comparison <- identical(round(selac.unrest, 3), -7066.477)
  print(comparison)
}, prof_output = paste0(profile_prefix,"_UNREST_rokasYeast.Rprof"))
htmlwidgets::saveWidget(prof_obj, paste0(profile_prefix,"_UNREST_rokasYeast.Rprofvis"))
})
try({
  prof_obj <- profvis({
  
  ## Test 3
  set.seed(4)
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
  codon.index.matrix <- selac:::CreateCodonMutationMatrixIndex()
  selac_gamma <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5), 5)), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=TRUE, gamma.type="median", ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE)
  comparison <- identical(round(selac_gamma, 3), -6999.538)
  print(comparison)
}, prof_output =paste0(profile_prefix,"_UNRESTmedian_rokasYeast.Rprof"))

  htmlwidgets::saveWidget(prof_obj, paste0(profile_prefix,"_UNRESTmedian_rokasYeast.Rprofvis"))
})
try({
  prof_obj <- profvis({
  ## Test 4
  set.seed(4)
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
  selac_gamma <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5), 5)), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=TRUE, gamma.type="quadrature", ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1)
  comparison <- identical(round(selac_gamma, 3), -6998.618)
  print(comparison)
}, prof_output = paste0(profile_prefix,"_UNRESTquadrature_rokasYeast.Rprof"))
  htmlwidgets::saveWidget(prof_obj, paste0(profile_prefix,"_UNRESTquadrature_rokasYeast.Rprofvis"))
})