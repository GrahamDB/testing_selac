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
# setup_selac_for_profiling()
test_selac_hmm <- function(phy, 
                           fasta.file, 
                           nuc.model=c("JC", "GTR", "HKY", "UNREST"),
                           gamma.type=c("none", "median","quadrature","lognormal" ),
                           nCores=1){
  nuc.model=match.arg(nuc.model)
  gamma.type=match.arg(gamma.type)
  if(nuc.model == "HKY") stop("HKY model not implemented for GetLikelihoodSAC_CodonForManyCharGivenAllParams.")
  include.gamma = (gamma.type != "none")
  if(!include.gamma) gamma.type = "quadrature"
  
  tmp.gene <- read.dna(fasta.file, format="fasta")
  tmp.gene <- as.list(as.matrix(cbind(tmp.gene)))
  
  chars <- selac:::DNAbinToCodonNumeric(tmp.gene)
  codon.data <- chars[phy$tip.label,c(1,1+sample(ncol(chars)-1,10))]
  codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
  codon.data <- selac:::SitePattern(codon.data)
  
  
  
  codon.index.matrix = selac:::CreateCodonMutationMatrixIndexEvolveAA()
  model.params = hmm.params
  if(nuc.model != "UNREST")
    model.params=c(model.params,std.base.freq)
  model.params=c(model.params,std.nuc.params[[nuc.model]])
  
  lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA
  if(include.gamma)
    model.params=c(model.params,std.gamma.shape)
  model.params = c(model.params,std.sel.reg)
  
  lSAC.c4mc.full(log(model.params), 
                 codon.data=codon.data, 
                 phy=phy,  
                 codon.freq.by.aa=NULL, 
                 codon.freq.by.gene=codon.freq.by.gene, 
                 numcode=1, 
                 diploid=TRUE, 
                 aa.properties=NULL, 
                 volume.fixed.value=std.gamma, 
                 nuc.model=nuc.model, 
                 codon.index.matrix=codon.index.matrix, 
                 include.gamma=include.gamma, 
                 gamma.type=gamma.type,
                 ncats=4, 
                 k.levels=0, 
                 logspace=TRUE, 
                 verbose=TRUE, 
                 n.cores.by.gene.by.site=nCores,
                 estimate.importance=FALSE) -> res
  
  return(res)
}

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
load_rokasYeast <- function(){
  tree <- read.tree("rokasYeast.tre")
  phy <- drop.tip(tree, "Calb")
  
  yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
  yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
  
  chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
  codon.data <- chars[phy$tip.label,]
  list(input.key="rokasYeast",phy=phy,codon.data=codon.data)
}


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
hmm.params = c(C.Phi.q.Ne = 2,
               alpha=1.829272, 
               beta=0.101799)
std.sel.reg = 0.01
## Notes on nuc.mutation.params:
# used as rates value in selac:::CreateNucleotideMutationMatrix(rates, model, base.freqs)->res
# either: length(base.freqs) == 4 && sum(base.freqs) == 1
# or: is.null(base.freqs) == TRUE
# dim(res) == c(4,4)
# rowSums(res) == rep(1,4)
## selac:::CreateNucleotideMutationMatrix with JC model
# rates = rates[1] (ie just uses first value)
## selac:::CreateNucleotideMutationMatrix with GTR model
# rates = rates[1:5] (ie just uses first 5 values)
## selac:::CreateNucleotideMutationMatrix with HKY model
# rates = rates[1:2] (ie just uses first 2 values)
## selac:::CreateNucleotideMutationMatrix with UNREST model
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
  gamma.type=match.arg(gamma.type)
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
  lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams
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



run_profile <- function(src_data,nuc.model,gamma.model,seed,nCores){
  set.seed(seed)
  cat(sprintf("Start: %s_%s_%s_%s_%i_%i\n",
              src_data$input.key,
              nuc.model,
              gamma.model,
              selac_release,
              nCores,
              seed))
  profile_prefix=sprintf("%s_%s_%s_%s_%i_%i",
                         src_data$input.key,
                         nuc.model,
                         gamma.model,
                         selac_release,
                         nCores,
                         seed)
  model.LL=NA
  try({
    prof_obj <- profvis({
      model.LL=test_selac_std(src_data$phy,
                              src_data$codon.data,
                              nuc.model = nuc.model,
                              gamma.type = gamma.model,
                              nCores = nCores)
    }, prof_output = paste0(profile_prefix,".Rprof"))
    htmlwidgets::saveWidget(prof_obj, 
                            file=paste0(profile_prefix,".Rprofvis.html"))
  })
  cat(sprintf("End: %s_%s_%s_%s_%i_%i\tLL: %0.3f\n",
              src_data$input.key,
              nuc.model,
              gamma.model,
              selac_release,
              nCores,
              seed,
              model.LL))
  if(!file.exists(paste0(src_data$input.key,"_LL_log.csv")))
    cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
        file=paste0(src_data$input.key,"_LL_log.csv"),
        append = T)
  cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
              src_data$input.key,
              nuc.model,
              gamma.model,
              selac_release,
              nCores,
              seed,
              model.LL),
      file=paste0(src_data$input.key,"_LL_log.csv"),
      append = T)
  model.LL
}

run_simple_selac_optimize <- function(seed=sample.int(1e6,1),ref="v1.6.1-rc1"){
  setup_selac_for_profiling(ref=ref)
  profile_prefix=sprintf("%s_%s_%s_%s_%i_%i",
                         "selac19XX",
                         "GTR",
                         "noneXquadrature",
                         selac_release,
                         3,
                         seed)
  src.key="selac19XX"
  set.seed(seed)
  cat(sprintf("Start: %s\n",
              profile_prefix))
  tree<-read.tree('selac_paper_data/SalichosRokas.tre')
  
  result=list(loglik=NA)
  nuc.model = 'GTR'
  gamma.type="noneXquadrature"
  nCores=3
  try({
    prof_obj <- profvis({
  ## start.from.mle set to allow manual specification of fasta files
  # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
  # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
  result <- SelacOptimize(codon.data.path = 'tmp_data/', phy = tree, n.partitions=3,
                          edge.length = 'optimize', optimal.aa = 'none', data.type='codon',
                          codon.model = 'GY94', nuc.model = 'GTR', 
                          include.gamma = FALSE, gamma.type='quadrature', ncats = 4, numcode = 1,
                          diploid = FALSE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                          n.cores.by.gene  = 3, max.restarts = 1, max.evals=20)
    }, prof_output = paste0(profile_prefix,".Rprof"))
    htmlwidgets::saveWidget(prof_obj, 
                            file=paste0(profile_prefix,".Rprofvis.html"))
  })
  cat(sprintf("End: %s\tLL: %0.3f\n",
              profile_prefix,
              result$loglik))
  if(!file.exists(paste0(src.key,"_LL_log.csv")))
    cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
        file=paste0(src.key,"_LL_log.csv"),
        append = T)
  cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
              src.key,
              nuc.model,
              gamma.type,
              selac_release,
              nCores,
              seed,
              result$loglik),
      file=paste0(src.key,"_LL_log.csv"),
      append = T)
  save(result,file=sprintf('selac_paper_output/yeastSalRokSelacGTRG_quad_%s.Rdata',profile_prefix))
  result$loglik
}

run_full_selac_optimize <- function(seed=sample.int(1e6,1),ref="v1.6.1-rc1", nCores=3){
  setup_selac_for_profiling(ref=ref)
  profile_prefix=sprintf("%s_%s_%s_%s_%i_%i",
                         "selacFULLb",
                         "GTR",
                         "noneXquadrature",
                         selac_release,
                         nCores,
                         seed)
  src.key="selacFULLb"
  set.seed(seed)
  cat(sprintf("Start: %s\n",
              profile_prefix))
  tree<-read.tree('selac_paper_data/SalichosRokas.tre')
  
  result=list(loglik=NA)
  nuc.model = 'GTR'
  gamma.type="noneXquadrature"
  #nCores=3
  try({
    prof_obj <- profvis({
      ## start.from.mle set to allow manual specification of fasta files
      # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
      # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
      result <- SelacOptimize(codon.data.path = 'selac_paper_data/', phy = tree, 
                              edge.length = 'optimize', optimal.aa = 'none', data.type='codon',
                              codon.model = 'GY94', nuc.model = 'GTR', 
                              include.gamma = FALSE, gamma.type='quadrature', ncats = 4, numcode = 1,
                              diploid = FALSE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                              n.cores.by.gene  = nCores, max.restarts = 1, max.evals=20)
    }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.5)
    save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
    htmlwidgets::saveWidget(prof_obj, 
                            file=paste0(profile_prefix,".Rprofvis.html"))
  })
  cat(sprintf("End: %s\tLL: %0.3f\n",
              profile_prefix,
              result$loglik))
  if(!file.exists(paste0(src.key,"_LL_log.csv")))
    cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
        file=paste0(src.key,"_LL_log.csv"),
        append = T)
  cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
              src.key,
              nuc.model,
              gamma.type,
              selac_release,
              nCores,
              seed,
              result$loglik),
      file=paste0(src.key,"_LL_log.csv"),
      append = T)
  save(result,file=sprintf('selac_paper_output/yeastSalRokSelacGTRG_quad_%s.Rdata',profile_prefix))
  result$loglik
}

run_test_ecoli_optimize <- function(seed=sample.int(1e6,1),ref="v1.6.1-rc1", nCores=3){
  setup_selac_for_profiling(ref=ref)
  src.key="ecoliTEST"
  nuc.model = 'UNREST'
  gamma.type="quadrature"
  profile_prefix=sprintf("%s_%s_%s_%s_%i_%i",
                         src.key,
                         nuc.model,
                         gamma.type,
                         selac_release,
                         nCores,
                         seed)
  set.seed(seed)
  cat(sprintf("Start: %s\n",
              profile_prefix))
  tree<-read.tree('kosi07_data/kosi07_codonphyml_tree_TEM.newick')
  fasta.file="kosi07_data/aligned_KOSI07_TEM.fasta"
  output.file.name=sprintf('ecoli_output/%s_restart.Rdata',profile_prefix)
  result=list(loglik=NA)
  opt.aa.type <- "optimize"
  # random starting values
  starting.vals <- matrix(runif(n = 15, min = 0.01, max = 5), ncol = 15, nrow = 1)
  tree$edge.length <- runif(nrow(tree$edge), 0.01, 3)
  
  #nCores=3
  try({
    prof_obj <- profvis({
      ## start.from.mle set to allow manual specification of fasta files
      # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
      # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
      result <- SelacOptimize(codon.data.path = 'kosi07_data/', phy = tree, 
                              edge.length = 'optimize', optimal.aa = opt.aa.type, data.type='codon',
                              codon.model = 'selac', nuc.model = nuc.model, edge.linked=TRUE,
                              include.gamma = TRUE, gamma.type='quadrature', ncats = 4, numcode = 1,
                              diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                              n.cores.by.gene  = nCores, n.cores.by.gene.by.site=1,
                              max.restarts = 1, max.evals=20, max.tol=1e-2, max.iterations = 15,
                              fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE,
                              output.restart.filename=output.file.name)
      # output.restart.filename=output.file.name, start.from.mle = TRUE,
      # mle.matrix=starting.vals, tol.step=1, partition.order = fasta.file)
    }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.5)
    save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
    htmlwidgets::saveWidget(prof_obj, 
                            file=paste0(profile_prefix,".Rprofvis.html"))
  })
  cat(sprintf("End: %s\tLL: %0.3f\n",
              profile_prefix,
              result$loglik))
  if(!file.exists(paste0(src.key,"_LL_log.csv")))
    cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
        file=paste0(src.key,"_LL_log.csv"),
        append = T)
  cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
              src.key,
              nuc.model,
              gamma.type,
              selac_release,
              nCores,
              seed,
              result$loglik),
      file=paste0(src.key,"_LL_log.csv"),
      append = T)
  cat("SELAC Done. saving results\n")
  
  result$seed <- seed
  result$startingValues <- starting.vals
  result$startingTree <- tree
  
  save(result,file=sprintf('ecoli_output/%s_result.Rdata',profile_prefix))
  result$loglik
}



run_ecoli_profile_mode <- function(mode=c("SHORTTEST","TEST","SHORT",
                                          "SHORTTESTHMM","SHORTHMM",
                                          "FASTHMMTEST","HMMEVAL50","HMMEVALFULL",
                                          "FASTHMMDEBUG","FASTHMMSINGLEDEBUG"),
                                   seed=sample.int(1e6,1),
                                   codon.model=c("selac","none","GY94","YN98"),
                                   nuc.model=c("GTR","UNREST","JC"),
                                   ref="v1.6.1-rc1", 
                                   include.gamma=T,
                                   gamma.type=c("quadrature","median","lognormal","none"),
                                   nCores=1){
  setup_selac_for_profiling(ref=ref)
  mode=match.arg(mode)
  src.key=paste0("ecoli",mode)
  codon.model = match.arg(codon.model)
  nuc.model = match.arg(nuc.model)
  if(!include.gamma)
    { gamma.type="quadrature"; gamma.mode="none";}
  else {
    gamma.mode=gamma.type=match.arg(gamma.type)
  } 
  if(gamma.type=="none"){
    include.gamma=F
    gamma.type="quadrature"
    gamma.mode="none"
  }
    
  profile_prefix=sprintf("%s_%s_%s_%s_%s_%i_%i",
                         src.key,
                         codon.model,
                         nuc.model,
                         gamma.mode,
                         selac_release,
                         nCores,
                         seed)
  if(file.exists(sprintf('ecoli_output/%s_result.Rdata',profile_prefix))){
    try({
      load(file=sprintf('ecoli_output/%s_result.Rdata',profile_prefix))
      if(!is.null(result$loglik) && is.finite(result$loglik)) {
        cat(sprintf("Skip: %s\n",
                    profile_prefix))
        return(result$loglik)
      } 
    })
    cat(sprintf("Rebuilding: %s\n",
                profile_prefix))
  }
  set.seed(seed)
  cat(sprintf("Start: %s\n",
              profile_prefix))
  tree<-read.tree('kosi07_data/kosi07_codonphyml_tree_TEM.newick')
  fasta.file="kosi07_data/aligned_KOSI07_TEM.fasta"
  output.file.name=sprintf('ecoli_output/%s_restart.Rdata',profile_prefix)
  result=list(loglik=NA)
  opt.aa.type <- "optimize"
  # random starting values
  starting.vals <- matrix(runif(n = 15, min = 0.01, max = 5), ncol = 15, nrow = 1)
  tree$edge.length <- runif(nrow(tree$edge), 0.01, 3)
  
  #nCores=3
  if(mode=="TEST"){
  try({
    prof_obj <- profvis({
      ## start.from.mle set to allow manual specification of fasta files
      # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
      # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
      result <- SelacOptimize(codon.data.path = 'kosi07_data/', phy = tree, 
                              edge.length = 'optimize', optimal.aa = opt.aa.type, data.type='codon',
                              codon.model = codon.model, nuc.model = nuc.model, edge.linked=TRUE,
                              include.gamma = include.gamma, gamma.type=gamma.type, ncats = 4, numcode = 1,
                              diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                              n.cores.by.gene  = 1, n.cores.by.gene.by.site=nCores,
                              max.restarts = 1, max.evals=20, max.tol=1e-2, max.iterations = 15,
                              fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE,
                              output.restart.filename=output.file.name)
      # output.restart.filename=output.file.name, start.from.mle = TRUE,
      # mle.matrix=starting.vals, tol.step=1, partition.order = fasta.file)
    }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.5)
    save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
    # htmlwidgets::saveWidget(prof_obj, 
    #                         file=paste0(profile_prefix,".Rprofvis.html"))
  })
  } else if(mode=="SHORTTEST"){
    try({
      prof_obj <- profvis({
        ## start.from.mle set to allow manual specification of fasta files
        # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
        # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
        result <- SelacOptimize(codon.data.path = 'kosi07_data/', phy = tree, 
                                edge.length = 'optimize', optimal.aa = opt.aa.type, data.type='codon',
                                codon.model = codon.model, nuc.model = nuc.model, edge.linked=TRUE,
                                include.gamma = include.gamma, gamma.type=gamma.type, ncats = 4, numcode = 1,
                                diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                                n.cores.by.gene  = 1, n.cores.by.gene.by.site=nCores,
                                max.restarts = 1, max.evals=1, max.tol=1e-2, max.iterations = 1,
                                fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE,
                                output.restart.filename=output.file.name)
        # output.restart.filename=output.file.name, start.from.mle = TRUE,
        # mle.matrix=starting.vals, tol.step=1, partition.order = fasta.file)
      }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.5)
      save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
      # htmlwidgets::saveWidget(prof_obj, 
      #                         file=paste0(profile_prefix,".Rprofvis.html"))
    })
  }else if(mode=="SHORTTESTHMM"){
    # HMM code requires starting edge length < 0.5 and > 1e-8
    tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
    try({
      prof_obj <- profvis({
        ## start.from.mle set to allow manual specification of fasta files
        # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
        # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
        result <- SelacHMMOptimize(codon.data.path = 'kosi07_data/', phy = tree, 
                                edge.length = 'optimize', data.type='codon',
                                codon.model = codon.model, nuc.model = nuc.model, edge.linked=TRUE,
                                include.gamma = include.gamma, gamma.type=gamma.type, ncats = 4, numcode = 1,
                                diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = TRUE,
                                n.cores.by.gene  = 1, n.cores.by.gene.by.site=nCores,
                                max.restarts = 1, max.evals=1, max.tol=1e-2, 
                                fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE,
                                output.restart.filename=output.file.name, max.iterations=1)
        # output.restart.filename=output.file.name, start.from.mle = TRUE,
        # mle.matrix=starting.vals, tol.step=1, partition.order = fasta.file)
      }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.5)
      save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
      # htmlwidgets::saveWidget(prof_obj, 
      #                         file=paste0(profile_prefix,".Rprofvis.html"))
    })
  }else if(mode=="SHORT"){
    try({
      prof_obj <- profvis({
        ## start.from.mle set to allow manual specification of fasta files
        # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
        # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
        result <- SelacOptimize(codon.data.path = 'kosi07_data/', phy = tree, 
                                edge.length = 'optimize', optimal.aa = opt.aa.type, data.type='codon',
                                codon.model = codon.model, nuc.model = nuc.model, edge.linked=TRUE,
                                include.gamma = include.gamma, gamma.type=gamma.type, ncats = 4, numcode = 1,
                                diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                                n.cores.by.gene  = 1, n.cores.by.gene.by.site=nCores,
                                max.restarts = 1, max.evals=1, max.tol=1e-2, max.iterations = 1,
                                fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE,
                                output.restart.filename=output.file.name)
        # output.restart.filename=output.file.name, start.from.mle = TRUE,
        # mle.matrix=starting.vals, tol.step=1, partition.order = fasta.file)
      }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.05)
      save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
      # htmlwidgets::saveWidget(prof_obj, 
      #                         file=paste0(profile_prefix,".Rprofvis.html"))
    })
  } else  if(mode=="SHORTHMM"){
    # HMM code requires starting edge length < 0.5 and > 1e-8
    tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
    try({
      prof_obj <- profvis({
        ## start.from.mle set to allow manual specification of fasta files
        # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
        # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
        result <- SelacHMMOptimize(codon.data.path = 'kosi07_data/', phy = tree, 
                                edge.length = 'optimize',data.type='codon',
                                codon.model = codon.model, nuc.model = nuc.model, edge.linked=TRUE,
                                include.gamma = include.gamma, gamma.type=gamma.type, ncats = 4, numcode = 1,
                                diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                                n.cores.by.gene  = 1, n.cores.by.gene.by.site=nCores,
                                max.restarts = 1, max.evals=1, max.tol=1e-2,
                                fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE,
                                output.restart.filename=output.file.name, max.iterations=1)
        # output.restart.filename=output.file.name, start.from.mle = TRUE,
        # mle.matrix=starting.vals, tol.step=1, partition.order = fasta.file)
      }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.5)
      save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
      # htmlwidgets::saveWidget(prof_obj, 
      #                         file=paste0(profile_prefix,".Rprofvis.html"))
    })
  } else if(mode=="FASTHMMDEBUG") {
    
    if(!file.exists("ecoli_chars.RData")){
      tmp.gene <- read.dna(fasta.file, format="fasta")
      tmp.gene <- as.list(as.matrix(cbind(tmp.gene)))
      
      chars <- selac:::DNAbinToCodonNumeric(tmp.gene)
      save(chars,file="ecoli_chars.RData")
    } else {
      load(file="ecoli_chars.RData")
    }
    codon.data <- chars[tree$tip.label,c(1,1+sample(ncol(chars)-1,10))]
    codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
    codon.data <- selac:::SitePattern(codon.data)
    codon.index.matrix = selac:::CreateCodonMutationMatrixIndexEvolveAA()
    model.params = hmm.params
    if(nuc.model != "UNREST")
      model.params=c(model.params,std.base.freq)
    model.params=c(model.params,std.nuc.params[[nuc.model]])
    
    lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA
    if(include.gamma)
      model.params=c(model.params,std.gamma.shape)
    model.params = c(model.params,std.sel.reg)
    tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
    # try({
        lSAC.c4mc.full(log(model.params), 
                       codon.data=codon.data, 
                       phy=tree,  
                       codon.freq.by.aa=NULL, 
                       codon.freq.by.gene=codon.freq.by.gene, 
                       numcode=1, 
                       diploid=TRUE, 
                       aa.properties=NULL, 
                       volume.fixed.value=std.gamma, 
                       nuc.model=nuc.model, 
                       codon.index.matrix=codon.index.matrix, 
                       include.gamma=include.gamma, 
                       gamma.type=gamma.type,
                       ncats=4, 
                       k.levels=0, 
                       logspace=TRUE, 
                       verbose=TRUE, 
                       n.cores.by.gene.by.site=nCores,
                       estimate.importance=FALSE) -> result$loglik
    # })
  }else  if(mode=="FASTHMMSINGLEDEBUG") {
    
    if(!file.exists("ecoli_chars.RData")){
      tmp.gene <- read.dna(fasta.file, format="fasta")
      tmp.gene <- as.list(as.matrix(cbind(tmp.gene)))
      
      chars <- selac:::DNAbinToCodonNumeric(tmp.gene)
      save(chars,file="ecoli_chars.RData")
    } else {
      load(file="ecoli_chars.RData")
    }
    codon.data <- chars[tree$tip.label,c(1,1+sample(ncol(chars)-1,10))]
    codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
    codon.data <- selac:::SitePattern(codon.data)
    codon.index.matrix = selac:::CreateCodonMutationMatrixIndexEvolveAA()
    model.params = hmm.params
    if(nuc.model != "UNREST")
      model.params=c(model.params,std.base.freq)
    model.params=c(model.params,std.nuc.params[[nuc.model]])
    
    # lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA
    GetLikelihoodSAC_CodonForManyCharVaryingBySiteEvolvingAA <- function(codon.data, phy, Q_codon_array, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, aa.optim_array, codon_mutation_matrix, Ne, rates, numcode, diploid, n.cores.by.gene.by.site=1, verbose=FALSE){
      
      nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
      final.likelihood.vector <- rep(NA, nsites.unique)
      
      #We rescale the codon matrix only:
      diag(codon_mutation_matrix) = 0
      diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
      scale.factor <- -sum(diag(codon_mutation_matrix) * codon.freq.by.gene, na.rm=TRUE)
      codon_mutation_matrix_scaled = codon_mutation_matrix * (1/scale.factor)
      #Finish the Q_array codon mutation matrix multiplication here:
      if(diploid == TRUE){
        Q_codon_array = (2 * Ne) * codon_mutation_matrix_scaled * Q_codon_array
      }else{
        Q_codon_array = Ne * codon_mutation_matrix_scaled * Q_codon_array
      }
      diag(Q_codon_array) = 0
      diag(Q_codon_array) = -rowSums(Q_codon_array)
      #Put the na.rm=TRUE bit here just in case -- when the amino acid is a stop codon, there is a bunch of NaNs. Should be fixed now.
      #scale.factor <- -sum(Q_codon_array[DiagArray(dim(Q_codon_array))] * equilibrium.codon.freq, na.rm=TRUE)
      
      ## This is obviously not very elegant, but not sure how else to code it to store this stuff in this way -- WORK IN PROGRESS:
      #expQt <- GetExpQt(phy=phy, Q=Q_codon_array, scale.factor=NULL, rates=rates)
      #Generate matrix of root frequencies for each optimal AA:
      root.p_array <- codon.freq.by.gene
      #root.p_array <- t(root.p_array)
      #root.p_array <- root.p_array / rowSums(root.p_array)
      #rownames(root.p_array) <- .unique.aa
      phy.sort <- reorder(phy, "pruningwise")
      # Q_codon_array_vectored <- c(t(Q_codon_array)) # has to be transposed
      # Q_codon_array_vectored <- Q_codon_array_vectored[.non_zero_pos]
      anc.indices <- unique(phy.sort$edge[,1])
      if(verbose){ 
        MultiCoreLikelihoodBySite <- function(nsite.index){
          tmp <- selac:::GetLikelihoodSAC_CodonForSingleCharGivenOptimumHMMScoring(charnum=nsite.index, codon.data=codon.data$unique.site.patterns, 
                                                                           phy=phy.sort, Q_codon_array=Q_codon_array, 
                                                                           root.p=root.p_array, scale.factor=scale.factor, 
                                                                           anc.indices=anc.indices, return.all=FALSE)
          cat(".")
          return(tmp)
        }
        
      } else {
        MultiCoreLikelihoodBySite <- function(nsite.index){
          tmp <- selac:::GetLikelihoodSAC_CodonForSingleCharGivenOptimumHMMScoring(charnum=nsite.index, codon.data=codon.data$unique.site.patterns, phy=phy.sort, Q_codon_array_vectored=Q_codon_array_vectored, root.p=root.p_array, scale.factor=scale.factor, anc.indices=anc.indices, return.all=FALSE)
          return(tmp)
        }
      }
      final.likelihood.vector <- unlist(lapply(1:nsites.unique, MultiCoreLikelihoodBySite))
      if(verbose) cat("|\n")
      return(final.likelihood.vector)
    }
    
    if(include.gamma)
      model.params=c(model.params,std.gamma.shape)
    model.params = c(model.params,std.sel.reg)
    tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
    phy=tree; verbose=T;
    diploid=T; numcode=1
    x <- model.params
    k.levels =0 
    importance.of.aa.dist.in.selective.environment.change = 1
    
    rate.for.selective.environment.change = x[length(x)]
    x = x[-length(x)]
    aa.properties=NULL
    if(include.gamma == TRUE){
      shape = x[length(x)]
      x = x[-length(x)]
    }
    
    C.Phi.q.Ne <- x[1]
    C <- 4
    q <- 4e-7
    Ne <- 5e6
    Phi.q.Ne <- C.Phi.q.Ne / C
    Phi.Ne <- Phi.q.Ne / q
    Phi <- Phi.Ne / Ne
    alpha <- x[2]
    beta <- x[3]
    gamma <- 0.0003990333 #volume.fixed.value
    
    if(k.levels > 0){
      if(nuc.model == "JC") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        #During the early stages of the optimization process it will try weird values for the base frequencies.
        stopifnot(!any(base.freqs < 0))
        if(any(base.freqs < 0)){
          return(1000000)
        }
        nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        poly.params <- x[7:8]
      }
      if(nuc.model == "GTR") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        #During the early stages of the optimization process it will try weird values for the base frequencies.
        stopifnot(!any(base.freqs < 0))
        if(any(base.freqs < 0)){
          return(1000000)
        }
        nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
        poly.params <- x[7:8]
      }
      if(nuc.model == "UNREST") {
        nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model, base.freqs=NULL)
        poly.params <- x[4:5]
      }
    }else{
      if(nuc.model == "JC") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        #During the early stages of the optimization process it will try weird values for the base frequencies.
        stopifnot(!any(base.freqs < 0))
        if(any(base.freqs < 0)){
          return(1000000)
        }
        nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      }
      if(nuc.model == "GTR") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        #During the early stages of the optimization process it will try weird values for the base frequencies.
        stopifnot(!any(base.freqs < 0))
        if(any(base.freqs < 0)){
          return(1000000)
        }
        nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
      }
      if(nuc.model == "UNREST") {
        nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=NULL)
      }
    }
    
    nuc.mutation.rates.vector <- c(nuc.mutation.rates, rate.for.selective.environment.change)
    codon_mutation_matrix <- matrix(nuc.mutation.rates.vector[codon.index.matrix], dim(codon.index.matrix))
    codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
    nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
    nsites <- sum(codon.data$site.pattern.counts)
    
    if(include.gamma==TRUE){
      if(gamma.type == "median"){
        rates.k <- DiscreteGamma(shape=shape, ncats=ncats)
        weights.k <- rep(1/ncats, ncats)
      }
      if(gamma.type == "quadrature"){
        rates.and.weights <- LaguerreQuad(shape=shape, ncats=ncats)
        rates.k <- rates.and.weights[1:ncats]
        weights.k <- rates.and.weights[(ncats+1):(ncats*2)]
      }
      if(gamma.type == "lognormal"){
        rates.and.weights <- LogNormalQuad(shape=shape, ncats=ncats)
        rates.k <- rates.and.weights[1:ncats]
        weights.k <- rates.and.weights[(ncats+1):(ncats*2)]
      }
      #ttmmpp <- c(nuc.mutation.rates.vector, nsites.unique, nsites, C, Phi, rates.k, q, Ne, shape, importance.of.aa.dist.in.selective.environment.change)
      #writeLines(text = paste(ttmmpp), con = "~/Desktop/selac_parameter.txt", sep = "\t")
      final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites.unique)
      for(k.cat in sequence(ncats)){
        if(k.levels > 0){
          aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
        }else{
          aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
        }
        Q_codon_array <- selac:::FastCreateEvolveAACodonFixationProbabilityMatrix(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999, importance.of.aa.dist.in.selective.environment.change) #Cedric: added importance
        final.likelihood.mat[k.cat,] = GetLikelihoodSAC_CodonForManyCharVaryingBySiteEvolvingAA(codon.data, phy, Q_codon_array, 
          codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, codon_mutation_matrix=codon_mutation_matrix, 
          Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site, verbose=verbose)
      }
      likelihood <- sum(log(colSums(exp(final.likelihood.mat)*weights.k)) * codon.data$site.pattern.counts)
    }else{
      if(k.levels > 0){
        aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
      }else{
        aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
      }
      Q_codon_array <- selac:::FastCreateEvolveAACodonFixationProbabilityMatrix(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999, importance.of.aa.dist.in.selective.environment.change) #Cedric: added importance
      final.likelihood = GetLikelihoodSAC_CodonForManyCharVaryingBySiteEvolvingAA(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site, verbose=verbose)
      likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)
    }
    
    if(neglnl) {
      likelihood <- -1 * likelihood
    }
    if(verbose > 1) {
      results.vector <- c(likelihood, C*Phi*q, alpha, beta, gamma, Ne, ape::write.tree(phy))
      names(results.vector) <- c("likelihood", "C.Phi.q.Ne", "alpha", "beta", "gamma", "Ne", "phy")
      print(results.vector)
    }else if(verbose){
      results.vector <- c(likelihood, alpha, beta, gamma)
      names(results.vector) <- c("likelihood", "alpha", "beta", "gamma")
      print(results.vector)
      
    }
    if(is.na(likelihood) || is.nan(likelihood)){
      res <-1000000
    }else{
      res <- likelihood
    }
    result$loglik <-res
    # try({
    lSAC.c4mc.full(log(model.params), 
                   codon.data=codon.data, 
                   phy=tree,  
                   codon.freq.by.aa=NULL, 
                   codon.freq.by.gene=codon.freq.by.gene, 
                   numcode=1, 
                   diploid=TRUE, 
                   aa.properties=NULL, 
                   volume.fixed.value=std.gamma, 
                   nuc.model=nuc.model, 
                   codon.index.matrix=codon.index.matrix, 
                   include.gamma=include.gamma, 
                   gamma.type=gamma.type,
                   ncats=4, 
                   k.levels=0, 
                   logspace=TRUE, 
                   verbose=TRUE, 
                   n.cores.by.gene.by.site=nCores,
                   estimate.importance=FALSE) -> result$loglik
    # })
  }else if(mode=="FASTHMMTEST") {
    
    if(!file.exists("ecoli_chars.RData")){
      tmp.gene <- read.dna(fasta.file, format="fasta")
      tmp.gene <- as.list(as.matrix(cbind(tmp.gene)))
      
      chars <- selac:::DNAbinToCodonNumeric(tmp.gene)
      save(chars,file="ecoli_chars.RData")
    } else {
      load(file="ecoli_chars.RData")
    }
    codon.data <- chars[tree$tip.label,c(1,1+sample(ncol(chars)-1,10))]
    codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
    codon.data <- selac:::SitePattern(codon.data)
    
    
    
    codon.index.matrix = selac:::CreateCodonMutationMatrixIndexEvolveAA()
    model.params = hmm.params
    if(nuc.model != "UNREST")
      model.params=c(model.params,std.base.freq)
    model.params=c(model.params,std.nuc.params[[nuc.model]])
    
    lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA
    if(include.gamma)
      model.params=c(model.params,std.gamma.shape)
    model.params = c(model.params,std.sel.reg)
    # tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
    try({
      prof_obj <- profvis({
        lSAC.c4mc.full(log(model.params), 
                       codon.data=codon.data, 
                       phy=tree,  
                       codon.freq.by.aa=NULL, 
                       codon.freq.by.gene=codon.freq.by.gene, 
                       numcode=1, 
                       diploid=TRUE, 
                       aa.properties=NULL, 
                       volume.fixed.value=std.gamma, 
                       nuc.model=nuc.model, 
                       codon.index.matrix=codon.index.matrix, 
                       include.gamma=include.gamma, 
                       gamma.type=gamma.type,
                       ncats=4, 
                       k.levels=0, 
                       logspace=TRUE, 
                       verbose=TRUE, 
                       n.cores.by.gene.by.site=nCores,
                       estimate.importance=FALSE) -> result$loglik
      }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.05)
      save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
      # htmlwidgets::saveWidget(prof_obj, 
      #                         file=paste0(profile_prefix,".Rprofvis.html"))
    })
    
    
  }  else if(mode=="HMMEVAL50") {
    
    if(!file.exists("ecoli_chars.RData")){
      tmp.gene <- read.dna(fasta.file, format="fasta")
      tmp.gene <- as.list(as.matrix(cbind(tmp.gene)))
      
      chars <- selac:::DNAbinToCodonNumeric(tmp.gene)
      save(chars,file="ecoli_chars.RData")
    } else {
      load(file="ecoli_chars.RData")
    }
    codon.data <- chars[tree$tip.label,c(1,1+sample(ncol(chars)-1,50))]
    codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
    codon.data <- selac:::SitePattern(codon.data)
    
    
    
    codon.index.matrix = selac:::CreateCodonMutationMatrixIndexEvolveAA()
    model.params = hmm.params
    if(nuc.model != "UNREST")
      model.params=c(model.params,std.base.freq)
    model.params=c(model.params,std.nuc.params[[nuc.model]])
    
    lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA
    if(include.gamma)
      model.params=c(model.params,std.gamma.shape)
    model.params = c(model.params,std.sel.reg)
    # tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
    try({
      prof_obj <- profvis({
        lSAC.c4mc.full(log(model.params), 
                       codon.data=codon.data, 
                       phy=tree,  
                       codon.freq.by.aa=NULL, 
                       codon.freq.by.gene=codon.freq.by.gene, 
                       numcode=1, 
                       diploid=TRUE, 
                       aa.properties=NULL, 
                       volume.fixed.value=std.gamma, 
                       nuc.model=nuc.model, 
                       codon.index.matrix=codon.index.matrix, 
                       include.gamma=include.gamma, 
                       gamma.type=gamma.type,
                       ncats=4, 
                       k.levels=0, 
                       logspace=TRUE, 
                       verbose=TRUE, 
                       n.cores.by.gene.by.site=nCores,
                       estimate.importance=FALSE) -> result$loglik
      }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.05)
      save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
      # htmlwidgets::saveWidget(prof_obj, 
      #                         file=paste0(profile_prefix,".Rprofvis.html"))
    })
    
    
  } else if(mode=="HMMEVALFULL") {
    
    if(!file.exists("ecoli_chars.RData")){
      tmp.gene <- read.dna(fasta.file, format="fasta")
      tmp.gene <- as.list(as.matrix(cbind(tmp.gene)))
      
      chars <- selac:::DNAbinToCodonNumeric(tmp.gene)
      save(chars,file="ecoli_chars.RData")
    } else {
      load(file="ecoli_chars.RData")
    }
    codon.data <- chars[tree$tip.label,]
    codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
    codon.data <- selac:::SitePattern(codon.data)
    
    
    
    codon.index.matrix = selac:::CreateCodonMutationMatrixIndexEvolveAA()
    model.params = hmm.params
    if(nuc.model != "UNREST")
      model.params=c(model.params,std.base.freq)
    model.params=c(model.params,std.nuc.params[[nuc.model]])
    
    lSAC.c4mc.full <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA
    if(include.gamma)
      model.params=c(model.params,std.gamma.shape)
    model.params = c(model.params,std.sel.reg)
    # tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
    try({
      prof_obj <- profvis({
        lSAC.c4mc.full(log(model.params), 
                       codon.data=codon.data, 
                       phy=tree,  
                       codon.freq.by.aa=NULL, 
                       codon.freq.by.gene=codon.freq.by.gene, 
                       numcode=1, 
                       diploid=TRUE, 
                       aa.properties=NULL, 
                       volume.fixed.value=std.gamma, 
                       nuc.model=nuc.model, 
                       codon.index.matrix=codon.index.matrix, 
                       include.gamma=include.gamma, 
                       gamma.type=gamma.type,
                       ncats=4, 
                       k.levels=0, 
                       logspace=TRUE, 
                       verbose=TRUE, 
                       n.cores.by.gene.by.site=nCores,
                       estimate.importance=FALSE) -> result$loglik
      }, prof_output = paste0(profile_prefix,".Rprof"),interval=0.5)
      save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
      # htmlwidgets::saveWidget(prof_obj, 
      #                         file=paste0(profile_prefix,".Rprofvis.html"))
    })
    
    
  } else {
    cat(sprintf("Request for %s mode not understood.\n",as.character(mode)))
  }
  cat(sprintf("End: %s\tLL: %0.3f\n",
              profile_prefix,
              result$loglik))
  if(!file.exists(paste0(src.key,"_",codon.model,"_LL_log.csv")))
    cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
        file=paste0(src.key,"_",codon.model,"_LL_log.csv"),
        append = T)
  cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
              src.key,
              nuc.model,
              gamma.mode,
              selac_release,
              nCores,
              seed,
              result$loglik),
      file=paste0(src.key,"_",codon.model,"_LL_log.csv"),
      append = T)
  cat("SELAC Done. saving results\n")
  
  result$seed <- seed
  result$startingValues <- starting.vals
  result$startingTree <- tree
  
  save(result,file=sprintf('ecoli_output/%s_result.Rdata',profile_prefix))
  result$loglik
}


run_test_ecoli_optimize_no_profile <- function(seed=sample.int(1e6,1),ref="v1.6.1-rc1", nCores=1, auto.skip=T){
  setup_selac_for_profiling(ref=ref)
  src.key="ecoliDEBUG"
  nuc.model = 'UNREST'
  gamma.type="quadrature"
  profile_prefix=sprintf("%s_%s_%s_%s_%i_%i",
                         src.key,
                         nuc.model,
                         gamma.type,
                         selac_release,
                         nCores,
                         seed)
  if(file.exists(sprintf('ecoli_output/%s_result.Rdata',profile_prefix))){
    try({
    load(file=sprintf('ecoli_output/%s_result.Rdata',profile_prefix))
      if(!is.null(result$loglik) && is.finite(result$loglik)) {
        cat(sprintf("Skip: %s\n",
                    profile_prefix))
        return(result$loglik)
      } 
    })
    cat(sprintf("Rebuilding: %s\n",
                profile_prefix))
  }
  set.seed(seed)
  cat(sprintf("Start: %s\n",
              profile_prefix))
  tree<-read.tree('kosi07_data/kosi07_codonphyml_tree_TEM.newick')
  fasta.file="kosi07_data/aligned_KOSI07_TEM.fasta"
  output.file.name=sprintf('ecoli_output/%s_restart.Rdata',profile_prefix)
  result=list(loglik=NA)
  opt.aa.type <- "optimize"
  # random starting values
  starting.vals <- matrix(runif(n = 15, min = 0.01, max = 5), ncol = 15, nrow = 1)
  tree$edge.length <- runif(nrow(tree$edge), 0.01, 3)
  
  #nCores=3
  # try({
  #   prof_obj <- profvis({
      ## start.from.mle set to allow manual specification of fasta files
      # requires mle.matrix to be set, setting to start.from.mle==FALSE values for now
      # mle.matrix[1,] = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
      result <- SelacOptimize(codon.data.path = 'kosi07_data/', phy = tree, 
                              edge.length = 'optimize', optimal.aa = opt.aa.type, data.type='codon',
                              codon.model = 'selac', nuc.model = nuc.model, edge.linked=TRUE,
                              include.gamma = TRUE, gamma.type='quadrature', ncats = 4, numcode = 1,
                              diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = TRUE,
                              n.cores.by.gene  = nCores, n.cores.by.gene.by.site=1,
                              max.restarts = 1, max.evals=20, max.tol=1e-2, max.iterations = 15,
                              fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE,
                              output.restart.filename=output.file.name)
      # output.restart.filename=output.file.name, start.from.mle = TRUE,
      # mle.matrix=starting.vals, tol.step=1, partition.order = fasta.file)
    # }, prof_output = paste0(profile_prefix,".Rprof"),interval=1)
    # save(prof_obj, file=paste0(profile_prefix,".Rprofvis.RData"))
    # htmlwidgets::saveWidget(prof_obj, 
                            # file=paste0(profile_prefix,".Rprofvis.html"))
  # })
  cat(sprintf("End: %s\tLL: %0.3f\n",
              profile_prefix,
              result$loglik))
  if(!file.exists(paste0(src.key,"_LL_log.csv")))
    cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
        file=paste0(src.key,"_LL_log.csv"),
        append = T)
  cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
              src.key,
              nuc.model,
              gamma.type,
              selac_release,
              nCores,
              seed,
              result$loglik),
      file=paste0(src.key,"_LL_log.csv"),
      append = T)
  cat("SELAC Done. saving results\n")
  
  result$seed <- seed
  result$startingValues <- starting.vals
  result$startingTree <- tree
  
  save(result,file=sprintf('ecoli_output/%s_result.Rdata',profile_prefix))
  result$loglik
}
