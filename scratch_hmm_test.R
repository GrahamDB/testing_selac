#tmp to build FASTHMM and FASTHMMTEST
#base

if(F){
  # nCores = 1
  
  # For nuc.model=="GTR"
  # and gamma.mode=="none"
  codon.model="selac"
  include.gamma = F
  gamma.type = "quadrature"
  
  # HMM code requires starting edge length < 0.5 and > 1e-8
  tree$edge.length <- runif(nrow(tree$edge), 0.01, 0.45)
  result <- SelacHMMOptimize(codon.data.path = 'kosi07_data/',
                             phy = tree, 
                             edge.length = 'optimize',
                             data.type='codon',
                             codon.model = codon.model, 
                             nuc.model = nuc.model, 
                             edge.linked=TRUE,
                             include.gamma = include.gamma, 
                             gamma.type=gamma.type, 
                             ncats = 4, 
                             numcode = 1,
                             diploid = TRUE, 
                             k.levels = 0, 
                             aa.properties = NULL, 
                             verbose = FALSE,
                             n.cores.by.gene = 1, 
                             n.cores.by.gene.by.site=nCores,
                             max.restarts = 1, 
                             max.evals=1, 
                             max.tol=1e-2,
                             fasta.rows.to.keep=NULL, 
                             recalculate.starting.brlen=FALSE, 
                             output.by.restart=FALSE,
                             output.restart.filename=output.file.name)
}
if(F){
  std.params = c(C.Phi.q.Ne = 4*4e-7*.5*5e6,
                 alpha=1.829272, 
                 beta=0.101799)
  
  hmm.params = c(C.Phi.q.Ne = 2,
                 alpha=1.829272, 
                 beta=0.101799)
  std.gamma=0.0003990333
  std.base.freq = c(A=0.25,C=0.25,G=0.25)
  std.poly.params = c(NA,NA)
  std.gamma.shape = 5
  std.sel.reg = 0.01
}
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