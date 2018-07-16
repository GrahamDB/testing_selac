rm(list=ls())

seed <- as.numeric(Sys.time())
set.seed(seed)
ncores.per.site <- 1
ncores.per.gene <- 1

run <- 1

# read tree
tree <- read.tree(file = 'kosi07_codonphyml_tree_TEM.newick')

# get path to data
codon.data.path <- "aln/"
fasta.files <- c("aln/aligned_KOSI07_TEM.fasta")
out.path <- "results/"

opt.aa.type <- "optimize"
# random starting values
starting.vals <- matrix(runif(n = 15, min = 0.01, max = 5), ncol = 15, nrow = 1)
tree$edge.length <- runif(nrow(tree$edge), 0.01, 3)


output.file.name <- paste(out.path, "TEM_", opt.aa.type, "_wG_", run, ".Rda", sep="")
# run selac
cat("Setup Done. Starting SELAC\n")
start <- Sys.time()
result <- SelacOptimize(codon.data.path, n.partitions=1, phy=tree, data.type="codon", codon.model="selac", edge.length="optimize", 
                        edge.linked=TRUE, optimal.aa=opt.aa.type, nuc.model="UNREST", include.gamma=TRUE, gamma.type="quadrature", 
                        ncats=4, numcode=2, diploid=TRUE, k.levels=0, aa.properties=NULL, verbose=FALSE, n.cores.by.gene=ncores.per.gene, 
                        n.cores.by.gene.by.site=ncores.per.site, max.tol=1e-2, max.evals=10000, max.initial.cond=1, max.iterations = 15,
                        fasta.rows.to.keep=NULL, recalculate.starting.brlen=FALSE, output.by.restart=FALSE, conv.crit = 0.01,
                        output.restart.filename=output.file.name, start.from.mle = TRUE,
                        mle.matrix=starting.vals, tol.step=1, partition.order = fasta.files)
end <- Sys.time()
end - start
cat("SELAC Done. saving results\n")

result$seed <- seed
result$startingValues <- starting.vals
result$startingTree <- tree
save(list = "result", file = output.file.name)
