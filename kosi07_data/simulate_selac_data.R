rm(list=ls())
library(selac)


phy <- read.tree('../data/selac.yeast.tre')

set.size <- 10

gtr.params <- runif(6, 0.1, 2)
acgt.freq <- runif(4)
acgt.freq <- (acgt.freq / sum(acgt.freq))[-4]

numcode <- 1
codon.sets <- selac:::CreateCodonSets()
codon.set.translate <- apply(codon.sets, 2, n2s)
codon.name <- apply(codon.set.translate, 1, paste, collapse="")
codon.aa <- sapply(codon.name,selac:::TranslateCodon, numcode=numcode)
aa.names <- unique(codon.aa)
root.codon.freq.matrix <- matrix(0, nrow=length(aa.names), ncol=length(codon.aa))
rownames(root.codon.freq.matrix) <- aa.names
colnames(root.codon.freq.matrix) <- names(codon.aa)
for (aa.index in sequence(dim(root.codon.freq.matrix)[1])) {
	matching.codons <- which(codon.aa == rownames(root.codon.freq.matrix)[aa.index])
	root.codon.freq.matrix[aa.index, matching.codons] <- 1/length(matching.codons)
}

for(g in 1:set.size)
{
  cat("simulating gene ", g, "\n")
  gene.length <- round(runif(1, 100, 1000))
  optim.aa <- sample(x = c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"), 
                     size = gene.length, replace = T)
  
  #They are ordered as follows: C.q.phi, alpha, beta, Ne, base.freqs for A C G, and the rates for the nucleotide model.
  model.param <- c(runif(1, 10^-10, 10^-5), runif(2, 0.01, 1), 57852168, acgt.freq, gtr.params)
  names(model.param) <- c("qphic", "alpha", "beta", "Ne", "A.freq", "C.freq", "G.freq", 
                          "GTR.alpha", "GTR.beta", "GTR.gamma", "GTR.delta", "GTR.epsilon", "GTR.eta")
  
  
  nucl.data <- SelacSimulator(phy, pars = model.param, aa.optim_array = optim.aa, 
                 root.codon.array = root.codon.freq.matrix, numcode=1, aa.properties=NULL, 
                 nuc.model = "GTR", k.levels=0, diploid=TRUE)
  
  sim.name <- paste("../data/simulated/simSet_cedric/gene", g, sep = "") 
  write.dna(x = nucl.data, file = paste(sim.name, "fasta", sep = "."), format = "fasta", nbcol = 1, colw = 60)
  save(list = c("model.param", "optim.aa", "root.codon.freq.matrix", "phy"), file = paste(sim.name, "rda", sep = "."))
}
