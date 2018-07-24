#SelacOptimize extraction

SelacOptimizeDebug <- function(codon.data.path, 
                               n.partitions=NULL, 
                               phy, 
                               data.type="codon", 
                               codon.model="selac", 
                               edge.length="optimize", 
                               edge.linked=TRUE, 
                               optimal.aa="optimize", 
                               nuc.model="GTR", 
                               include.gamma=FALSE, 
                               gamma.type="quadrature", 
                               ncats=4, 
                               numcode=1, 
                               diploid=TRUE, 
                               k.levels=0, 
                               aa.properties=NULL, 
                               verbose=FALSE, 
                               n.cores.by.gene=1, n.cores.by.gene.by.site=1, 
                               max.tol=1e-3, max.tol.edges=1e-3, max.evals=1000000, max.restarts=3, 
                               user.optimal.aa=NULL, fasta.rows.to.keep=NULL, recalculate.starting.brlen=TRUE, 
                               output.by.restart=TRUE, output.restart.filename="restartResult", 
                               user.supplied.starting.param.vals=NULL, tol.step=1, 
                               optimizer.algorithm="NLOPT_LN_SBPLX", 
                               start.from.mle=FALSE, mle.matrix=NULL, partition.order=NULL, max.iterations=6) {
  
  if(!data.type == "codon" & !data.type == "nucleotide"){
    stop("Check that your data type input is correct. Options are codon or nucleotide", call.=FALSE)
  }
  if(!codon.model == "none" & !codon.model == "GY94" & !codon.model == "YN98" & !codon.model == "FMutSel0" & !codon.model == "FMutSel" & !codon.model == "selac"){
    stop("Check that your codon model is correct. Options are GY94, FMutSel0, or selac", call.=FALSE)
  }
  if(!edge.length == "optimize" & !edge.length == "fixed"){
    stop("Check that you have a supported edge length option. Options are optimize or fixed.", call.=FALSE)
  }
  if(!optimal.aa == "optimize" & !optimal.aa == "majrule" & !optimal.aa == "averaged" & !optimal.aa == "none" & !optimal.aa == "user"){
    stop("Check that you have a supported optimal amino acid option. Options are optimize, majrule, none, or user", call.=FALSE)
  }
  if(!nuc.model == "JC" & !nuc.model == "GTR" & !nuc.model == "UNREST"){
    stop("Check that you have a supported nucleotide substitution model. Options are JC, GTR, or UNREST.", call.=FALSE)
  }
  if(!gamma.type == "quadrature" & !gamma.type == "median" & !gamma.type == "lognormal"){
    stop("Check that you have a supported gamma type. Options are quadrature after Felsenstein 2001 or median after Yang 1994 or lognormal.", call.=FALSE)
  }
  
  if(!is.null(user.optimal.aa)){
    if(is.list(user.optimal.aa) == FALSE){
      stop("User-supplied optimal amino acids must be input as a list.", call.=FALSE)
    }
  }
  
  if(start.from.mle == TRUE){
    partitions <- partition.order
  }else{
    partitions <- system(paste("ls -1 ", codon.data.path, "*.fasta", sep=""), intern=TRUE)
  }
  
  if(is.null(n.partitions)){
    n.partitions <- length(partitions)
  }else{
    n.partitions = n.partitions
  }
  
  if(n.partitions<n.cores.by.gene) {
    warning(paste0("You have ", n.partitions, " partition (set with the n.partitions argument) but are asking to run across ", n.cores.by.gene, " cores, so ", n.cores.by.gene - n.partitions, " cores will not be used"))
  }
  
  cat(paste("Using", n.cores.by.gene * n.cores.by.gene.by.site, "total processors", sep=" "), "\n")
  
  cat("Initializing data and model parameters...", "\n")
  
  site.pattern.data.list <- as.list(numeric(n.partitions))
  site.pattern.count.list <- as.list(numeric(n.partitions))
  nsites.vector <- c()
  if(optimal.aa == "none"){
    if(data.type == "nucleotide"){
      empirical.base.freq.list <- as.list(numeric(n.partitions))
      starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
      for (partition.index in sequence(n.partitions)) {
        gene.tmp <- read.dna(partitions[partition.index], format='fasta')
        if(!is.null(fasta.rows.to.keep)){
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
        }else{
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
        }
        starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type=data.type, recalculate.starting.brlen=recalculate.starting.brlen)$edge.length
        nucleotide.data <- DNAbinToNucleotideNumeric(gene.tmp)
        nucleotide.data <- nucleotide.data[phy$tip.label,]
        nsites.vector = c(nsites.vector, dim(nucleotide.data)[2] - 1)
        empirical.base.freq <- as.matrix(nucleotide.data[,-1])
        empirical.base.freq <- table(empirical.base.freq, deparse.level = 0)/sum(table(empirical.base.freq, deparse.level = 0))
        empirical.base.freq.list[[partition.index]] <- as.vector(empirical.base.freq[1:4])
        nucleotide.data <- SitePattern(nucleotide.data, includes.optimal.aa=FALSE)
        site.pattern.data.list[[partition.index]] = nucleotide.data$unique.site.patterns
        site.pattern.count.list[[partition.index]] = nucleotide.data$site.pattern.counts
      }
    }else{
      codon.freq.by.gene.list <- as.list(numeric(n.partitions))
      empirical.aa.freq.list <- as.list(numeric(n.partitions))
      starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
      for (partition.index in sequence(n.partitions)) {
        gene.tmp <- read.dna(partitions[partition.index], format='fasta')
        if(!is.null(fasta.rows.to.keep)){
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
        }else{
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
        }
        starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type=data.type, recalculate.starting.brlen=recalculate.starting.brlen)$edge.length
        codon.data <- DNAbinToCodonNumeric(gene.tmp)
        codon.data <- codon.data[phy$tip.label,]
        nsites.vector = c(nsites.vector, dim(codon.data)[2] - 1)
        aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=numcode)
        aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
        empirical.aa.freq.list[[partition.index]] <- GetAAFreqsByGene(codon.data[,-1], aa.optim, numcode=numcode)
        codon.freq.by.gene.list[[partition.index]] <- GetCodonFreqsByGene(codon.data[,-1])
        codon.data <- SitePattern(codon.data, includes.optimal.aa=FALSE)
        site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
        site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
      }
    }
  }else{
    codon.freq.by.aa.list <- as.list(numeric(n.partitions))
    codon.freq.by.gene.list <- as.list(numeric(n.partitions))
    starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
    aa.optim.list <- as.list(numeric(n.partitions))
    aa.optim.full.list <- as.list(numeric(n.partitions))
    for (partition.index in sequence(n.partitions)) {
      gene.tmp <- read.dna(partitions[partition.index], format='fasta')
      if(!is.null(fasta.rows.to.keep)){
        gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
      }else{
        gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
      }
      starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type=data.type, recalculate.starting.brlen=recalculate.starting.brlen)$edge.length
      codon.data <- DNAbinToCodonNumeric(gene.tmp)
      codon.data <- codon.data[phy$tip.label,]
      nsites.vector = c(nsites.vector, dim(codon.data)[2] - 1)
      aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=numcode)
      if(optimal.aa == "user"){
        aa.optim <- user.optimal.aa[[partition.index]]
        aa.optim.full.list[[partition.index]] <- aa.optim
      }else{
        aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
        aa.optim.full.list[[partition.index]] <- aa.optim
      }
      codon.freq.by.aa.list[[partition.index]] <- GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=numcode)
      codon.freq.by.gene.list[[partition.index]] <- GetCodonFreqsByGene(codon.data[,-1])
      aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
      colnames(aa.optim.frame.to.add) <- colnames(codon.data)
      codon.data <- rbind(codon.data, aa.optim.frame.to.add)
      codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
      site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
      site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
      aa.optim.list[[partition.index]] = codon.data$optimal.aa
    }
  }
  
  opts <- list("algorithm" = optimizer.algorithm, "maxeval" = max.evals, "ftol_rel" = max.tol)
  opts.edge <- list("algorithm" = optimizer.algorithm, "maxeval" = max.evals, "ftol_rel" = max.tol.edges)
  
  
  results.final <- c()
  if(nuc.model == "JC"){
    nuc.ip = NULL
    max.par.model.count = 0
    parameter.column.names <- c()
  }
  if(nuc.model == "GTR"){
    nuc.ip = rep(1, 5)
    max.par.model.count = 5
    parameter.column.names <- c("C_A", "G_A", "T_A", "G_C", "T_C")
  }
  if(nuc.model == "UNREST"){
    nuc.ip = rep(1, 11)
    max.par.model.count = 11
    parameter.column.names <- c("C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
  }
  
  if(optimal.aa=="none") {
    if(data.type == "nucleotide"){
      codon.index.matrix = NA
      if(include.gamma == TRUE){
        ip = c(1,nuc.ip)
        upper = c(5, rep(21, length(ip)-1))
        lower = rep(-21, length(ip))
        max.par.model.count = max.par.model.count + 1
        parameter.column.names <- c("shape.gamma", parameter.column.names)
      }else{
        ip = nuc.ip
        upper = rep(21, length(ip))
        lower = rep(-21, length(ip))
      }
      index.matrix = matrix(0, n.partitions, length(ip))
      index.matrix[1,] = 1:ncol(index.matrix)
      ip.vector = ip
      upper.vector = upper
      lower.vector = lower
      if(n.partitions > 1){
        for(partition.index in 2:n.partitions){
          ip.vector = c(ip.vector, ip)
          upper.vector = c(upper.vector, upper)
          lower.vector = c(lower.vector, lower)
          index.matrix.tmp = numeric(max.par.model.count)
          index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
          index.matrix[partition.index,] <- index.matrix.tmp
        }
      }
      number.of.current.restarts <- 1
      best.lik <- 1000000
      while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing analysis...", sep=""), "\n")
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        if(edge.length == "optimize"){
          cat("       Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #opts.edge <- opts
          upper.edge <- rep(log(10), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("       Optimizing model parameters", "\n")
        ParallelizedOptimizedByGene <- function(n.partition){
          optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsLarge, ub=upper.vector[1:dim(mle.pars.mat)[2]], lb=lower.vector[1:dim(mle.pars.mat)[2]], opts=opts, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
          tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          return(tmp.pars)
        }
        results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
        results.final <- NULL
        results.final$objective <- sum(parallelized.parameters[,1])
        results.final$solution <- c(t(parallelized.parameters[,-1]))
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
        print(results.final$objective)
        print(mle.pars.mat)
        
        current.likelihood <- results.final$objective
        cat(paste("Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number<=max.iterations){
          cat(paste("Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          if(edge.length == "optimize"){
            cat("       Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            
            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("       Optimizing model parameters", "\n")
          opts.params <- opts
          opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
          
          ParallelizedOptimizedByGene <- function(n.partition){
            optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsLarge, ub=upper.vector[1:dim(mle.pars.mat)[2]], lb=lower.vector[1:dim(mle.pars.mat)[2]], opts=opts.params, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
            tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
            return(tmp.pars)
          }
          results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
          results.final <- NULL
          results.final$objective <- sum(parallelized.parameters[,1])
          results.final$solution <- c(t(parallelized.parameters[,-1]))
          mle.pars.mat <- index.matrix
          mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
          print(results.final$objective)
          print(mle.pars.mat)
          
          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp,file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
        }
        number.of.current.restarts <- number.of.current.restarts + 1
      }
      
      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
      cat("Finished. Summarizing results...", "\n")
      colnames(mle.pars.mat) <- parameter.column.names
      
      if(edge.length == "optimize"){
        np <- max(index.matrix) + length(phy$edge.length)
      }else{
        np <- max(index.matrix)
      }
      obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, data.type=data.type, codon.model=codon.model, nsites=nsites.vector, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.evals=max.evals)
      class(obj) = "selac"
    }else{
      if(codon.model == "GY94"){
        max.par.model.count <- 2
        ip = c(1,1)
        parameter.column.names <- c("V", "kappa")
        upper = rep(log(99), length(ip))
        lower = rep(-21, length(ip))
        
        codon.index.matrix = NA
        
        index.matrix = matrix(0, n.partitions, length(ip))
        index.matrix[1,] = 1:ncol(index.matrix)
        ip.vector = ip
        upper.vector = upper
        lower.vector = lower
        if(n.partitions > 1){
          for(partition.index in 2:n.partitions){
            #ip.vector = c(ip.vector, 1)
            #upper.vector = c(upper.vector, log(99))
            #lower.vector = c(lower.vector, -10)
            #index.matrix.tmp = numeric(max.par.model.count)
            #index.matrix.tmp[2] = 2
            #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
            #index.matrix[partition.index,] <- index.matrix.tmp
            index.matrix[partition.index,] <- 1:ncol(index.matrix)
          }
        }
      }
      if(codon.model == "YN98"){
        max.par.model.count <- 2
        ip = c(1,1)
        parameter.column.names <- c("omega", "kappa")
        upper = rep(log(99), length(ip))
        lower = rep(-21, length(ip))
        
        codon.index.matrix = NA
        
        index.matrix = matrix(0, n.partitions, length(ip))
        index.matrix[1,] = 1:ncol(index.matrix)
        ip.vector = ip
        upper.vector = upper
        lower.vector = lower
        if(n.partitions > 1){
          for(partition.index in 2:n.partitions){
            #ip.vector = c(ip.vector, 1)
            #upper.vector = c(upper.vector, log(99))
            #lower.vector = c(lower.vector, -10)
            #index.matrix.tmp = numeric(max.par.model.count)
            #index.matrix.tmp[2] = 2
            #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
            #index.matrix[partition.index,] <- index.matrix.tmp
            index.matrix[partition.index,] <- 1:ncol(index.matrix)
          }
        }
      }
      if(codon.model == "FMutSel0"){
        empirical.aa.freq.unlist <- matrix(unlist(empirical.aa.freq.list), ncol = 21, byrow = TRUE)
        empirical.aa.freq <- colSums(empirical.aa.freq.unlist)/ sum(colSums(empirical.aa.freq.unlist))
        fitness.pars <- GetFitnessStartingValues(codon.freqs=empirical.aa.freq)[-c(17,21)]
        aa.ordered <- .unique.aa
        aa.ordered <- aa.ordered[-c(17,21)]
        if(nuc.model == "UNREST"){
          max.par.model.count <- max.par.model.count + 1 + 19
          ip = c(0.4, nuc.ip, fitness.pars)
          parameter.column.names <- c("omega", parameter.column.names, paste("fitness", aa.ordered, sep="_"))
          upper = rep(log(99), length(ip))
          lower = rep(-10, length(ip))
        }else{
          max.par.model.count <- max.par.model.count + 3 + 1 + 19
          ip = c(0.4, .25, .25, .25, nuc.ip, fitness.pars)
          parameter.column.names <- c("omega", "freqA", "freqC", "freqG", parameter.column.names, paste("fitness", aa.ordered, sep="_"))
          upper = c(log(99), 0, 0, 0, rep(log(99), length(ip)-4))
          lower = rep(-10, length(ip))
        }
        
        codon.index.matrix = NA
        
        index.matrix = matrix(0, n.partitions, length(ip))
        index.matrix[1,] = 1:ncol(index.matrix)
        ip.vector = ip
        upper.vector = upper
        lower.vector = lower
        if(n.partitions > 1){
          for(partition.index in 2:n.partitions){
            #ip.vector = c(ip.vector, 0.4)
            #upper.vector = c(upper.vector, log(99))
            #lower.vector = c(lower.vector, -10)
            #index.matrix.tmp = numeric(max.par.model.count)
            #index.matrix.tmp[2:max.par.model.count] = 2:max.par.model.count
            #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
            #index.matrix[partition.index,] <- index.matrix.tmp
            index.matrix[partition.index,] <- 1:ncol(index.matrix)
          }
        }
      }
      if(codon.model == "FMutSel"){
        empirical.codon.freq.unlist <- matrix(unlist(codon.freq.by.gene.list), ncol = 64, byrow = TRUE)
        empirical.codon.freq <- colSums(empirical.codon.freq.unlist)/ sum(colSums(empirical.codon.freq.unlist))
        fitness.pars <- GetFitnessStartingValues(codon.freqs=empirical.codon.freq, n.pars=64)
        codon.ordered <- .codon.name
        codon.ordered <- codon.ordered[-c(49,51,57,64)]
        if(nuc.model == "UNREST"){
          max.par.model.count <- max.par.model.count + 1 + 60
          ip = c(0.4, nuc.ip, fitness.pars)
          parameter.column.names <- c("omega", parameter.column.names, paste("fitness", codon.ordered, sep="_"))
          upper = c(rep(log(99), length(ip)-3))
          lower = rep(-10, length(ip))
        }else{
          max.par.model.count <- max.par.model.count + 3 + 1 + 60
          ip = c(0.4, .25, .25, .25, nuc.ip, fitness.pars)
          parameter.column.names <- c("omega", "freqA", "freqC", "freqG", parameter.column.names, paste("fitness", codon.ordered, sep="_"))
          upper = c(log(99), 0, 0, 0, rep(log(99), length(ip)-4))
          lower = rep(-10, length(ip))
        }
        
        codon.index.matrix = NA
        
        index.matrix = matrix(0, n.partitions, length(ip))
        index.matrix[1,] = 1:ncol(index.matrix)
        ip.vector = ip
        upper.vector = upper
        lower.vector = lower
        if(n.partitions > 1){
          for(partition.index in 2:n.partitions){
            #ip.vector = c(ip.vector, 0.4)
            #upper.vector = c(upper.vector, log(99))
            #lower.vector = c(lower.vector, -10)
            #index.matrix.tmp = numeric(max.par.model.count)
            #index.matrix.tmp[2:max.par.model.count] = 2:max.par.model.count
            #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
            #index.matrix[partition.index,] <- index.matrix.tmp
            index.matrix[partition.index,] <- 1:ncol(index.matrix)
          }
        }
      }
      
      number.of.current.restarts <- 1
      best.lik <- 10000000
      while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing analysis...", sep=""), "\n")
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        print(mle.pars.mat)
        if(edge.length == "optimize"){
          cat("       Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #opts.edge <- opts
          upper.edge <- rep(log(10), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("       Optimizing model parameters", "\n")
        #ParallelizedOptimizedByGene <- function(n.partition){
        optim.by.gene <- nloptr(x0=log(ip.vector), eval_f = OptimizeModelParsLarge, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
        #tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
        #return(tmp.pars)
        #}
        #results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        #results.set <- lapply(1:n.partitions, ParallelizedOptimizedByGene)
        #parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
        #results.final <- NULL
        #results.final$objective <- sum(parallelized.parameters[,1])
        #results.final$solution <- c(t(parallelized.parameters[,-1]))
        results.final$objective <- optim.by.gene$objective
        results.final$solution <- optim.by.gene$solution
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
        print(results.final$objective)
        print(mle.pars.mat)
        
        current.likelihood <- results.final$objective
        cat(paste("Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number<=max.iterations){
          cat(paste("Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          if(edge.length == "optimize"){
            cat("       Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            
            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("       Optimizing model parameters", "\n")
          opts.params <- opts
          opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
          print(length(results.final$solution))
          print(results.final$solution)
          #ParallelizedOptimizedByGene <- function(n.partition){
          optim.by.gene <- nloptr(x0=results.final$solution, eval_f = OptimizeModelParsLarge, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
          #tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          # return(tmp.pars)
          #}
          #results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          #results.set <- lapply(1:n.partitions, ParallelizedOptimizedByGene)
          #parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
          #results.final <- NULL
          #results.final$objective <- sum(parallelized.parameters[,1])
          #results.final$solution <- c(t(parallelized.parameters[,-1]))
          results.final$objective <- optim.by.gene$objective
          results.final$solution <- optim.by.gene$solution
          mle.pars.mat <- index.matrix
          mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
          print(results.final$objective)
          print(mle.pars.mat)
          
          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, empirical.aa.freqs=empirical.aa.freq.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp,file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
        }
        number.of.current.restarts <- number.of.current.restarts + 1
      }
      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
      cat("Finished. Summarizing results...", "\n")
      colnames(mle.pars.mat) <- parameter.column.names
      
      if(edge.length == "optimize"){
        np <- max(index.matrix) + length(phy$edge.length)
      }else{
        np <- max(index.matrix)
      }
      obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, data.type=data.type, codon.model=codon.model, nsites=nsites.vector, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, empirical.aa.freqs=empirical.aa.freq.list, max.tol=max.tol, max.evals=max.evals)
      class(obj) = "selac"
    }
  }
  if(optimal.aa=="majrule" | optimal.aa=="optimize" | optimal.aa=="averaged" | optimal.aa=="user") {
    codon.index.matrix = CreateCodonMutationMatrixIndex()
    cpv.starting.parameters <- GetAADistanceStartingParameters(aa.properties=aa.properties)
    if(max.restarts > 1){
      selac.starting.vals <- matrix(0, max.restarts+1, 3)
      selac.starting.vals[,1] <- runif(n = max.restarts+1, min = (10^-10)*5e6, max = (10^-5)*5e6)
      selac.starting.vals[,2] <- runif(n = max.restarts+1, min = 0.01, max = 3)
      selac.starting.vals[,3] <- runif(n = max.restarts+1, min = 0.01, max = 1)
    }else{
      if(is.null(user.supplied.starting.param.vals)){
        selac.starting.vals <- matrix(c(2, 1.8292716544, 0.1017990371), 1, 3)
        selac.starting.vals <- rbind(selac.starting.vals, c(2, 1.8292716544, 0.1017990371))
      }else{
        selac.starting.vals <- matrix(c(user.supplied.starting.param.vals[1], user.supplied.starting.param.vals[2], user.supplied.starting.param.vals[3]), 1, 3)
        selac.starting.vals <- rbind(selac.starting.vals, c(user.supplied.starting.param.vals[1], user.supplied.starting.param.vals[2], user.supplied.starting.param.vals[3]))
      }
    }
    if(include.gamma == TRUE){
      #Gamma variation is turned ON:
      if(nuc.model == "JC"){
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "shape.gamma")
          upper = c(log(50),  21, 21, 0, 0, 0, 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 1
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "shape.gamma")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10, 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 2 + 1
        }
      }
      if(nuc.model == "GTR") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma")
          upper = c(log(50), 21, 21, 0, 0, 0, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 1
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 2	+ 1
        }
      }
      if(nuc.model == "UNREST") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma")
          upper = c(log(50), 21, 21, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 1
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma")
          upper = c(log(50), 21, 21, 10, 10, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 2 + 1
        }
      }
      index.matrix = matrix(0, n.partitions, max.par.model.count)
      index.matrix[1,] = 1:ncol(index.matrix)
      ip.vector = ip
      if(n.partitions > 1){
        # Gamma variation is turned ON:
        for(partition.index in 2:n.partitions){
          if(nuc.model == "JC"){
            #ip.vector = c(ip.vector, ip[1], ip[8])
            if(start.from.mle == TRUE){
              ip.vector = c(ip.vector, mle.matrix[partition.index,1])
            }else{
              ip.vector = c(ip.vector, ip[1])
            }
          }else{
            if(nuc.model == "GTR"){
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                #index.matrix.tmp[c(2:11)] = c(2:11)
                #ip.vector = c(ip.vector, ip[1], ip[12])
                index.matrix.tmp[c(2:12)] = c(2:12)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, mle.matrix[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }else{
                #index.matrix.tmp[c(2:13)] = c(2:13)
                #ip.vector = c(ip.vector, ip[1], ip[14])
                index.matrix.tmp[c(2:14)] = c(2:14)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, mle.matrix[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }
            }else{
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                #index.matrix.tmp[c(2:14)] = c(2:14)
                #ip.vector = c(ip.vector, ip[1], ip[15])
                index.matrix.tmp[c(2:15)] = c(2:15)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, mle.matrix[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }else{
                #index.matrix.tmp[c(2:16)] = c(2:16)
                #ip.vector = c(ip.vector, ip[1], ip[17])
                index.matrix.tmp[c(2:17)] = c(2:17)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, ip[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }
            }
          }
          index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
          index.matrix[partition.index,] <- index.matrix.tmp
        }
      }
    }else{
      # Gamma variation is turned OFF:
      if(nuc.model == "JC"){
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25)
          }
          parameter.column.names <- c("C.q.phi.Ne",  "alpha", "beta", "freqA", "freqC", "freqG")
          upper = c(log(50), 21, 21, 0, 0, 0)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 0
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 0 + 2
        }
      }
      if(nuc.model == "GTR") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C")
          upper = c(log(50), 21, 21, 0, 0, 0, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 0
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 0 + 2
        }
      }
      if(nuc.model == "UNREST") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], nuc.ip)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
          upper = c(log(50), 21, 21, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 0
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
          upper = c(log(50), 21, 21, 10, 10, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 0 + 2
        }
      }
      index.matrix = matrix(0, n.partitions, max.par.model.count)
      index.matrix[1,] = 1:ncol(index.matrix)
      ip.vector = ip
      if(n.partitions > 1){
        for(partition.index in 2:n.partitions){
          if(nuc.model == "JC"){
            if(start.from.mle == TRUE){
              ip.vector = c(ip.vector, mle.matrix[partition.index,1])
            }else{
              ip.vector = c(ip.vector, ip[1])
            }
          }else{
            if(nuc.model == "GTR"){
              if(start.from.mle == TRUE){
                ip.vector = c(ip.vector, mle.matrix[partition.index,1])
              }else{
                ip.vector = c(ip.vector, ip[1])
              }
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                index.matrix.tmp[c(2:11)] = c(2:11)
              }else{
                index.matrix.tmp[c(2:13)] = c(2:13)
              }
            }else{
              if(start.from.mle == TRUE){
                ip.vector = c(ip.vector, mle.matrix[partition.index,1])
              }else{
                ip.vector = c(ip.vector, ip[1])
              }
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                index.matrix.tmp[c(2:14)] = c(2:14)
              }else{
                index.matrix.tmp[c(2:16)] = c(2:16)
              }
            }
          }
          index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
          index.matrix[partition.index,] <- index.matrix.tmp
        }
      }
    }
    
    #THIS IS FOR THERE IS A SEPARATE GAMMA PER GENE:
    #if(include.gamma == TRUE){
    #    index.matrix.red <- t(matrix(1:(n.partitions*2), 2, n.partitions))
    #}else{
    #    index.matrix.red <- t(matrix(1:n.partitions, 1, n.partitions))
    #}
    
    #This is so we can break out alpha, beta, GTR, and gamma which are shared among ALL genes:
    index.matrix.red <- t(matrix(1:n.partitions, 1, n.partitions))
    
    if(optimal.aa == "optimize"){
      number.of.current.restarts <- 1
      aa.optim.original <- aa.optim.list
      best.lik <- 1000000
      while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing random restart ", number.of.current.restarts,"...", sep=""), "\n")
        aa.optim.list <- aa.optim.original
        cat("       Doing first pass using majority-rule optimal amino acid...", "\n")
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        print(mle.pars.mat)
        if(edge.length == "optimize"){
          cat("              Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #phy$edge.length <- colMeans(starting.branch.lengths) / (1/selac.starting.vals[number.of.current.restarts, 2])
          #opts.edge <- opts
          upper.edge <- rep(log(50), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("              Optimizing model parameters", "\n")
        
        #if(include.gamma == TRUE){
        #    alpha.beta.gtr <- mle.pars.mat[1,c(2:(max.par.model.count-1))]
        #    upper.bounds.shared <- upper[c(2:(max.par.model.count-1))]
        #    lower.bounds.shared <- lower[c(2:(max.par.model.count-1))]
        #}else{
        alpha.beta.gtr <- mle.pars.mat[1,c(2:max.par.model.count)]
        upper.bounds.shared <- upper[c(2:max.par.model.count)]
        lower.bounds.shared <- lower[c(2:max.par.model.count)]
        #}
        
        ParallelizedOptimizedByGene <- function(n.partition){
          #if(include.gamma == TRUE){
          #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
          #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
          #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
          #}else{
          tmp.par.mat <- as.matrix(mle.pars.mat[,1])
          upper.bounds.gene <- upper[1]
          lower.bounds.gene <- lower[1]
          #}
          optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          return(tmp.pars)
        }
        results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        #if(include.gamma == TRUE){
        #The number of columns is 3: [1] log-likelihood, [2] C.q.phi, [3] phi gamma:
        #parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
        #}else{
        #The number of columns is 2: [1] log-likelihood, [2] C.q.phi:
        parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
        #}
        results.final <- NULL
        results.final$objective <- sum(parallelized.parameters[,1])
        results.final$solution <- c(t(parallelized.parameters[,-1]))
        mle.pars.mat.red <- index.matrix.red
        mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
        print(mle.pars.mat.red)
        optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
        results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
        alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
        #if(include.gamma == TRUE){
        #    mle.pars.mat <- c()
        #    for(row.index in 1:dim(mle.pars.mat.red)[1]){
        #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
        #    }
        #}else{
        mle.pars.mat <- c()
        for(row.index in 1:dim(mle.pars.mat.red)[1]){
          mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
        }
        #}
        print(results.final$objective)
        print(mle.pars.mat)
        
        current.likelihood <- results.final$objective
        cat(paste("       Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number <= max.iterations){
          cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          cat("              Optimizing amino acids", "\n")
          aa.optim.list <- as.list(numeric(n.partitions))
          ParallelizedOptimizeAAByGene <- function(n.partition){
            gene.tmp <- read.dna(partitions[n.partition], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            codon.data <- DNAbinToCodonNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            tmp.aa.optim.full <- GetOptimalAAPerSite(x=log(mle.pars.mat[n.partition,]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, neglnl=TRUE, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
            return(tmp.aa.optim.full)
          }
          aa.optim.full.list <- mclapply(1:n.partitions, ParallelizedOptimizeAAByGene, mc.cores=n.cores.by.gene)
          for(partition.index in sequence(n.partitions)) {
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            codon.data <- DNAbinToCodonNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            codon.freq.by.aa.list[[partition.index]] <- GetCodonFreqsByAA(codon.data[,-1], aa.optim.full.list[[partition.index]], numcode=numcode)
            aa.optim.frame.to.add <- matrix(c("optimal", aa.optim.full.list[[partition.index]]), 1, dim(codon.data)[2])
            colnames(aa.optim.frame.to.add) <- colnames(codon.data)
            codon.data <- rbind(codon.data, aa.optim.frame.to.add)
            codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
            site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
            site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
            aa.optim.list[[partition.index]] = codon.data$optimal.aa
          }
          if(edge.length == "optimize"){
            cat("              Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            
            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("              Optimizing model parameters", "\n")
          
          ParallelizedOptimizedByGene <- function(n.partition){
            #if(include.gamma == TRUE){
            #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
            #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
            #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
            #}else{
            tmp.par.mat <- as.matrix(mle.pars.mat[,1])
            upper.bounds.gene <- upper[1]
            lower.bounds.gene <- lower[1]
            #}
            opts.params <- opts
            opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f=OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts.params, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
            return(tmp.pars)
          }
          results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          #if(include.gamma == TRUE){
          #The number of columns is 3: [1] log-likelihood, [2] C.q.phi.Ne, [3] phi gamma:
          #parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
          #}else{
          #The number of columns is 2: [1] log-likelihood, [2] C.q.phi.Ne:
          parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
          #}
          results.final <- NULL
          results.final$objective <- sum(parallelized.parameters[,1])
          results.final$solution <- c(t(parallelized.parameters[,-1]))
          mle.pars.mat.red <- index.matrix.red
          mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
          
          optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
          alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
          #if(include.gamma == TRUE){
          #    mle.pars.mat <- c()
          #    for(row.index in 1:dim(mle.pars.mat.red)[1]){
          #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
          #    }
          #}else{
          mle.pars.mat <- c()
          for(row.index in 1:dim(mle.pars.mat.red)[1]){
            mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
          }
          #}
          print(results.final$objective)
          print(mle.pars.mat)
          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=aa.optim.full.list, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp, file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
          best.aa.optim.list <- aa.optim.full.list
          best.codon.freq.by.aa <- codon.freq.by.aa.list
          best.codon.freq.by.gene <- codon.freq.by.gene.list
        }
        number.of.current.restarts <- number.of.current.restarts + 1
        ip.vector[c(index.matrix[,1])] <- selac.starting.vals[number.of.current.restarts, 1]
        ip.vector[2:3] <- selac.starting.vals[number.of.current.restarts, 2:3]
        aa.optim.list <- aa.optim.original
      }
      selac.starting.vals <- best.ip
      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      aa.optim.full.list <- best.aa.optim.list
      codon.freq.by.aa.list <- best.codon.freq.by.aa
      codon.freq.by.gene.list <- best.codon.freq.by.gene
      
      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
    }else{
      if(optimal.aa == "averaged"){
        aa.optim.list = NULL
      }
      number.of.current.restarts <- 1
      best.lik <- 1000000
      while(number.of.current.restarts < (max.restarts+1)){
        if(optimal.aa == "user"){
          cat(paste("Finished. Performing random restart ", number.of.current.restarts," using user-supplied optimal amino acids...", sep=""), "\n")
        }else{
          cat(paste("Finished. Performing random restart ", number.of.current.restarts," using majority-rule optimal amino acids...", sep=""), "\n")
        }
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        cat("       Doing first pass...", "\n")
        print(mle.pars.mat)
        if(edge.length == "optimize"){
          cat("              Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #phy$edge.length <- colMeans(starting.branch.lengths) / (1/selac.starting.vals[number.of.current.restarts, 2])
          #opts.edge <- opts
          upper.edge <- rep(log(50), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("              Optimizing model parameters", "\n")
        
        #if(include.gamma == TRUE){
        #    alpha.beta.gtr <- mle.pars.mat[1,c(2:(max.par.model.count-1))]
        #    upper.bounds.shared <- upper[c(2:(max.par.model.count-1))]
        #    lower.bounds.shared <- lower[c(2:(max.par.model.count-1))]
        #}else{
        alpha.beta.gtr <- mle.pars.mat[1,c(2:max.par.model.count)]
        upper.bounds.shared <- upper[c(2:max.par.model.count)]
        lower.bounds.shared <- lower[c(2:max.par.model.count)]
        #}
        
        ParallelizedOptimizedByGene <- function(n.partition){
          #if(include.gamma == TRUE){
          #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
          #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
          #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
          #}else{
          tmp.par.mat <- as.matrix(mle.pars.mat[,1])
          upper.bounds.gene <- upper[1]
          lower.bounds.gene <- lower[1]
          #}
          optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          return(tmp.pars)
        }
        results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        #if(include.gamma == TRUE){
        #The number of columns is 3: [1] log-likelihood, [2] C.q.phi.Ne, [3] phi gamma:
        #    parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
        #}else{
        #The number of columns is 2: [1] log-likelihood, [2] C.q.phi.Ne:
        parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
        #}
        results.final <- NULL
        results.final$objective <- sum(parallelized.parameters[,1])
        results.final$solution <- c(t(parallelized.parameters[,-1]))
        mle.pars.mat.red <- index.matrix.red
        mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
        print(mle.pars.mat.red)
        optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
        results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
        alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
        #if(include.gamma == TRUE){
        #    mle.pars.mat <- c()
        #    for(row.index in 1:dim(mle.pars.mat.red)[1]){
        #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
        #    }
        #}else{
        mle.pars.mat <- c()
        for(row.index in 1:dim(mle.pars.mat.red)[1]){
          mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
        }
        #}
        print(results.final$objective)
        print(mle.pars.mat)
        
        current.likelihood <- results.final$objective
        cat(paste("       Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number <= max.iterations){
          cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          if(edge.length == "optimize"){
            cat("              Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            
            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("              Optimizing model parameters", "\n")
          
          ParallelizedOptimizedByGene <- function(n.partition){
            #if(include.gamma == TRUE){
            #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
            #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
            #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
            #}else{
            tmp.par.mat <- as.matrix(mle.pars.mat[,1])
            upper.bounds.gene <- upper[1]
            lower.bounds.gene <- lower[1]
            #}
            opts.params <- opts
            opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f=OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts.params, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
            return(tmp.pars)
          }
          results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          #if(include.gamma == TRUE){
          #The number of columns is 3: [1] log-likelihood, [2] C.q.phi, [3] phi gamma:
          #    parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
          #}else{
          #The number of columns is 2: [1] log-likelihood, [2] C.q.phi:
          parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
          #}
          results.final <- NULL
          results.final$objective <- sum(parallelized.parameters[,1])
          results.final$solution <- c(t(parallelized.parameters[,-1]))
          mle.pars.mat.red <- index.matrix.red
          mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
          
          optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
          alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
          #if(include.gamma == TRUE){
          #   mle.pars.mat <- c()
          #   for(row.index in 1:dim(mle.pars.mat.red)[1]){
          #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
          #    }
          #}else{
          mle.pars.mat <- c()
          for(row.index in 1:dim(mle.pars.mat.red)[1]){
            mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
          }
          #}
          print(results.final$objective)
          print(mle.pars.mat)
          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=aa.optim.full.list, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp, file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
          best.codon.freq.by.aa <- codon.freq.by.aa.list
          best.codon.freq.by.gene <- codon.freq.by.gene.list
        }
        number.of.current.restarts <- number.of.current.restarts + 1
        print(ip.vector)
        ip.vector[c(index.matrix[,1])] <- selac.starting.vals[number.of.current.restarts, 1]
        ip.vector[2:3] <- selac.starting.vals[number.of.current.restarts, 2:3]
        print(ip.vector)
      }
      selac.starting.vals <- best.ip
      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      codon.freq.by.aa.list <- best.codon.freq.by.aa
      codon.freq.by.gene.list <- best.codon.freq.by.gene
      
      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
    }
    cat("Finished. Summarizing results...", "\n")
    colnames(mle.pars.mat) <- parameter.column.names
    
    if(edge.length == "optimize"){
      if(optimal.aa == "user" | optimal.aa == "majrule" | optimal.aa == "averaged"){
        np <- max(index.matrix) + length(phy$edge.length)
      }else{
        np <- max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)
      }
    }else{
      if(optimal.aa == "user" | optimal.aa == "majrule" | optimal.aa == "averaged"){
        np <- max(index.matrix)
      }else{
        np <- max(index.matrix) + sum(nsites.vector)
      }
    }
    
    #Counting parameters: Do we count the nsites too? Yup.
    obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=aa.optim.full.list, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=selac.starting.vals)
    class(obj) = "selac"
  }
  return(obj)
}


