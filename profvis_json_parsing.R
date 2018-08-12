#profvis json parsing functions
library(foreach)
library(lattice)
library(iterators)

# str(prof.json.data <- read_json("selacFULL_GTR_noneXquadrature_v1.6.1-rc1_3_830124.Rprofvis.json",simplifyVector = T), list.len = 4 )
# prof.data<- as.data.frame(prof.json.data$x$message$prof)
# length(unique(with(prof.data,paste(label,filenum,linenum,sep="#"))))
# length(unique(prof.data[prof.data$meminc!=0,"time"]))

read_profvis_json <- function(file.name){
  library(jsonlite)
  prof.json.data <- read_json("selacFULL_GTR_noneXquadrature_v1.6.1-rc1_3_830124.Rprofvis.json",simplifyVector = T)
  list(prof_data=as.data.frame(prof.json.data$x$message$prof),
       interval=prof.json.data$x$message$interval,
       files=prof.json.data$x$message$files)
}

label_key <- function(prof.data){
  with(prof.data,paste(label,filename,linenum,sep="#"))
}
file_key <- function(prof.data){
  na.omit(unique(prof.data[,c("filenum","filename")]))
}

meminc.ts <- function(prof.data){ 
  key=unique(prof.data$time)
  res <- data.frame(time=key,mem.total=0,mem.delta=0,mem.alloc=0,mem.release=0,row.names = paste0("T",key))
  alloc=unique(prof.data[,c("time","memalloc")])
  inc=prof.data[prof.data$meminc!=0,c("time","meminc")]
  res[ paste0("T",alloc$time),"mem.total"]=alloc$memalloc
  res[ paste0("T",inc$time),"mem.delta"]=inc$meminc
  res$mem.alloc[res$mem.delta>0]=res$mem.delta[res$mem.delta>0]
  res$mem.release[res$mem.delta<0]=-res$mem.delta[res$mem.delta<0]
  
  return(res)
}

flatten_recursion <- function(prof.data){
  res <- unique(prof.data[,c("time","label","filename","linenum")])
  res$time_key= paste0("T",res$time)
    res$label_key<-label_key(res)
  res
}

summ_prof.data <- function(prof.data){
  blocks <- flatten_recursion(prof.data)
  mem <- meminc.ts(prof.data )
  label_data <- unique(prof.data[,c("label","filename","linenum","filenum")])
  row.names(label_data)<-label_key(label_data)
  label_data$t.count=0
  label_data$mem.alloc=0
  label_data$mem.release=0
  
  tmp <- tapply(X = blocks$time_key,blocks$label_key,length)
  label_data[names(tmp),"t.count"]=tmp
  tmp <- tapply(X = mem[blocks$time_key,"mem.alloc"],blocks$label_key,sum)
  label_data[names(tmp),"mem.alloc"]=tmp
  tmp <- tapply(X = mem[blocks$time_key,"mem.release"],blocks$label_key,sum)
  label_data[names(tmp),"mem.release"]=tmp
  label_data$t.rel = label_data$t.count / nrow(mem)
  label_data$mem.alloc.rel = label_data$mem.alloc / sum(mem$mem.alloc)
  label_data$mem.release.rel = label_data$mem.release / sum(mem$mem.release)
  
  label_data
}

rprof_proc <- function(raw_file){
  profile.summary <- list(func=custum_summaryRprof(filename = raw_file,memory = "both",lines = "hide",basenames = 3),
                          lines=custum_summaryRprof(filename = raw_file,memory = "both",lines = "show",basenames = 3),
                          mixed=custum_summaryRprof(filename = raw_file,memory = "both",lines = "both",basenames = 3))
  save(profile.summary,file=paste0(raw_file,"_summary.RData"))
  return(profile.summary$func$by.self)
}
# dir(pattern="selacFULLb_GTR_noneXquadrature_744c6c8_1_[0-9]*.Rprof")

unionALL <- function(...){
  items<-list(...)
  len<-length(items)
  
  if(len==0) return(logical(0))
  if(len==1) return(items[[1]])
  if(len==2) return(union(items[[1]],items[[2]]))
  union(do.call(unionALL,items[1:(len%/%2)]),
        do.call(unionALL,items[-(1:(len%/%2))]))
}

intersectALL <- function(...){
  items<-list(...)
  len<-length(items)
  
  if(len==0) return(logical(0))
  if(len==1) return(items[[1]])
  if(len==2) return(intersect(items[[1]],items[[2]]))
  intersect(do.call(intersectALL,items[1:(len%/%2)]),
        do.call(intersectALL,items[-(1:(len%/%2))]))
}

custum_summaryRprof <-function (filename = "Rprof.out", chunksize = 5000, 
                                memory = c("none", "both", "tseries", "stats"), 
                                lines = c("hide", "show", "both"), 
                                index = 2, diff = TRUE, exclude = NULL, basenames = 1) {
  con <- file(filename, "rt")
  on.exit(close(con))
  firstline <- readLines(con, n = 1L)
  if (!length(firstline)) 
    stop(gettextf("no lines found in %s", sQuote(filename)), 
         domain = NA)
  sample.interval <- as.numeric(strsplit(firstline, "=")[[1L]][2L])/1e+06
  memory.profiling <- substr(firstline, 1L, 6L) == "memory"
  line.profiling <- grepl("line profiling", firstline)
  if (line.profiling) 
    filenames <- character(0)
  memory <- match.arg(memory)
  if (memory != "none" && !memory.profiling) 
    stop("profile does not contain memory information")
  if (memory == "tseries") 
    return(Rprof_memory_summary(filename = con, chunksize = chunksize, 
                                label = index, diff = diff, exclude = exclude, sample.interval = sample.interval))
  else if (memory == "stats") 
    return(Rprof_memory_summary(filename = con, chunksize = chunksize, 
                                aggregate = index, diff = diff, exclude = exclude, 
                                sample.interval = sample.interval))
  lines <- match.arg(lines)
  if (lines != "hide" && !line.profiling) 
    stop("profile does not contain line information")
  fnames <- NULL
  ucounts <- NULL
  fcounts <- NULL
  memcounts <- NULL
  umem <- NULL
  repeat ({
    chunk <- readLines(con, n = chunksize)
    if (line.profiling) {
      filenamelines <- grep("^#File [0-9]+: ", chunk)
      if (length(filenamelines)) {
        fnum <- as.integer(sub("^#File ([0-9]+): .*", 
                               "\\1", chunk[filenamelines]))
        filenames[fnum] <- sub("^#File [0-9]+: ", "", 
                               chunk[filenamelines])
        if (basenames) {
          dirnames <- dirname(filenames[fnum])
          filenames[fnum] <- basename(filenames[fnum])
          for (i in seq_len(basenames - 1)) {
            tail <- basename(dirnames)
            filenames[fnum] <- ifelse(tail == ".", filenames[fnum], 
                                      paste0(tail, "/", filenames[fnum]))
            parent <- dirname(dirnames)
            dirnames <- ifelse(dirnames == parent, ".", 
                               parent)
          }
        }
        chunk <- chunk[-filenamelines]
      }
    }
    if (length(chunk) == 0L) 
      break
    if (memory.profiling) {
      memprefix <- attr(regexpr(":[0-9]+:[0-9]+:[0-9]+:[0-9]+:", 
                                chunk), "match.length")
      if (memory == "both") {
        memstuff <- substr(chunk, 2L, memprefix - 1L)
        memcounts <- pmax(apply(sapply(strsplit(memstuff, 
                                                ":"), as.numeric), 1, diff), 0)
        if (!is.matrix(memcounts)) 
          memcounts <- matrix(memcounts, nrow = 1)
        memcounts <- c(0, rowSums(cbind(memcounts[, 
                                                  1L:2L, drop = FALSE] * 8, memcounts[, 3L, 
                                                                                      drop = FALSE])))
        rm(memstuff)
      }
      chunk <- substr(chunk, memprefix + 1L, nchar(chunk, 
                                                   "c"))
      if (any((nc <- nchar(chunk, "c")) == 0L)) {
        chunk <- chunk[nc > 0L]
        memcounts <- memcounts[nc > 0L]
      }
    }
    chunk <- strsplit(chunk, " ")
    if (line.profiling) 
      chunk <- lapply(chunk, function(x) {
        locations <- !startsWith(x, "\"")
        # if (lines != "hide") {
          fnum <- sub("#.*", "", x[locations])
          lnum <- sub(".*#", "", x[locations])
          x[locations] <- paste0(filenames[as.integer(fnum)], 
                                 "#", lnum)
        # }
        if(lines != "show" && any(which(locations)>1)){
          lab_locations=which(locations[-1])
          x[lab_locations] = paste0(x[lab_locations+1],"#",x[lab_locations])
        }
        
        switch(lines, hide = x <- x[!locations], show = x <- x[locations])
        if (length(x)) 
          x
        else "<no location>"
      })
    newfirsts <- sapply(chunk, "[[", 1L)
    newuniques <- lapply(chunk, unique)
    ulen <- lengths(newuniques)
    newuniques <- unlist(newuniques)
    new.utable <- table(newuniques)
    new.ftable <- table(factor(newfirsts, levels = names(new.utable)))
    if (memory == "both") 
      new.umem <- rowsum(memcounts[rep.int(seq_along(memcounts), 
                                           ulen)], newuniques)
    fcounts <- rowsum(c(as.vector(new.ftable), fcounts), 
                      c(names(new.ftable), fnames))
    ucounts <- rowsum(c(as.vector(new.utable), ucounts), 
                      c(names(new.utable), fnames))
    if (memory == "both") 
      umem <- rowsum(c(new.umem, umem), c(names(new.utable), 
                                          fnames))
    fnames <- sort(unique(c(fnames, names(new.utable))))
  })
  firstnum <- fcounts * sample.interval
  uniquenum <- ucounts * sample.interval
  index1 <- order(-firstnum, -uniquenum)
  index2 <- order(-uniquenum, -firstnum)
  if (lines == "show") {
    filename <- sub("#.*$", "", fnames)
    linenum <- rep(0, length(filename))
    hasline <- filename != fnames
    linenum[hasline] <- as.numeric(sub("^.*#", "", fnames[hasline]))
    index3 <- order(filename, linenum)
  }
  firstpct <- round(100 * firstnum/sum(firstnum), 2)
  uniquepct <- round(100 * uniquenum/sum(firstnum), 2)
  digits <- ifelse(sample.interval < 0.01, 3L, 2L)
  firstnum <- round(firstnum, digits)
  uniquenum <- round(uniquenum, digits)
  if (memory == "both") 
    memtotal <- round(umem/1048576, 1)
  rval <- data.frame(firstnum, firstpct, uniquenum, uniquepct)
  names(rval) <- c("self.time", "self.pct", "total.time", 
                   "total.pct")
  rownames(rval) <- fnames
  if (memory == "both") 
    rval$mem.total <- memtotal
  by.self <- rval[index1, ]
  by.self <- by.self[by.self[, 1L] > 0, ]
  by.total <- rval[index2, c(3L, 4L, if (memory == "both") 5L, 
                             1L, 2L)]
  result <- list(by.self = by.self, by.total = by.total)
  if (lines == "show") 
    result <- c(result, list(by.line = rval[index3, ]))
  c(result, sample.interval = sample.interval, sampling.time = sum(fcounts) * 
      sample.interval)
}

comparison.lattice.layer <- 
  function(src.key, nuc.model, gamma.type, selac_release, nCores, seed,
           min.self.pct=NULL, min.self.time=NULL, 
           min.total.pct=NULL, min.total.time=NULL) {
    profile_prefix=sprintf("%s_%s_%s_%s_%i_%i",
                           src.key, nuc.model, gamma.type,
                           selac_release, nCores, seed)
    
    if(!file.exists(paste0(profile_prefix,".Rprof")) ||
       grepl(pattern = "^ecoli", x=src.key) 
       && !file.exists(paste0("ecoli_output/",profile_prefix,"_result.Rdata"))){
      res <- numeric(0)
      dim(res) <- c(0,5,1,1,1,1)
      dimnames(res) <- list(label=character(0),
                            measure=c("self.time", "self.pct",
                                      "total.time", "total.pct", 
                                      "mem.total"),
                            src=paste(src.key,nuc.model,gamma.type,sep="_"),
                                            ver=selac_release, mc=paste0("p",nCores),
                                            seed=paste0("S",seed))
      return(res)
    }
    if(file.exists(paste0("ecoli_output/",profile_prefix,"_result.Rdata"))){
      load(paste0("ecoli_output/",profile_prefix,"_result.Rdata"))
      if(is.null(result) || !is.finite(result$loglik)){
        res <- numeric(0)
        dim(res) <- c(0,5,1,1,1,1)
        dimnames(res) <- list(label=character(0),
                              measure=c("self.time", "self.pct",
                                        "total.time", "total.pct", 
                                        "mem.total"),
                              src=paste(src.key,nuc.model,gamma.type,sep="_"),
                              ver=selac_release, mc=paste0("p",nCores),
                              seed=paste0("S",seed))
        return(res)
      }
    }
    if(file.exists(paste0(profile_prefix,"_Rprofsummary.RData")) &&
       (file.mtime(paste0(profile_prefix,"_Rprofsummary.RData")) > file.mtime(paste0(profile_prefix,".Rprof")) )  ){
      load(file=paste0(profile_prefix,"_Rprofsummary.RData"))
    } else {
      rprof_summary <- NULL
      try({
      rprof_summary <- custum_summaryRprof(paste0(profile_prefix,".Rprof"),
                                           lines="hide", basenames = 3, memory = "both")[["by.self"]]
      save(rprof_summary, file=paste0(profile_prefix,"_Rprofsummary.RData"))
      })
      if(is.null(rprof_summary)){
        res <- numeric(0)
        dim(res) <- c(0,5,1,1,1,1)
        dimnames(res) <- list(label=character(0),
                              measure=c("self.time", "self.pct",
                                        "total.time", "total.pct", 
                                        "mem.total"),
                              src=paste(src.key,nuc.model,gamma.type,sep="_"),
                              ver=selac_release, mc=paste0("p",nCores),
                              seed=paste0("S",seed))
        return(res)
      }
    }
    if(!is.null(min.self.time))
      rprof_summary <- rprof_summary[rprof_summary$self.time>=min.self.time,]
    if(!is.null(min.self.pct))
      rprof_summary <- rprof_summary[rprof_summary$self.pct>=min.self.pct,]
    if(!is.null(min.total.time))
      rprof_summary <- rprof_summary[rprof_summary$total.time>=min.total.time,]
    if(!is.null(min.total.pct))
      rprof_summary <- rprof_summary[rprof_summary$total.pct>=min.total.pct,]
    tmp<-dimnames(rprof_summary)
    rprof_summary <-as.matrix(rprof_summary)
    dim(rprof_summary) <- c(dim(rprof_summary),rep(1,4))
    names(tmp)=c("label","measure")
    tmp$label <- backtrack_equivalent_version(tmp$label)
    dimnames(rprof_summary) <- c(tmp,list(src=paste(src.key,nuc.model,gamma.type,sep="_"),
                                          ver=selac_release, mc=paste0("p",nCores),
                                          seed=paste0("S",seed)))
    rprof_summary
  }
library(parallel)
comparison.lattice <- 
  function(src.key, nuc.model, gamma.type, selac_release, nCores, seed, ...){
    keys=expand.grid(src.key, nuc.model, gamma.type, selac_release, nCores, seed)
    # print(keys)
    args=list(...)
    mcmapply(comparison.lattice.layer,
           src.key=keys[,1], nuc.model=keys[,2], gamma.type=keys[,3], selac_release=keys[,4], nCores=keys[,5], seed=keys[,6],
           MoreArgs = args, SIMPLIFY = FALSE,mc.cores=4,mc.preschedule = F) -> res_values
    dim_set <- sapply(c("label","measure","src","ver","mc","seed"),
                      function(y)  do.call(unionALL, 
                                           lapply(res_values,
                                                  function(x) dimnames(x)[[y]])),simplify=F)
      #do.call(unionALL,lapply(res,function(x) dimnames(x)[["label"]] ))
    res = rep(NA,do.call(prod,lapply(dim_set,length)))
    dim(res) <- sapply(dim_set,length)
    dimnames(res) <- dim_set
    for(res_layer in res_values)
      res[as.matrix(do.call(expand.grid,dimnames(res_layer)))]<-res_layer
    
    return(res)
  }


## comparing accross versions
#asterion@gdb438:~/git/selac$ git diff -r "ca13013" -r "744c6c8" | grep '^@@\|^diff --git'
#diff --git a/DESCRIPTION b/DESCRIPTION
#@@ -1,8 +1,8 @@
#diff --git a/R/selac.R b/R/selac.R
#@@ -1214,27 +1214,44 @@ GetLikelihoodSAC_CodonForManyCharVaryingBySite <- function(codon.data, phy, Q_co
#@@ -3215,6 +3232,98 @@ GetMaxName <- function(x) {
#@@ -3240,7 +3349,7 @@ GetExpQt <- function(phy, Q, scale.factor, rates=NULL){
#asterion@gdb438:~/git/selac$ git diff -r "ca13013" -r "dd94866" | grep '^@@\|^diff --git'
#diff --git a/DESCRIPTION b/DESCRIPTION
#@@ -1,8 +1,8 @@
#diff --git a/R/selac.R b/R/selac.R
#@@ -1120,21 +1120,18 @@ GetLikelihoodSAC_CodonForSingleCharGivenOptimum <- function(charnum=1, codon.dat
#@@ -1214,27 +1211,44 @@ GetLikelihoodSAC_CodonForManyCharVaryingBySite <- function(codon.data, phy, Q_co
#@@ -3215,6 +3229,98 @@ GetMaxName <- function(x) {
#@@ -3240,7 +3346,7 @@ GetExpQt <- function(phy, Q, scale.factor, rates=NULL){
#asterion@gdb438:~/git/selac$ git diff -r "744c6c8" -r "dd94866" | grep '^@@\|^diff --git'
#diff --git a/R/selac.R b/R/selac.R
#@@ -1120,21 +1120,18 @@ GetLikelihoodSAC_CodonForSingleCharGivenOptimum <- function(charnum=1, codon.dat

# 744c6c8/R/selac.R#[[1-1122]] == dd94866/R/selac.R#[[1-1122]]       :0
# 744c6c8/R/selac.R#[[1138-]] == dd94866/R/selac.R#[[1135-]]         :3
# ca13013/R/selac.R#[[1-1122]] == dd94866/R/selac.R#[[1-1122]]
# ca13013/R/selac.R#[[1138-1216]] == dd94866/R/selac.R#[[1135-1213]]
# ca13013/R/selac.R#[[1238-3217]] == dd94866/R/selac.R#[[1252-3231]]
# ca13013/R/selac.R#[[3218-3242]] == dd94866/R/selac.R#[[3324-3348]]
# ca13013/R/selac.R#[[3244-]] == dd94866/R/selac.R#[[3350-]]
# ca13013/R/selac.R#[[1-1216]] == 744c6c8/R/selac.R#[[1-1216]]       :0
# ca13013/R/selac.R#[[1238-3217]] == 744c6c8/R/selac.R#[[1255-3234]] :-17
# ca13013/R/selac.R#[[3218-3242]] == 744c6c8/R/selac.R#[[3327-3351]] :-109
# ca13013/R/selac.R#[[3244-]] == 744c6c8/R/selac.R#[[3353-]]         :-109
backtrack_equivalent_version <- function(label_names){
  res <- label_names
  # print(res)
  check.key<-which(grepl("#",label_names) & grepl("^GrahamDB",label_names) )
  if(length(check.key) == 0) return(res)
  if(length(check.key) == 1) {
    tmp <- strsplit(label_names[check.key],"#")[[1]]
    label_data <- data.frame(tmp[1],as.integer(tmp[2]),tmp[3],stringsAsFactors=F)
  }else{
  label_data <- as.data.frame(t(sapply(strsplit(label_names[check.key],"#"), as.character)),stringsAsFactors=F)
  label_data[[2]]<-as.integer(label_data[[2]])
  }
  # print(label_data)
  tmp <-grepl("selac.R",label_data[[1]]) 
  label_data[[1]][!tmp] <- sub("^GrahamDB-selac-.......","GrahamDB-selac-ca13013",label_data[[1]][!tmp])
  
  #e90938f..ab3e84e
  #a ab3e84e to e90938f, file ab3e84e/R/selac.R
  tmp <-grepl("ab3e84e/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1229 & tmp
  slot2 <- label_data[[2]] >= 1231 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-ab3e84e","GrahamDB-selac-e90938f",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 0
  
  #ca4a198..e90938f
  #a e90938f to ca4a198, file e90938f/R/selac.R
  tmp <-grepl("e90938f/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1227 & tmp
  slot2 <- label_data[[2]] >= 1235 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-e90938f","GrahamDB-selac-ca4a198",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + -6

  #f1c3d1e..ca4a198
  #a ca4a198 to f1c3d1e, file ca4a198/R/selac.R
  tmp <-grepl("ca4a198/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1221 & tmp
  slot2 <- label_data[[2]] >= 1225 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-ca4a198","GrahamDB-selac-f1c3d1e",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + -2
  
  #86cc7e9..f1c3d1e
  #a f1c3d1e to 86cc7e9, file f1c3d1e/R/selac.R
  tmp <-grepl("f1c3d1e/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1134 & tmp
  slot2 <- label_data[[2]] >= 1136 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-f1c3d1e","GrahamDB-selac-86cc7e9",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 0
  
  #5eec8b6..86cc7e9
  #a 86cc7e9 to 5eec8b6, file 86cc7e9/R/selac.R
  tmp <-grepl("86cc7e9/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1119 & tmp
  slot2 <- label_data[[2]] >= 1124 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-86cc7e9","GrahamDB-selac-5eec8b6",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 0
  
  #2a6c393..5eec8b6
  #a 5eec8b6 to 2a6c393, file 5eec8b6/R/selac.R
  tmp <-grepl("5eec8b6/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1118 & tmp
  slot2 <- label_data[[2]] >= 1120 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-5eec8b6","GrahamDB-selac-2a6c393",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 0
  
  #ade1e58..2a6c393
  #a 2a6c393 to ade1e58, file 2a6c393/R/selac.R
  tmp <-grepl("2a6c393/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1077 & tmp
  slot2 <- label_data[[2]] >= 1079 & label_data[[2]] <= 1102 & tmp
  slot3 <- label_data[[2]] >= 1144 & label_data[[2]] <= 1206 & tmp
  slot4 <- label_data[[2]] >= 1214 & tmp
  slots <- slot1|slot2|slot3|slot4
  label_data[[1]][slots] <- sub("^GrahamDB-selac-2a6c393","GrahamDB-selac-ade1e58",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 0
  label_data[[2]][slot3] <- label_data[[2]][slot3] + -38
  label_data[[2]][slot4] <- label_data[[2]][slot4] + -38
  
  #5c98a1f..ade1e58
  #a ade1e58 to 5c98a1f, file ade1e58/R/selac.R
  tmp <-grepl("ade1e58/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1187 & tmp
  slot2 <- label_data[[2]] >= 1189 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-ade1e58","GrahamDB-selac-5c98a1f",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 0
  
  #9cc4412..5c98a1f
  #a 5c98a1f to 9cc4412, file 5c98a1f/R/selac.R
  tmp <-grepl("5c98a1f/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1453 & tmp
  slot2 <- label_data[[2]] >= 1463 & label_data[[2]] <= 1583 & tmp
  slot3 <- label_data[[2]] >= 1593 & label_data[[2]] <= 1662 & tmp
  slot4 <- label_data[[2]] >= 1672 & label_data[[2]] <= 1702 & tmp
  slot5 <- label_data[[2]] >= 1712 & label_data[[2]] <= 1763 & tmp
  slot6 <- label_data[[2]] >= 1773 & tmp
  slots <- slot1|slot2|slot3|slot4|slot5|slot6
  label_data[[1]][slots] <- sub("^GrahamDB-selac-5c98a1f","GrahamDB-selac-9cc4412",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + -5
  label_data[[2]][slot3] <- label_data[[2]][slot3] + -10
  label_data[[2]][slot4] <- label_data[[2]][slot4] + -15
  label_data[[2]][slot5] <- label_data[[2]][slot5] + -20
  label_data[[2]][slot6] <- label_data[[2]][slot6] + -25
  
  #000ff0d..9cc4412
  #a 9cc4412 to 000ff0d, file 9cc4412/R/selac.R
  tmp <-grepl("9cc4412/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 5122 & tmp
  slot2 <- label_data[[2]] >= 5124 & label_data[[2]] <= 5541 & tmp
  slot3 <- label_data[[2]] >= 5543 & tmp
  slots <- slot1|slot2|slot3
  label_data[[1]][slots] <- sub("^GrahamDB-selac-9cc4412","GrahamDB-selac-000ff0d",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 0
  label_data[[2]][slot3] <- label_data[[2]][slot3] + 0
  
  #e94b31e..000ff0d
  #a 000ff0d to e94b31e, file 000ff0d/R/selac.R
  tmp <-grepl("000ff0d/R/selac.R",label_data[[1]])
  slot1 <- label_data[[2]] >= 0 & label_data[[2]] <= 1084 & tmp
  slot2 <- label_data[[2]] >= 1103 & tmp
  slots <- slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-000ff0d","GrahamDB-selac-e94b31e",label_data[[1]][slots])
  label_data[[2]][slot1] <- label_data[[2]][slot1] + 0
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 7
  
  tmp <-grepl("e94b31e/R/selac.R",label_data[[1]]) 
  slot1 <- label_data[[2]] <= 1146 & tmp
  slot1 <- label_data[[2]] >= 1148 & label_data[[2]] <= 1178 & tmp | slot1
  slot2 <- label_data[[2]] >= 1197 & label_data[[2]] <= 1443 & tmp
  slot2 <- label_data[[2]] >= 1445 & label_data[[2]] <= 1453 & tmp | slot2
  slot2 <- label_data[[2]] >= 1455 & tmp | slot2
  slots=slot1|slot2
  label_data[[1]][slots] <- sub("^GrahamDB-selac-e94b31e","GrahamDB-selac-9e90977",label_data[[1]][slots])
  label_data[[2]][slot2] <- label_data[[2]][slot2] - 13
  
  tmp <-grepl("9e90977/R/selac.R",label_data[[1]]) 
  slot1 <- label_data[[2]] <= 3358 & tmp
  slot2 <- label_data[[2]] >= 3363 & label_data[[2]] <= 3385 & tmp
  slot3 <- label_data[[2]] >= 3393 & tmp
  slots=slot1|slot2|slot3
  label_data[[1]][slots] <- sub("^GrahamDB-selac-9e90977","GrahamDB-selac-8f13a5a",label_data[[1]][slots])
  label_data[[2]][slot2] <- label_data[[2]][slot2] - 4
  label_data[[2]][slot3] <- label_data[[2]][slot3] - 6
  
  tmp <-grepl("8f13a5a/R/selac.R",label_data[[1]]) 
  slot1 <- label_data[[2]] <= 3443 & tmp
  slot2 <- label_data[[2]] >= 3452 & tmp
  label_data[[1]][slot1 | slot2] <- sub("^GrahamDB-selac-8f13a5a","GrahamDB-selac-eb6213a",label_data[[1]][slot1 | slot2])
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 2
  
  tmp <-grepl("eb6213a/R/selac.R",label_data[[1]]) 
  slot1 <- label_data[[2]] <= 3233 & tmp
  slot2 <- label_data[[2]] >= 3422 & tmp
  label_data[[1]][slot1 | slot2] <- sub("^GrahamDB-selac-eb6213a","GrahamDB-selac-38a4c36",label_data[[1]][slot1 | slot2])
  label_data[[2]][slot2] <- label_data[[2]][slot2] - 101
  
  tmp <-grepl("38a4c36/R/selac.R",label_data[[1]]) 
  slot1 <- label_data[[2]] <= 1136 & tmp
  slot2 <- label_data[[2]] >= 1139 & label_data[[2]] <= 3242 & tmp
  slot3 <- label_data[[2]] >= 3244 & label_data[[2]] <= 3252 & tmp
  slot4 <- label_data[[2]] >= 3254 & tmp
  slots=slot1|slot2|slot3|slot4
  label_data[[1]][slots] <- sub("^GrahamDB-selac-38a4c36","GrahamDB-selac-dd94866",label_data[[1]][slots])
  label_data[[2]][slot2|slot3|slot4] <- label_data[[2]][slot2|slot3|slot4] -1
  
  tmp <-grepl("dd94866/R/selac.R",label_data[[1]]) 
  slot1 <- label_data[[2]] <= 1122 & tmp
  slot2 <- label_data[[2]] >= 1135 & tmp
  label_data[[1]][slot1 | slot2] <- sub("^GrahamDB-selac-dd94866","GrahamDB-selac-744c6c8",label_data[[1]][slot1 | slot2])
  label_data[[2]][slot2] <- label_data[[2]][slot2] + 3
  
  tmp <-grepl("744c6c8/R/selac.R",label_data[[1]]) 
  slot1 <- label_data[[2]] <= 1216 & tmp
  slot2 <- label_data[[2]] >= 1255 & label_data[[2]] <= 3234 & tmp
  slot3 <- label_data[[2]] >= 3324 & label_data[[2]] <= 3351 & tmp
  slot4 <- label_data[[2]] >= 3353 & tmp
  label_data[[1]][slot1|slot2|slot3|slot4] <- 
    sub("^GrahamDB-selac-744c6c8","GrahamDB-selac-ca13013",
        label_data[[1]][slot1|slot2|slot3|slot4])
  label_data[[2]][slot2] <- label_data[[2]][slot2] - 17
  label_data[[2]][slot3 | slot4] <- label_data[[2]][slot3 |slot4] - 109
  # print(label_data)
  res[check.key] <- paste(label_data[[1]],label_data[[2]],label_data[[3]],sep="#")
  res
}
if(F){
  comparison.lattice("selacFULLb","GTR","noneXquadrature",c("v1.6.1-rc1","744c6c8","dd94866"),1,c(3000:3009),min.self.pct = 8)
}
if(F){
  system.time({comparison.lattice("selacFULLb","GTR","noneXquadrature",
                                  c("v1.6.1-rc1","744c6c8","dd94866"),1,
                                  c(3000:3020),min.total.time = 1000) -> test_result_mat_full;})
  dim(test_result_mat3)
  system.time({comparison.lattice("selacFULLb","GTR","noneXquadrature",c("v1.6.1-rc1","744c6c8","dd94866"),1,c(3000:3020),min.total.time = 1000) -> test_result_mat_full;}); dim(test_result_mat3)
  dim(test_result_mat_full)
  save(test_result_mat_full,file="profile_yeast_summary.RData"); apply(test_result_mat_full,1:5,mean,na.rm=T)[,c(1,3),1,,1]
  system.time({comparison.lattice("selacFULLb","GTR","noneXquadrature",c("v1.6.1-rc1","744c6c8","dd94866"),
                                  1,c(3000:3020),min.total.time = 100) -> test_result_mat_full;}); 
  dim(test_result_mat_full); save(test_result_mat_full,file="profile_yeast_summary.RData"); 
  apply(test_result_mat_full,1:5,mean,na.rm=T)[,c(1,3),1,,1] -> test_result_means; 
  test_result_means[test_result_means[,1,1]>700|test_result_means[,1,2]>700|test_result_means[,1,3]>700  ,,]
  sapply(rownames(test_result_means),function(x) any(test_result_means[x,1,] ))
  warnings()
  sapply(rownames(test_result_means),function(x) any(test_result_means[x,1,]>700 ))
  which(sapply(rownames(test_result_means),function(x) any(test_result_means[x,1,]>700 )))
  which(sapply(rownames(test_result_means),function(x) any(test_result_means[x,1,]>700 )))->foo
  test_result_means[foo,1,]
  test_result_means[foo,2,]
  which(sapply(rownames(test_result_means),function(x) any(test_result_means[x,2,]>700 )))->bar
  test_result_means[bar,2,]
  test_result_means[bar,1,]
  test_result_means[bar,2,]
  test_result_means[bar,1,]
  write.csv(test_result_means[bar,1,],file="profile_yeast_21runMeanTime_self.csv")
  which(sapply(rownames(test_result_means),function(x) any(test_result_means[x,2,]>600 )))->bar
  test_result_means[bar,1,]
  write.csv(test_result_means[bar,1,],file="profile_yeast_21runMeanTime_self.csv")
  write.csv(test_result_means[bar,2,],file="profile_yeast_21runMeanTime_total.csv")
  bar
  test_result_means[names(bar),1,]
  all.equal(test_result_means[names(bar),1,],test_result_means[bar,1,])
  apply(test_result_mat_full,1:5,mean,na.rm=T)[,c(5),1,,1] -> test_result_memory;
  test_result_memory[names(bar),]
  write.csv(test_result_memory[names(bar),],file="profile_yeast_21runMeanTime_memory.csv")
  apply(test_result_mat_full,1:5,sd,na.rm=T)[,c(1,3),1,,1] -> test_result_time_sd;
  apply(test_result_mat_full,1:5,sd,na.rm=T)[,c(5),1,,1] -> test_result_memory_sd;
  test_result_time_sd[names(bar),1,]
  test_result_time_sd[names(bar),2,]
  apply(test_result_mat_full,1:5,function(x) sum(is.finite(x)) )[,c(1,5),1,,1] -> test_result_counts;
  all.equal(test_result_counts[,1,],test_result_counts[,2,])
  apply(test_result_mat_full,1:5,function(x) sum(is.finite(x)) )[,1,1,,1] -> test_result_counts;
  test_result_counts
  write.csv(test_result_counts[names(bar),],file="profile_yeast_21run_counts.csv")
  
  write.csv(test_result_counts[,],file="profile_yeast_21run_counts_all.csv")
  write.csv(apply(test_result_mat_full,1:5,mean,na.rm=T)[,,1,,1],file="profile_yeast_21run_means_all.csv")
  write.csv(apply(test_result_mat_full,1:5,sd,na.rm=T)[,,1,,1],file="profile_yeast_21run_stDevs_all.csv")
  
  apply(test_result_counts>10,1,any)
  names(which(apply(test_result_counts>10,1,any)))-> atleast10
  names(which(sapply(rownames(test_result_means),function(x) any(test_result_means[x,2,]>600 ))))->totalTime600
  
  apply(test_result_mat_full,1:5,mean,na.rm=T)[intersect(atleast10,totalTime600),,1,,1]->test_result_imp_means
  # aperm(test_result_imp_means,c(1,3,2))
  write.csv( aperm(test_result_imp_means,c(1,3,2)) ,file="profile_yeast_21run_means_all_alt.csv")
  apply(test_result_mat_full,1:5,sd,na.rm=T)[intersect(atleast10,totalTime600),,1,,1]->test_result_imp_sds
  write.csv( aperm(test_result_imp_sds,c(1,3,2)) ,file="profile_yeast_21run_stDevs_all_alt.csv")
}
if(F){
  system.time({comparison.lattice("ecoliTEST","UNREST","quadrature",
                                  c("v1.6.1-rc1","744c6c8","dd94866"),1,
                                  c(4002),min.total.time = 10) -> test_result_mat_full;})
  system.time({comparison.lattice("ecoliTEST","UNREST","quadrature",
                                  c("v1.6.1-rc1","744c6c8","dd94866"),1,
                                  c(4002),min.total.time = 10) -> test_result_mat_fullb;})
  print(all.equal(test_result_mat_full,test_result_mat_fullb))
  # apply(test_result_mat_full,1:5,function(x) sum(is.finite(x)) )[,1,1,,1] -> test_result_counts;
  # names(which(apply(test_result_counts>10,1,any)))-> atleast10
  print(test_result_mat_full[,1,1,,1,1])
  length(names(which(sapply(rownames(test_result_mat_full),function(x) any(test_result_mat_full[x,3,1,,1,1]>600 ))))->totalTime600)
  tmp_foo<-aperm(test_result_mat_full[totalTime600,,1,,1,1],c(1,3,2))
  print(tmp_foo[order(apply(tmp_foo[,,"total.time"],1,max,na.rm=T),decreasing = T),, ])
}

if(F){
  system.time({comparison.lattice("ecoliSHORT_selac","GTR",
                                  c("none","median",  "quadrature"),
                                  c("v1.6.1-rc1","744c6c8","dd94866","d60d4c6","38a4c36"),1,
                                  c(5100:5129)) -> test_result_mat_full;})
  system.time({comparison.lattice("ecoliSHORT_selac","GTR",
                                  c("none","median",  "quadrature"),
                                  c("v1.6.1-rc1","744c6c8","dd94866","d60d4c6","38a4c36"),1,
                                  c(5100:5129)) -> test_result_mat_fullb;})
  print(all.equal(test_result_mat_full,test_result_mat_fullb))
  
  
  apply(test_result_mat_full,1:5,mean,na.rm=T)[,,,,1] -> test_result_means; 
  apply(test_result_mat_full,1:5,function(x) sum(is.finite(x)) )[,1,,,1] -> test_result_counts;
  names(which(apply(test_result_counts>10,1,any)))-> atleast10
  # dimnames(which(apply(test_result_counts>10,2:3,any)))-> reg.atleast10
  print(test_result_mat_full[,1,1,,1,1])
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>600 ))))->totalTime600)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>300 ))))->totalTime300)
  tmp_foo<-aperm(test_result_means[totalTime600,,,],c(1,4,2,3))
  print(tmp_foo[order(apply(tmp_foo[,,"total.time","ecoliSHORT_selac_GTR_none"],1,max,na.rm=T),decreasing = T),
                ,,"ecoliSHORT_selac_GTR_none" ])
}


if(F){
  system.time({comparison.lattice("ecoliSHORT_selac",c("GTR","UNREST"),
                                  c("none","median",  "quadrature","lognormal"),
                                  c("v1.6.1-rc1","744c6c8","dd94866","38a4c36",
                                    "8f13a5a","9e90977"),1,
                                  c(5100:5129)) -> test_result_mat_full;})
  system.time({comparison.lattice("ecoliSHORT_selac",c("GTR","UNREST"),
                                  c("none","median",  "quadrature","lognormal"),
                                  c("v1.6.1-rc1","744c6c8","dd94866","38a4c36",
                                    "8f13a5a","9e90977"),1,
                                  c(5100:5129)) -> test_result_mat_fullb;})
  print(all.equal(test_result_mat_full,test_result_mat_fullb))
  dimnames(test_result_mat_full)[[4]][-1]<-paste0("r",dimnames(test_result_mat_full)[[4]][-1])
  dim(test_result_mat_full)
  apply(test_result_mat_full,1:5,mean,na.rm=T)[,,,,1] -> test_result_means; 
  apply(test_result_mat_full,1:5,function(x) sum(is.finite(x)) )[,1,,,1] -> test_result_counts;
  length(names(which(apply(test_result_counts>5,1,any)))-> atleast10)
  dim(test_result_means)
  dim(test_result_means[atleast10,,,])
  # dimnames(which(apply(test_result_counts>10,2:3,any)))-> reg.atleast10
  
  print(head(test_result_mat_full[,1,1,,1,1],n=20))
  print(head(test_result_means[,3,1,],n=20))
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>600 ))))->totalTime600)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>300 ))))->totalTime300)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,"total.pct",,,,1]>20 ))))->totalPercent20)
  tmp_foo<-aperm(test_result_means[totalPercent20,,,],c(1,4,2,3))
  print(tmp_foo[order(apply(tmp_foo[,,"total.time","ecoliSHORT_selac_GTR_none"],1,max,na.rm=T),decreasing = T),
                ,,"ecoliSHORT_selac_GTR_none" ])
  tmp_bar <- aperm(test_result_means[totalPercent20,,,],c(1,3,4,2))
  tmp_bar_n <- aperm(test_result_counts[totalPercent20,,],c(1,2,3))
  tmp_bar_n <- tmp_bar_n[order(apply(tmp_bar[,,,"total.time"],1,max,na.rm=T),decreasing = T),, ]
  tmp_bar <- tmp_bar[order(apply(tmp_bar[,,,"total.time"],1,max,na.rm=T),decreasing = T),,, ]
  tmp_bar_names <- dimnames(tmp_bar)
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[1]],tmp_bar_names[[2]])
  tmp_bar_rownames <- paste(tmp_bar_names2[[1]],tmp_bar_names2[[2]],sep="." )
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[3]],tmp_bar_names[[4]])
  tmp_bar_colnames <-paste(tmp_bar_names2[[1]],tmp_bar_names2[[2]],sep="." )
  dim(tmp_bar_n)=c(prod(dim(tmp_bar_n)[1:2]),prod(dim(tmp_bar_n)[-(1:2)]))
  dim(tmp_bar)=c(prod(dim(tmp_bar)[1:2]),prod(dim(tmp_bar)[-(1:2)]))
  dim(tmp_bar)
  rownames(tmp_bar)=tmp_bar_rownames
  rownames(tmp_bar_n)=tmp_bar_rownames
  colnames(tmp_bar)=tmp_bar_colnames
  colnames(tmp_bar_n)=tmp_bar_names[[3]]
  names(dimnames(tmp_bar))<- paste(names(tmp_bar_names)[c(1,3)],names(tmp_bar_names)[c(2,4)],sep=".")
  write.csv(tmp_bar,file="selac_revision_comparison.csv",  na = "")
  write.csv(tmp_bar_n,file="selac_revision_comparison_counts.csv",  na = "")
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[3]],tmp_bar_names[[4]])
  colnames(tmp_bar) <-paste(tmp_bar_names2[[2]],tmp_bar_names2[[1]],sep="." )
  colnames(tmp_bar_n) <-paste("samples",colnames(tmp_bar_n),sep="." )
  write.csv(cbind(tmp_bar,tmp_bar_n),file="selac_revision_comparison_combo.csv",  na = "")
  save.image(file="selac_timings_201807261551.RData")
}



if(F){
  revision.key=c("ab3e84e", "5c98a1f")
  system.time({comparison.lattice(paste0("ecoli",c("FASTHMMTEST","HMMEVAL50"),"_selac"),
                                  c("GTR","UNREST"),
                                  c("none","median",  "quadrature","lognormal"),
                                  revision.key,1,
                                  c(3010:3045)) -> test_result_mat_full;})
  system.time({comparison.lattice(paste0("ecoli",c("FASTHMMTEST","HMMEVAL50"),"_selac"),
                                  c("GTR","UNREST"),
                                  c("none","median",  "quadrature","lognormal"),
                                  revision.key,1,
                                  c(3010:3045)) -> test_result_mat_fullb;})
  print(all.equal(test_result_mat_full,test_result_mat_fullb))
  dimnames(test_result_mat_full)[[4]]<-paste0("r",dimnames(test_result_mat_full)[[4]])
  dim(test_result_mat_full)
  apply(test_result_mat_full,1:5,mean,na.rm=T)[,,,,1] -> test_result_means; 
  apply(test_result_mat_full,1:5,function(x) sum(is.finite(x)) )[,1,,,1] -> test_result_counts;
  length(names(which(apply(test_result_counts>5,1,any)))-> atleast10)
  dim(test_result_means)
  dim(test_result_means[atleast10,,,])
  # dimnames(which(apply(test_result_counts>10,2:3,any)))-> reg.atleast10
  hmm_fast <- paste0("ecoliFASTHMMTEST_selac_UNREST_",c("none","median",  "quadrature","lognormal"))
  hmm_std <- paste0("ecoliHMMEVAL50_selac_UNREST_",c("none","median",  "quadrature","lognormal"))
  keep_versions <- paste0("r",c("ab3e84e", "5c98a1f"))
  length(names(which(apply(test_result_counts[,hmm_fast,]>0,1,any)))-> hmm_fast_lines)
  length(names(which(apply(test_result_counts[,hmm_std ,]>0,1,any)))-> hmm_std_lines)
  
  print(head(test_result_mat_full[hmm_fast_lines,1,hmm_fast[1],keep_versions,1,1],n=20))
  print(head(test_result_mat_full[hmm_std_lines ,1,hmm_std[1], keep_versions,1,1],n=20))
  print(head(test_result_means[hmm_fast_lines,3,hmm_fast[1],keep_versions],n=20))
  print(head(test_result_means[hmm_std_lines, 3,hmm_std[1], keep_versions],n=20))
  print(head(test_result_means[hmm_fast_lines,4,hmm_fast[1],keep_versions],n=20))
  print(head(test_result_means[hmm_std_lines, 4,hmm_std[1], keep_versions],n=20))
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>600 ))))->totalTime600)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>300 ))))->totalTime300)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,"total.pct",,,,1]>20 ))))->totalPercent20)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,"total.pct",,,,1]>1 ))))->totalPercent01)
  tmp_foo<-aperm(test_result_means[intersect(totalPercent20,hmm_fast_lines),,hmm_fast,keep_versions],
                 c(1,4,2,3))
  print(tmp_foo[order(apply(tmp_foo[,,"total.time","ecoliFASTHMMTEST_selac_UNREST_none"],1,max,na.rm=T),decreasing = T),
                ,,"ecoliFASTHMMTEST_selac_UNREST_none" ])
  tmp_foo<-aperm(test_result_means[intersect(totalPercent20,hmm_std_lines),,hmm_std,keep_versions],
                 c(1,4,2,3))
  print(tmp_foo[order(apply(tmp_foo[,,"total.time","ecoliHMMEVAL50_selac_UNREST_none"],1,max,na.rm=T),decreasing = T),
                ,,"ecoliHMMEVAL50_selac_UNREST_none" ])
  tmp_bar <- aperm(test_result_means[totalPercent01,,c(hmm_fast,hmm_std),keep_versions],c(1,3,4,2))
  tmp_bar_n <- aperm(test_result_counts[totalPercent01,c(hmm_fast,hmm_std),keep_versions],c(1,2,3))
  tmp_bar_n <- tmp_bar_n[order(apply(tmp_bar[,,,"total.time"],1,max,na.rm=T),decreasing = T),, ]
  tmp_bar <- tmp_bar[order(apply(tmp_bar[,,,"total.time"],1,max,na.rm=T),decreasing = T),,, ]
  tmp_bar_names <- dimnames(tmp_bar)
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[1]],tmp_bar_names[[2]])
  tmp_bar_rownames <- paste(tmp_bar_names2[[1]],tmp_bar_names2[[2]],sep="." )
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[3]],tmp_bar_names[[4]])
  tmp_bar_colnames <-paste(tmp_bar_names2[[1]],tmp_bar_names2[[2]],sep="." )
  dim(tmp_bar_n)=c(prod(dim(tmp_bar_n)[1:2]),prod(dim(tmp_bar_n)[-(1:2)]))
  dim(tmp_bar)=c(prod(dim(tmp_bar)[1:2]),prod(dim(tmp_bar)[-(1:2)]))
  dim(tmp_bar)
  rownames(tmp_bar)=tmp_bar_rownames
  rownames(tmp_bar_n)=tmp_bar_rownames
  colnames(tmp_bar)=tmp_bar_colnames
  colnames(tmp_bar_n)=tmp_bar_names[[3]]
  names(dimnames(tmp_bar))<- paste(names(tmp_bar_names)[c(1,3)],names(tmp_bar_names)[c(2,4)],sep=".")
  write.csv(tmp_bar,file="selac_revision_hmm_comparison.csv",  na = "")
  write.csv(tmp_bar_n,file="selac_revision_hmm_comparison_counts.csv",  na = "")
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[3]],tmp_bar_names[[4]])
  colnames(tmp_bar) <-paste(tmp_bar_names2[[2]],tmp_bar_names2[[1]],sep="." )
  colnames(tmp_bar_n) <-paste("samples",colnames(tmp_bar_n),sep="." )
  write.csv(cbind(tmp_bar,tmp_bar_n),file="selac_revision_hmm_comparison_combo.csv",  na = "")
  save.image(file=format(Sys.time(),format = "selac_timings_hmm_%Y%m%d%H%M.RData"))
}


if(F){
  revision.key=c("ab3e84e", "5c98a1f")
  system.time({comparison.lattice(paste0("ecoli",c("FASTHMMTEST","HMMEVAL50"),"_selac"),
                                  c("GTR","UNREST"),
                                  c("none","median",  "quadrature","lognormal"),
                                  revision.key,1,
                                  c(3010:3045)) -> test_result_mat_full;})
  system.time({comparison.lattice(paste0("ecoli",c("FASTHMMTEST","HMMEVAL50"),"_selac"),
                                  c("GTR","UNREST"),
                                  c("none","median",  "quadrature","lognormal"),
                                  revision.key,1,
                                  c(3010:3045)) -> test_result_mat_fullb;})
  print(all.equal(test_result_mat_full,test_result_mat_fullb))
  dimnames(test_result_mat_full)[[4]]<-paste0("r",dimnames(test_result_mat_full)[[4]])
  dim(test_result_mat_full)
  apply(test_result_mat_full,1:5,mean,na.rm=T)[,,,,1] -> test_result_means; 
  apply(test_result_mat_full,1:5,function(x) sum(is.finite(x)) )[,1,,,1] -> test_result_counts;
  length(names(which(apply(test_result_counts>5,1,any)))-> atleast10)
  dim(test_result_means)
  dim(test_result_means[atleast10,,,])
  # dimnames(which(apply(test_result_counts>10,2:3,any)))-> reg.atleast10
  hmm_fast <- paste0("ecoliFASTHMMTEST_selac_UNREST_",c("none","median",  "quadrature","lognormal"))
  hmm_std <- paste0("ecoliHMMEVAL50_selac_UNREST_",c("none","median",  "quadrature","lognormal"))
  keep_versions <- paste0("r",c("ab3e84e", "5c98a1f"))
  length(names(which(apply(test_result_counts[,hmm_fast,]>0,1,any)))-> hmm_fast_lines)
  length(names(which(apply(test_result_counts[,hmm_std ,]>0,1,any)))-> hmm_std_lines)
  
  print(head(test_result_mat_full[hmm_fast_lines,1,hmm_fast[2],keep_versions,1,1],n=20))
  print(head(test_result_mat_full[hmm_std_lines ,1,hmm_std[2], keep_versions,1,1],n=20))
  print(head(test_result_means[hmm_fast_lines,3,hmm_fast[2],keep_versions],n=20))
  print(head(test_result_means[hmm_std_lines, 3,hmm_std[2], keep_versions],n=20))
  print(head(test_result_means[hmm_fast_lines,4,hmm_fast[2],keep_versions],n=20))
  print(head(test_result_means[hmm_std_lines, 4,hmm_std[2], keep_versions],n=20))
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>600 ))))->totalTime600)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,3,,,,1]>300 ))))->totalTime300)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,"total.pct",,,,1]>20 ))))->totalPercent20)
  length(names(which(sapply(rownames(test_result_mat_full),
                            function(x) any(test_result_mat_full[x,"total.pct",,,,1]>1 ))))->totalPercent01)
  tmp_foo<-aperm(test_result_means[intersect(totalPercent20,hmm_fast_lines),,hmm_fast,keep_versions],
                 c(1,4,2,3))
  print(tmp_foo[order(apply(tmp_foo[,,"total.time","ecoliFASTHMMTEST_selac_UNREST_lognormal"],1,max,na.rm=T),decreasing = T),
                ,,"ecoliFASTHMMTEST_selac_UNREST_lognormal" ])
  tmp_foo<-aperm(test_result_means[intersect(totalPercent20,hmm_std_lines),,hmm_std,keep_versions],
                 c(1,4,2,3))
  print(tmp_foo[order(apply(tmp_foo[,,"total.time","ecoliHMMEVAL50_selac_UNREST_lognormal"],1,max,na.rm=T),decreasing = T),
                ,,"ecoliHMMEVAL50_selac_UNREST_lognormal" ])
  tmp_bar <- aperm(test_result_means[totalPercent01,,c(hmm_fast,hmm_std),keep_versions],c(1,3,4,2))
  tmp_bar_n <- aperm(test_result_counts[totalPercent01,c(hmm_fast,hmm_std),keep_versions],c(1,2,3))
  tmp_bar_n <- tmp_bar_n[order(apply(tmp_bar[,,,"total.time"],1,max,na.rm=T),decreasing = T),, ]
  tmp_bar <- tmp_bar[order(apply(tmp_bar[,,,"total.time"],1,max,na.rm=T),decreasing = T),,, ]
  tmp_bar_names <- dimnames(tmp_bar)
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[1]],tmp_bar_names[[2]])
  tmp_bar_rownames <- paste(tmp_bar_names2[[1]],tmp_bar_names2[[2]],sep="." )
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[3]],tmp_bar_names[[4]])
  tmp_bar_colnames <-paste(tmp_bar_names2[[1]],tmp_bar_names2[[2]],sep="." )
  dim(tmp_bar_n)=c(prod(dim(tmp_bar_n)[1:2]),prod(dim(tmp_bar_n)[-(1:2)]))
  dim(tmp_bar)=c(prod(dim(tmp_bar)[1:2]),prod(dim(tmp_bar)[-(1:2)]))
  dim(tmp_bar)
  rownames(tmp_bar)=tmp_bar_rownames
  rownames(tmp_bar_n)=tmp_bar_rownames
  colnames(tmp_bar)=tmp_bar_colnames
  colnames(tmp_bar_n)=tmp_bar_names[[3]]
  names(dimnames(tmp_bar))<- paste(names(tmp_bar_names)[c(1,3)],names(tmp_bar_names)[c(2,4)],sep=".")
  write.csv(tmp_bar,file="selac_revision_hmm_scan_comparison.csv",  na = "")
  write.csv(tmp_bar_n,file="selac_revision_hmm_scan_comparison_counts.csv",  na = "")
  tmp_bar_names2 <- expand.grid(tmp_bar_names[[3]],tmp_bar_names[[4]])
  colnames(tmp_bar) <-paste(tmp_bar_names2[[2]],tmp_bar_names2[[1]],sep="." )
  colnames(tmp_bar_n) <-paste("samples",colnames(tmp_bar_n),sep="." )
  write.csv(cbind(tmp_bar,tmp_bar_n),file="selac_revision_hmm_scan_comparison_combo.csv",  na = "")
  
  str(site_pattern_counts <- read.csv("site_pattern_counts.csv",row.names = 1))
  test_results_table <- na.omit(as.data.frame.table(test_result_mat_full,responseName = "Value"))

  # xyplot(log10(Freq)~factor(seed), groups = factor(src,levels = hmm_fast),
  #        data=as.data.frame.table(tapply(test_results_table$Value,test_results_table[c(2,3,4,6)],max),stringsAsFactors = F),
  #        subset=measure=="total.time",auto.key=list(columns=2))
  test_results_max<-as.data.frame.table(tapply(test_results_table$Value,
                                               test_results_table[c(2,3,4,6)],max, na.rm=T),
                                        responseName = "MAX",
                                        stringsAsFactors = F)
  test_results_max$mode=paste0("ecoli",c("FASTHMMTEST","HMMEVAL50","HMMEVALFULL"))
  test_results_max$mode[grepl("^ecoliFASTHMMTEST",test_results_max$src)] = "ecoliFASTHMMTEST"
  test_results_max$mode[grepl("^ecoliHMMEVAL50",test_results_max$src)] = "ecoliHMMEVAL50"
  test_results_max$mode[grepl("^ecoliHMMEVALFULL",test_results_max$src)] = "ecoliHMMEVALFULL"
  test_results_max$gamma=""
  test_results_max$gamma[grepl("none$",test_results_max$src)] = "none"
  test_results_max$gamma[grepl("median$",test_results_max$src)] = "median"
  test_results_max$gamma[grepl("quadrature$",test_results_max$src)] = "quadrature"
  test_results_max$gamma[grepl("lognormal$",test_results_max$src)] = "lognormal"
  site_pattern_counts$seed.key = paste0("S",site_pattern_counts$seed )
  
  summary(site_pattern_counts_uniq <- unique(site_pattern_counts[c(1,4,10,8,9)]))
  row.names(site_pattern_counts_uniq) <- with(site_pattern_counts_uniq,paste(block, gamma,seed.key,sep="_"))
  str(site_pattern_counts_uniq[with(test_results_max,paste(mode, gamma,seed,sep="_"))  ,4:5])
  test_results_max_iter <- cbind(test_results_max, site_pattern_counts_uniq[with(test_results_max,paste(mode, gamma,seed,sep="_"))  ,4:5])
  
  xyplot(MAX~neval,groups=ver,subset = measure=="total.time",data=test_results_max_iter, 
         auto.key=list(columns=2))  
  xyplot(MAX~jitter(neval),groups=ver,subset = measure=="total.time",data=test_results_max_iter,
         auto.key=list(columns=2))
  max_row_common <- test_results_max_iter[test_results_max_iter$measure=="mem.total",c(2:4,6:9)]
  max_row_keys <- with(max_row_common,paste(src, ver,seed,sep="_"))
  max_row_keys <- sub(pattern="^ecoliHMMEVAL50_selac_", replacement = "eH50_",x = max_row_keys)
  max_row_keys <- sub(pattern="^ecoliHMMEVALFULL_selac_", replacement = "eHMM_",x = max_row_keys)
  max_row_keys <- sub(pattern="^ecoliFASTHMMTEST_selac_", replacement = "eH10_",x = max_row_keys)
  max_row_keys <- sub(pattern="_UNREST", replacement = "U",x = max_row_keys)
  max_row_keys <- sub(pattern="_GTR", replacement = "G",x = max_row_keys)
  max_row_keys <- sub(pattern="_none", replacement = "n",x = max_row_keys)
  max_row_keys <- sub(pattern="_quadrature", replacement = "q",x = max_row_keys)
  max_row_keys <- sub(pattern="_median", replacement = "m",x = max_row_keys)
  max_row_keys <- sub(pattern="_lognormal", replacement = "l",x = max_row_keys)
  max_row_keys <- sub(pattern="_r([0-9a-f]{3})[0-9a-f]{4}", replacement = "\\1",x = max_row_keys)
  max_row_keys <- sub(pattern="_S", replacement = "S",x = max_row_keys)
  rownames(max_row_common) <- max_row_keys
  test_results_max_row <- 
    foreach(key=max_row_keys,
            data_block=iter(as.matrix(test_results_max_iter$MAX),by="row",chunksize=5),
            measure_block=iter(as.matrix(test_results_max_iter$measure),by="row",chunksize=5),
            .combine=rbind) %:%
    foreach(data=iter(data_block,by="row",chunksize=1),
            measure=iter(measure_block,by="row",chunksize=1),
            .combine=cbind) %do% matrix(data=data,dimnames = list(key,measure))
    # {
    #           cat("+")
    #           res<-data.frame(X=data,row.names = key); colnames(res)<-names(res) <- measure; res
    #         }
  test_results_max_row <- cbind(max_row_common,test_results_max_row)
  xyplot(total.time~jitter(neval),groups=ver,data=test_results_max_row,
         auto.key=list(columns=2))
  
  xyplot(mem.total~jitter(neval),groups=ver,data=test_results_max_row,
         auto.key=list(columns=2))
  save.image(file=format(Sys.time(),format = "selac_timings_hmm_scan_%Y%m%d%H%M.RData"))
}
