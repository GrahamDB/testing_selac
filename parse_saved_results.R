

validation_table_row <- function(mode=c("SHORTTEST","TEST","SHORT"),
                                 seed=sample.int(1e6,1),
                                 codon.model=c("selac","none","GY94","YN98"),
                                 nuc.model=c("GTR","UNREST","JC"),
                                 selac_release="v1.6.1-rc1", 
                                 include.gamma=T,
                                 gamma.type=c("quadrature","median","lognormal","none"),
                                 nCores=1){
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
      if(!is.null(result$loglik)) {
        
        return(data.frame(SRC=src.key,
                          Nuc.Model=nuc.model,
                          Gamma.model=gamma.mode,
                          Revision=selac_release,
                          nCores=as.integer(nCores),
                          seed=as.integer(seed),
                          model.LL=result$loglik,
                          stringsAsFactors = F))
      } 
    })
    cat(sprintf("Rebuilding: %s\n",
                profile_prefix))
  }
  return(data.frame(SRC=character(),
                    Nuc.Model=character(),
                    Gamma.model=character(),
                    Revision=character(),
                    nCores=integer(),
                    seed=integer(),
                    model.LL=numeric(),
                    stringsAsFactors = F))
}

scan_for_completed_files <- function(mode=c("SHORTTEST","TEST","SHORT"),
                                     seed=5100:5129,
                                     codon.model=c("selac","none","GY94","YN98"),
                                     nuc.model=c("GTR","UNREST","JC"),
                                     selac_release=c("v1.6.1-rc1","744c6c8","38a4c36","d60d4c6","dd94866" ),
                                     gamma.mode=c("quadrature","median","lognormal","none"),
                                     nCores=1){
  mode=match.arg(mode)
  src.key=paste0("ecoli",mode)
  codon.model = match.arg(codon.model)
  profile_regex=sprintf("^(%s)_(%s)_(%s)_(%s)_(%s)_(%s)_(%s)_result.Rdata",
                        paste0(src.key,collapse = "|"),
                        paste0(codon.model,collapse = "|"),
                         paste0(nuc.model,collapse="|"),
                         paste0(gamma.mode,collapse="|"),
                         paste0(selac_release,collapse="|"),
                         paste0(nCores,collapse="|"),
                         paste0(seed,collapse="|"))
  print(profile_regex)
  junk=dir(path="ecoli_output/",pattern = profile_regex
        )
  tmp=t(sapply(strsplit(x=junk,split = "_"),as.character))[,3:7]
  colnames(tmp)<- c("nuc.model","gamma.type","selac_release","nCores","seed")
  tmp<-as.data.frame(tmp,stringsAsFactors = F)
  tmp<-cbind(mode=mode,codon.model=codon.model,tmp,stringsAsFactors = F)
  tmp$nCores=as.integer(tmp$nCores)
  tmp$seed=as.integer(tmp$seed)
  return(tmp)
}
build_validation_table <- function(scan_keys){
  res <- do.call(rbind,do.call(mapply,c(list(FUN=validation_table_row),as.list(scan_keys),list(USE.NAMES=F,SIMPLIFY=F))))
  cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n")
  # print(colnames(res))
  # str(res)
  with(res,cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
                       SRC,
                       Nuc.Model,
                       Gamma.model,
                       Revision,
                       nCores,
                       seed,
                       model.LL),sep=""))
  return(invisible(res))
}

rebuild_validation_table <- function(mode=c("SHORTTEST","TEST","SHORT"),
                                     codon.model=c("selac","none","GY94","YN98")){
  mode=match.arg(mode)
  src.key=paste0("ecoli",mode)
  codon.model = match.arg(codon.model)
  scan_keys<-scan_for_completed_files(mode=mode,codon.model = codon.model ,seed="[0-9]+")
  res <- do.call(rbind,do.call(mapply,c(list(FUN=validation_table_row),as.list(scan_keys),list(USE.NAMES=F,SIMPLIFY=F))))
  cat("SRC,Nuc.Model,Gamma.model,Revision,nCores,seed,model.LL\n",
      file=paste0(src.key,"_",codon.model,"_LL_log.csv"),
      append = F)
  # print(colnames(res))
  # str(res)
  with(res,cat(sprintf("\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%0.3f\n",
                       SRC,
                       Nuc.Model,
                       Gamma.model,
                       Revision,
                       nCores,
                       seed,
                       model.LL),sep="",
               file=paste0(src.key,"_",codon.model,"_LL_log.csv"),
               append = T))
  return(invisible(res))
}


write_validation_table_direct <- function(x){
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
}

load_ecoli_profile_mode <- function(mode=c("SHORTTEST","TEST","SHORT"),
                                   seed=sample.int(1e6,1),
                                   codon.model=c("selac","none","GY94","YN98"),
                                   nuc.model=c("GTR","UNREST","JC"),
                                   ref="v1.6.1-rc1", 
                                   include.gamma=T,
                                   gamma.type=c("quadrature","median","lognormal","none"),
                                   nCores=1){
  # setup_selac_for_profiling(ref=ref)
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
}