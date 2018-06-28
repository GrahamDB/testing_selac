# setup
# source(".Rprofile")
# .First()
# dep_imports<-strsplit("expm, MASS, parallel, phangorn, seqinr, statmod, zoo, RColorBrewer",", ")[[1]]
# install.packages(c("ape","deSolve","nloptr","deSolve","nnet"), 
#                  type="source", INSTALL_opts="--with-keep.source")
# install.packages(path.expand(file.path("~","git","selac_1.6.1.tar.gz")) , repos=NULL,
#                  type="source", INSTALL_opts="--with-keep.source")
# 
# install.packages(dep_imports,type="source", INSTALL_opts="--with-keep.source")
if(length(ls(pattern="^selac_release$"))==0)
  selac_release=NULL
if(length(ls(pattern="^user_path$"))==0)
  user_path=.libPaths()
print(user_path)
library(parallel,lib.loc = user_path)
library(devtools,lib.loc = user_path)
setup_selac_for_profiling <- function(ref="v1.6.1-rc1"){
  if(!is.null(selac_release)){
    ref=selac_release
  } else {
    selac_release <<- ref
  }
  print(selac_release)
  library(parallel)
  local_path<-file.path(getwd(),
                        R.Version()$platform,
                        paste0(R.Version()$major,
                               ".",
                               strsplit(R.Version()$minor,".",fixed = T)[[1]][1]),"common");
  ref_path <- file.path(getwd(),
                        R.Version()$platform,
                        paste0(R.Version()$major,
                               ".",
                               strsplit(R.Version()$minor,".",fixed = T)[[1]][1]),
                        paste0("selac_",ref));
  
  if(!dir.exists(ref_path))
    dir.create(ref_path,recursive = T,mode = "0755")
  if(!dir.exists(local_path))
    dir.create(local_path,recursive = T,mode = "0755")
  
  .libPaths(c(local_path,ref_path))
  
  if(! "selac" %in% installed.packages()[,"Package"]){
    cat("\nSelac not found, installing.\n\n")
    library(parallel,lib.loc = user_path)
    library(devtools,lib.loc = user_path)
    library(httr,lib.loc = user_path)
    library(curl,lib.loc = user_path)
    need_packs=c("expm", "MASS", "phangorn", "seqinr", "statmod", 
                 "zoo", "RColorBrewer","ape","deSolve","nloptr","deSolve","nnet")
    if(any(!need_packs %in% installed.packages()[,"Package"] )){
      
      .libPaths(c(local_path))
      install.packages(c("expm", "MASS", "phangorn", "seqinr", "statmod", 
                         "zoo", "RColorBrewer","ape","deSolve","nloptr","deSolve","nnet"
      )[!need_packs %in% installed.packages()[,"Package"]],
      type="source", INSTALL_opts="--with-keep.source")
    }
    .libPaths(c(ref_path,local_path))
    install_github("GrahamDB/selac",ref=ref, args="--with-keep.source")
    if(!"selac" %in% installed.packages()[,"Package"])
      stop("Failed to install selac")
  }
  library(selac)
}
