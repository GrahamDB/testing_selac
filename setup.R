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
setup_selac_for_profiling <- function(ref="v1.6.1-rc1"){
  if(!is.null(selac_release)){
    ref=selac_release
  } else {
    selac_release <<- ref
  }
  local_path<-file.path(getwd(),
                        R.Version()$platform,
                        paste0(R.Version()$major,
                               ".",
                               strsplit(R.Version()$minor,".",fixed = T)[[1]][1]),
                        paste0("selac_",ref));
  
  if(!dir.exists(local_path))
    dir.create(local_path,recursive = T,mode = "0755")
  
  .libPaths(local_path)
  
  if(!require(selac)){
    install.packages(c("expm", "MASS", "parallel", "phangorn", "seqinr", "statmod", 
      "zoo", "RColorBrewer","ape","deSolve","nloptr","deSolve","nnet"),type="source", 
      INSTALL_opts="--with-keep.source")
    install_github("GrahamDB/selac",ref=ref)
    if(!require(selac))
      stop("Failed to install selac")
  }
}
