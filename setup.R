# setup
source(".Rprofile")
.First()
dep_imports<-strsplit("expm, MASS, parallel, phangorn, seqinr, statmod, zoo, RColorBrewer",", ")[[1]]
install.packages(c("ape","deSolve","nloptr","deSolve","nnet"), 
                 type="source", INSTALL_opts="--with-keep.source")
install.packages(path.expand(file.path("~","git","selac_1.6.1.tar.gz")) , repos=NULL,
                 type="source", INSTALL_opts="--with-keep.source")

install.packages(dep_imports,type="source", INSTALL_opts="--with-keep.source")
