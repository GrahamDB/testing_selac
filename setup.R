# setup
source(".Rprofile")
.First()

install.packages(c("ape","deSolve","nloptr","deSolve") )
install.packages(c("nnet") )
install.packages("/home/asterion/git/selac_1.6.1.tar.gz" , repos=NULL)
strsplit("expm, MASS, parallel, phangorn, seqinr, statmod, zoo, RColorBrewer",", ")[[1]][!strsplit("expm, MASS, parallel, phangorn, seqinr, statmod, zoo, RColorBrewer",", ")[[1]] %in% rownames(installed.packages())]-> temp
cat(deparse(temp),sep="\n")
install.packages(c("expm", "phangorn", "statmod", "zoo", "RColorBrewer"))
