library(expm)
library(microbenchmark)


test_short_exmpAtv <- function(n.mat=3,n.branches=4) {
  mat=matrix(runif(n=n.mat*n.mat),nrow=n.mat); 
  #system.time({ res.1 = ex})
  vec.branches=lapply(seq_len(n.branches),function(x) matrix(runif(n=n.mat),nrow=n.mat))
  t.branches=lapply(seq_len(n.branches),function(x) runif(n=1))
  m2.times=system.time({
    m2.res<-lapply(seq_len(n.branches), 
                   function(x) expAtv(mat,vec.branches[[x]],t.branches[[x]] )$eAtv )})
  m1.times=system.time({
    m1.res<-lapply(seq_len(n.branches), 
                   function(x) (expm(mat*t.branches[[x]]) %*% vec.branches[[x]] )[,1]) })
  #str(m1.res)
  all.times=rbind(m1=m1.times,m2=m2.times)
  list(all.equal(m1.res,m2.res),all.times,apply(all.times,MARGIN = 2,diff),log(all.times))
}


test_codon_exmpAtv <- function(n.mat=3,n.branches=4, n.codon=100) {
  mat=matrix(runif(n=n.mat*n.mat),nrow=n.mat); 
  #system.time({ res.1 = ex})
  vec.branches=lapply(seq_len(n.codon),
                      function(y)
                        lapply(seq_len(n.branches),
                      function(x) {
                        res<-matrix(runif(n=n.mat),nrow=n.mat);
  res/sum(res)}))
  t.branches=lapply(seq_len(n.branches),function(x) runif(n=1))
  m2.times=system.time({
    m2.res<-lapply(seq_len(n.codon),
                   function(y) lapply(seq_len(n.branches), 
                   function(x) expAtv(mat,vec.branches[[y]][[x]],t.branches[[x]] )$eAtv ))})
  m1.times=system.time({
    m1.mat=lapply(seq_len(n.branches), 
                  function(x) (expm(mat*t.branches[[x]])))
    m1.res<-lapply(seq_len(n.codon),
                   function(y) lapply(seq_len(n.branches), 
                   function(x) (m1.mat[[x]] %*% vec.branches[[y]][[x]] )[,1]) )})
  #str(m1.res)
  all.times=rbind(m1=m1.times,m2=m2.times)
  list(all.equal(m1.res,m2.res),all.times,apply(all.times,MARGIN = 2,diff),log(all.times))
}

test_codon_exmpAtv_W77 <- function(n.mat=3,n.branches=4, n.codon=100) {
  mat=matrix(runif(n=n.mat*n.mat),nrow=n.mat); 
  #system.time({ res.1 = ex})
  vec.branches=lapply(seq_len(n.codon),
                      function(y)
                        lapply(seq_len(n.branches),
                               function(x) {
                                 res<-matrix(runif(n=n.mat),nrow=n.mat);
                                 res/sum(res)}))
  t.branches=lapply(seq_len(n.branches),function(x) runif(n=1))
  m2.times=system.time({
    m2.res<-lapply(seq_len(n.codon),
                   function(y) lapply(seq_len(n.branches), 
                                      function(x) expAtv(mat,vec.branches[[y]][[x]],t.branches[[x]] )$eAtv ))})
  m1.times=system.time({
    m1.mat=lapply(seq_len(n.branches), 
                  function(x) (expm(mat*t.branches[[x]],method = "Ward77")))
    m1.res<-lapply(seq_len(n.codon),
                   function(y) lapply(seq_len(n.branches), 
                                      function(x) (m1.mat[[x]] %*% vec.branches[[y]][[x]] )[,1]) )})
  #str(m1.res)
  all.times=rbind(m1=m1.times,m2=m2.times)
  list(all.equal(m1.res,m2.res),all.times,apply(all.times,MARGIN = 2,diff),log(all.times))
}

test_codon_exmpAtv_alt <- function(n.mat=3,n.branches=4, n.codon=100) {
  mat=matrix(runif(n=n.mat*n.mat),nrow=n.mat); 
  #system.time({ res.1 = ex})
  vec.branches=lapply(seq_len(n.codon),
                      function(y)
                        lapply(seq_len(n.branches),
                               function(x) {
                                 res<-matrix(runif(n=n.mat),nrow=n.mat);
                                 res/sum(res)}))
  t.branches=lapply(seq_len(n.branches),function(x) runif(n=1))
  m2.times=system.time({
    m2.res<-lapply(seq_len(n.branches),
                   function(x) lapply(seq_len(n.codon), 
                                      function(y) expAtv(mat,vec.branches[[y]][[x]],t.branches[[x]] )$eAtv ))})
  m1.times=system.time({
    m1.res<-lapply(seq_len(n.branches),
                   function(x) {tmp.mat=expm(mat*t.branches[[x]]);
                   lapply(seq_len(n.codon), 
                                      function(y) (tmp.mat %*% vec.branches[[y]][[x]] )[,1]);} )})
  #str(m1.res)
  all.times=rbind(m1=m1.times,m2=m2.times)
  list(all.equal(m1.res,m2.res),all.times,apply(all.times,MARGIN = 2,diff),log(all.times))
}