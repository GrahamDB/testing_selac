#profvis json parsing functions

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


