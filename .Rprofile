# working directory init

.First <-function(){
  local_path<-file.path(getwd(),
                        R.Version()$platform,
                        paste0(R.Version()$major,
                               ".",
                               strsplit(R.Version()$minor,".",fixed = T)[[1]][1]))
  if(!dir.exists(local_path))
    dir.create(local_path,recursive = T,mode = "0755")
  .libPaths(local_path)
  options(Ncpus=12)
}